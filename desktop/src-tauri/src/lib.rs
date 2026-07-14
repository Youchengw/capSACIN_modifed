//! capSACIN Studio — Tauri 2 backend.
//!
//! Owns the bundled Python sidecar, correlates JSON-line responses with the
//! request that produced them, and exposes file/dialog helpers to React.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

use tauri::{AppHandle, Emitter, Manager};
use tauri_plugin_shell::process::{CommandChild, CommandEvent};
use tauri_plugin_shell::ShellExt;
use tokio::sync::oneshot;

type PendingResult = Result<serde_json::Value, String>;

struct SidecarState {
    child: Mutex<Option<CommandChild>>,
    pending: Mutex<HashMap<String, oneshot::Sender<PendingResult>>>,
    expected_termination: Mutex<bool>,
}

static NEXT_REQ_ID: AtomicU64 = AtomicU64::new(1);

fn fail_all_pending(app: &AppHandle, message: &str) {
    let state = app.state::<SidecarState>();
    let pending = std::mem::take(&mut *state.pending.lock().unwrap());
    for (_, sender) in pending {
        let _ = sender.send(Err(message.to_string()));
    }
}

fn start_sidecar(app: &AppHandle) -> Result<CommandChild, String> {
    let sidecar_cmd = app
        .shell()
        .sidecar("capsacin-sidecar")
        .map_err(|e| format!("Sidecar not found: {e}"))?;

    let (mut rx, child) = sidecar_cmd
        .spawn()
        .map_err(|e| format!("Failed to spawn sidecar: {e}"))?;

    let app_handle = app.clone();
    std::thread::spawn(move || {
        while let Some(event) = rx.blocking_recv() {
            match event {
                CommandEvent::Stdout(line) => {
                    let text = String::from_utf8_lossy(&line);
                    for raw_line in text.lines() {
                        let trimmed = raw_line.trim();
                        if trimmed.is_empty() {
                            continue;
                        }
                        let Ok(parsed) = serde_json::from_str::<serde_json::Value>(trimmed) else {
                            continue;
                        };

                        match parsed["type"].as_str().unwrap_or("") {
                            "progress" => {
                                let _ = app_handle.emit("sidecar:progress", &parsed);
                            }
                            "result" => {
                                let _ = app_handle.emit("sidecar:result", &parsed);
                                if let Some(request_id) = parsed["id"].as_str() {
                                    let sender = app_handle
                                        .state::<SidecarState>()
                                        .pending
                                        .lock()
                                        .unwrap()
                                        .remove(request_id);
                                    if let Some(sender) = sender {
                                        let status = parsed["status"].as_str().unwrap_or("error");
                                        let result = if status == "ok" {
                                            Ok(parsed
                                                .get("data")
                                                .cloned()
                                                .unwrap_or(serde_json::Value::Null))
                                        } else if status == "cancelled" {
                                            Err("Operation cancelled".to_string())
                                        } else {
                                            Err(parsed["error"]["message"]
                                                .as_str()
                                                .unwrap_or("Unknown sidecar error")
                                                .to_string())
                                        };
                                        let _ = sender.send(result);
                                    }
                                }
                            }
                            _ => {}
                        }
                    }
                }
                CommandEvent::Stderr(line) => {
                    eprintln!("[sidecar] {}", String::from_utf8_lossy(&line));
                }
                CommandEvent::Error(error) => {
                    fail_all_pending(&app_handle, &format!("Sidecar error: {error}"));
                }
                CommandEvent::Terminated(status) => {
                    let expected = {
                        let state = app_handle.state::<SidecarState>();
                        let mut expected = state.expected_termination.lock().unwrap();
                        std::mem::replace(&mut *expected, false)
                    };
                    if !expected {
                        fail_all_pending(&app_handle, &termination_message(&status));
                    }
                    *app_handle.state::<SidecarState>().child.lock().unwrap() = None;
                    break;
                }
                _ => {}
            }
        }
    });

    Ok(child)
}

fn termination_message(status: &tauri_plugin_shell::process::TerminatedPayload) -> String {
    if status.signal == Some(9) {
        "The analysis engine was stopped by macOS (signal 9), usually because of memory pressure. Close memory-heavy apps or try a smaller structure, then retry; the engine will restart automatically.".to_string()
    } else {
        format!("Analysis engine terminated unexpectedly: {status:?}")
    }
}

async fn call_sidecar(
    app: &AppHandle,
    operation: &str,
    params: serde_json::Value,
) -> PendingResult {
    let request_id = format!("req-{}", NEXT_REQ_ID.fetch_add(1, Ordering::Relaxed));
    let request = serde_json::json!({
        "id": request_id,
        "operation": operation,
        "params": params,
    });
    let request_line = request.to_string() + "\n";
    let (sender, receiver) = oneshot::channel();

    {
        let state = app.state::<SidecarState>();
        state
            .pending
            .lock()
            .unwrap()
            .insert(request_id.clone(), sender);

        let mut child_guard = state.child.lock().unwrap();
        if child_guard.is_none() {
            match start_sidecar(app) {
                Ok(child) => *child_guard = Some(child),
                Err(error) => {
                    drop(child_guard);
                    state.pending.lock().unwrap().remove(&request_id);
                    return Err(error);
                }
            }
        }

        if let Some(child) = child_guard.as_mut() {
            if let Err(error) = child.write(request_line.as_bytes()) {
                state.pending.lock().unwrap().remove(&request_id);
                return Err(format!("Write to sidecar failed: {error}"));
            }
        }
    }

    match tokio::time::timeout(std::time::Duration::from_secs(1800), receiver).await {
        Ok(Ok(result)) => result,
        Ok(Err(_)) => Err("Sidecar response channel closed".to_string()),
        Err(_) => {
            app.state::<SidecarState>()
                .pending
                .lock()
                .unwrap()
                .remove(&request_id);
            Err("Operation timed out after 30 minutes".to_string())
        }
    }
}

fn kill_sidecar(app: &AppHandle) -> Result<(), String> {
    let state = app.state::<SidecarState>();
    // Mark the exit before sending SIGKILL. Otherwise the termination event
    // can race ahead and turn an intentional Cancel into an unexpected-error
    // message for the pending request.
    *state.expected_termination.lock().unwrap() = true;
    let child = state.child.lock().unwrap().take();
    if child.is_none() {
        *state.expected_termination.lock().unwrap() = false;
    }

    // Resolve callers as cancelled before killing the process so the UI never
    // receives the raw signal-9 termination event for a user cancellation.
    fail_all_pending(app, "Operation cancelled");
    child
        .map(|child| {
            child
                .kill()
                .map_err(|error| format!("Failed to stop sidecar: {error}"))
        })
        .unwrap_or(Ok(()))
}

#[tauri::command]
async fn inspect_structure(app: AppHandle, input_path: String) -> PendingResult {
    call_sidecar(
        &app,
        "inspect_structure",
        serde_json::json!({ "input_path": input_path }),
    )
    .await
}

#[tauri::command]
async fn prepare_preview(app: AppHandle, params: serde_json::Value) -> PendingResult {
    let mut values = params.as_object().cloned().unwrap_or_default();
    values.remove("weight");
    call_sidecar(&app, "prepare_preview", serde_json::Value::Object(values)).await
}

#[tauri::command]
async fn run_slice(app: AppHandle, params: serde_json::Value) -> PendingResult {
    call_sidecar(&app, "run_slice", params).await
}

#[tauri::command]
async fn cancel_operation(app: AppHandle) -> Result<(), String> {
    kill_sidecar(&app)
}

#[tauri::command]
async fn open_pdb_dialog(app: AppHandle) -> Result<Option<String>, String> {
    use tauri_plugin_dialog::DialogExt;
    let path = app
        .dialog()
        .file()
        .add_filter("PDB Files", &["pdb", "ent"])
        .blocking_pick_file();
    Ok(path.map(|path| path.to_string()))
}

#[tauri::command]
async fn save_pdb_dialog(app: AppHandle, source_path: String) -> Result<Option<String>, String> {
    use tauri_plugin_dialog::DialogExt;
    let path = app
        .dialog()
        .file()
        .add_filter("PDB Files", &["pdb"])
        .blocking_save_file();
    if let Some(destination) = path {
        std::fs::copy(&source_path, destination.to_string())
            .map_err(|error| format!("Failed to copy PDB: {error}"))?;
        Ok(Some(destination.to_string()))
    } else {
        Ok(None)
    }
}

#[tauri::command]
fn read_text_file(path: String) -> Result<String, String> {
    let path = Path::new(&path);
    let extension = path
        .extension()
        .and_then(|value| value.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    if extension != "cif" && extension != "mmcif" {
        return Err("Only viewer mmCIF files can be read".to_string());
    }
    let metadata = std::fs::metadata(path)
        .map_err(|error| format!("Could not inspect viewer file: {error}"))?;
    if metadata.len() > 256 * 1024 * 1024 {
        return Err("Viewer file exceeds the 256 MB safety limit".to_string());
    }
    std::fs::read_to_string(path).map_err(|error| format!("Could not read viewer file: {error}"))
}

#[tauri::command]
fn get_structure_path(app: AppHandle, name: String) -> Result<String, String> {
    if !name
        .chars()
        .all(|character| character.is_ascii_alphanumeric() || character == '-' || character == '_')
    {
        return Err("Invalid built-in structure name".to_string());
    }
    let resource_dir = app
        .path()
        .resource_dir()
        .unwrap_or_else(|_| PathBuf::from("."));
    let pdb_path = resource_dir.join("pdb").join(format!("{name}.pdb"));
    if pdb_path.exists() {
        Ok(pdb_path.to_string_lossy().to_string())
    } else {
        Err(format!("Built-in structure not found: {name}"))
    }
}

#[tauri::command]
fn list_builtin_structures(app: AppHandle) -> Result<Vec<String>, String> {
    let resource_dir = app
        .path()
        .resource_dir()
        .unwrap_or_else(|_| PathBuf::from("."));
    let pdb_dir = resource_dir.join("pdb");
    let mut names = Vec::new();
    if let Ok(entries) = std::fs::read_dir(pdb_dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.extension().is_some_and(|extension| extension == "pdb") {
                if let Some(stem) = path.file_stem() {
                    names.push(stem.to_string_lossy().to_string());
                }
            }
        }
    }
    names.sort();
    Ok(names)
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_fs::init())
        .plugin(tauri_plugin_shell::init())
        .manage(SidecarState {
            child: Mutex::new(None),
            pending: Mutex::new(HashMap::new()),
            expected_termination: Mutex::new(false),
        })
        .setup(|_| {
            // Allows the build script to exercise full Tauri/plugin startup
            // without leaving a GUI process running in CI.
            if std::env::var_os("CAPSACIN_STARTUP_CHECK").is_some() {
                std::process::exit(0);
            }
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![
            inspect_structure,
            prepare_preview,
            run_slice,
            cancel_operation,
            open_pdb_dialog,
            save_pdb_dialog,
            read_text_file,
            get_structure_path,
            list_builtin_structures,
        ])
        .on_window_event(|window, event| {
            if let tauri::WindowEvent::Destroyed = event {
                let _ = kill_sidecar(window.app_handle());
            }
        })
        .run(tauri::generate_context!())
        .expect("error while running capSACIN Studio");
}
