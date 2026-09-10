"""Entry point for `python -m sidecar`."""

if __name__ == "__main__":
    # tqdm creates multiprocessing locks. In a frozen executable the resource
    # tracker must be dispatched before it can import and start another server.
    from multiprocessing import freeze_support

    freeze_support()

    from sidecar.main import main

    main()
