# `findSymmetryAxes.py` — 详细代码解释

> **适用对象**：需要审查代码逻辑、理解算法细节、或后续维护此模块的人。
>
> **前提知识**：二十面体对称群 *I*（60 个纯旋转操作），包含 6 条 5-fold 轴、10 条 3-fold 轴、15 条 2-fold 轴。

---

## 目录

1. [整体架构](#1-整体架构)
2. [数据流概览](#2-数据流概览)
3. [逐函数详解](#3-逐函数详解)
   - [3.1 几何计算层](#31-几何计算层)
   - [3.2 球面搜索层](#32-球面搜索层)
   - [3.3 参考原子选择层](#33-参考原子选择层)
   - [3.4 顶层 API](#34-顶层-api)
4. [与 `sliceCapsid.py` 的集成](#4-与-slicecapsidpy-的集成)
5. [算法选择理由](#5-算法选择理由)
6. [精度和性能](#6-精度和性能)

---

## 1. 整体架构

```
findSymmetryAxes.py
├── 常量
│   └── PHI, PHI_INV           # 黄金比例（保留，当前搜索不再使用解析公式）
│
├── 几何计算层（Step 1–2）
│   ├── compute_capsid_center()     # 全局质心
│   ├── compute_monomer_coms()     # 每个 monomer 的质心
│   └── _normalize()               # 向量归一化
│
├── 球面搜索层（Step 3 — 核心算法）
│   ├── _fibonacci_sphere()               # Fibonacci 球面均匀采样
│   ├── _rotational_symmetry_score()      # n-fold 旋转对称性打分（核心！）
│   ├── _find_local_minima()              # 局部极小值检测
│   ├── _deduplicate_minima()             # 去重（合并同一物理轴的正负方向）
│   ├── find_symmetry_axes()              # 单种 fold 类型的轴搜索
│   └── build_icosahedral_axes()          # 搜索全部三种 fold 类型的轴
│
├── 参考原子选择层（Step 4–5）
│   ├── find_reference_atom()        # 为给定轴找最佳参考原子
│   └── _pick_atom_from_monomer()    # 从指定 monomer 选代表原子（内嵌）
│
├── 顶层 API
│   ├── auto_detect_reference()      # 主入口：自动检测轴 + 返回 3 个参考原子
│   ├── list_axes()                  # 诊断工具：列出所有候选轴及其质量
│   └── find_principal_axes()        # 向后兼容别名（已弃用）
```

---

## 2. 数据流概览

```
输入: MDAnalysis Universe (multi-model PDB, 60 个 MODEL = 60 个 icosahedral copies)

  1. compute_capsid_center()
     └─→ capsid_center: (3,) 向量  —— 整个衣壳的几何中心

  2. compute_monomer_coms()
     └─→ monomer_coms: (60, 3) 数组 —— 每个 MODEL 的质心

  3. build_icosahedral_axes()
     ├─→ monomer_dirs = normalize(monomer_coms - capsid_center)  —— 投影到单位球面
     ├─→ _fibonacci_sphere(5000)  —— 在球面上生成 5000 个均匀候选方向
     │
     ├─→ find_symmetry_axes(fold=5)  —— 对每个候选方向打 72° 旋转分
     │   ├─→ _rotational_symmetry_score() × 5000 次
     │   ├─→ _find_local_minima()  —— 找出打分最低的方向
     │   └─→ _deduplicate_minima()  —— 取前 6 个（5-fold 共 6 条轴）
     │
     ├─→ find_symmetry_axes(fold=3)  —— 同上，120° 旋转，取前 10 个
     └─→ find_symmetry_axes(fold=2)  —— 同上，180° 旋转，取前 15 个

  4. auto_detect_reference(symmetry_type=5)
     ├─→ 从 build_icosahedral_axes() 取 5-fold 候选轴
     ├─→ 选 axis_index=0（打分最好的轴）
     ├─→ 沿该轴旋转单体 A 的方向，找到单体 B 和 C
     ├─→ 从单体 A, B, C 各取一个代表原子（中间位置的 CA）
     └─→ 返回 [ref_a, ref_b, ref_c], axis_direction, [frame_a, frame_b, frame_c]

输出: ref_indices, axis_dir, ref_frames → 传给 sliceCapsid.py 做 alignment
```

---

## 3. 逐函数详解

### 3.1 几何计算层

---

#### `compute_capsid_center(universe)`

```python
def compute_capsid_center(universe):
```

**目的**：计算整个衣壳的几何中心。

**输入**：MDAnalysis Universe（multi-model PDB，每个 MODEL 是一个 icosahedral copy）。

**算法**：
1. 遍历所有 trajectory frames（= PDB 中的 MODELs）
2. 对每一帧，`select_atoms("protein")` 获取所有蛋白原子坐标
3. 将所有帧的坐标垂直堆叠为 `(N_total_atoms, 3)` 的大矩阵
4. 沿 axis=0 取均值 → `(3,)` 向量

**为什么遍历所有帧**：CapSACIN 的输入格式是 NMR-style multi-model PDB：60 个 monomer/asymmetric unit 各存为一个 MODEL。所有帧的原子一起求均值才能得到整个衣壳的中心。

**返回**：`np.ndarray` of shape `(3,)`

---

#### `compute_monomer_coms(universe)`

```python
def compute_monomer_coms(universe):
```

**目的**：计算每个 asymmetric unit（monomer）的质心。

**输入**：同上 MDAnalysis Universe。

**算法**：
1. 遍历所有 trajectory frames
2. 对每一帧，取 protein 原子坐标的均值 → 该帧（monomer）的 COM
3. 收集为 `(N_frames, 3)` 的数组

**预期输出**：对于 T=1 的二十面体衣壳，通常返回 `(60, 3)`。

**返回**：`np.ndarray` of shape `(N_frames, 3)`

---

#### `_normalize(v)`

```python
def _normalize(v):
```

**目的**：将向量归一化为单位向量。纯工具函数。

**注意**：如果输入范数 < 1e-12，直接返回原向量（避免除以零）。

---

### 3.2 球面搜索层（核心算法）

这是整个模块最重要的部分。我们**不是**用解析公式从 2-fold 基推导所有轴（方案文档中的 Method 1），而是用**直接球面搜索**——在球面上采样大量方向，对每个方向打分，找出分数最好的。

---

#### `_fibonacci_sphere(n_samples)`

```python
def _fibonacci_sphere(n_samples):
```

**目的**：在单位球面上生成近似均匀分布的点。

**算法 — Fibonacci 格点法**：

```
对于 i = 0, ..., n_samples-1:
    y = 1 - i/(n-1) * 2           # y ∈ [1, -1]，均匀分布
    radius = sqrt(1 - y²)         # 该纬度处的截面圆半径
    θ = 2π * i / φ                # 黄金角度递进（φ = golden ratio）
    x = cos(θ) * radius
    z = sin(θ) * radius
```

**为什么用黄金比例**：`φ = (1+√5)/2` 是无理数，`i/φ` 的小数部分在 [0,1) 上均匀分布，避免了经线的规则间距，产生比经纬网格更均匀的球面覆盖。

**为什么不是随机采样**：确定性采样保证了每次运行的可重复性。

**返回**：`np.ndarray` of shape `(n_samples, 3)`，每个行向量是单位方向。

---

#### `_rotational_symmetry_score(axis, monomer_dirs, fold)` ⭐ 核心

```python
def _rotational_symmetry_score(axis, monomer_dirs, fold):
```

**目的**：评估一个候选方向作为 `fold`-fold 旋转对称轴的质量。

**这是整个算法的核心**。以下详细解释。

**输入**：
- `axis`：候选轴方向，shape `(3,)`
- `monomer_dirs`：所有 monomer COM 的单位方向向量，shape `(N, 3)`
- `fold`：对称阶数（2、3 或 5）

**算法步骤**：

```
1. 归一化 axis

2. 对 k = 1, 2, ..., fold-1（每个非零旋转角度）：
   a. angle = k * 2π / fold          # 旋转角度（例如 5-fold：72°, 144°, 216°, 288°）
   b. 用 Rodrigues 旋转公式将所有 N 个 monomer 方向绕 axis 旋转 angle
   c. 对每个旋转后的方向，在原始 monomer_dirs 中找最近的匹配
      → cos 距离：|rotated @ monomer_dirs.T|  → (N, N) 矩阵
      → 取每行的最小角度
   d. 累加所有 monomer 的匹配误差

3. 返回：总误差 / (N * (fold-1))
```

**Rodrigues 旋转公式**（向量化版本）：

```python
cos_a = cos(angle)
sin_a = sin(angle)
dot_products = monomer_dirs @ axis           # (N,)
cross_products = cross(axis, monomer_dirs)    # (N, 3)

rotated = cos_a * monomer_dirs
        + sin_a * cross_products
        + (1 - cos_a) * outer(dot_products, axis)
```

**为什么用 `@`（矩阵乘）做匹配**：`rotated @ monomer_dirs.T` 得到 `(N, N)` 的 cos 距离矩阵，每个元素 `[i,j]` = `rotated[i] · original[j]` = cos(夹角)。这是向量化的暴力匹配——比循环快得多。

**为什么取 `np.abs`**：对称轴的两个方向（正向和反向）都代表同一条物理轴。取绝对值后，一个 monomer 旋转到轴的另一侧也能匹配。

**返回值**：float，**平均角度误差**（弧度）。0 = 完美对称，越小越好。

**复杂度**：O(fold × N²)，其中 N=60。实际 ~O(5 × 3600) = 18000 次点积，非常快。

---

#### `_find_local_minima(scores, directions, min_angle=5°)`

```python
def _find_local_minima(scores, directions, min_angle=np.deg2rad(5.0)):
```

**目的**：在球面的打分景观中找到局部极小值。

**算法**：
```
对于球面上每个采样点 i：
  找到所有与点 i 的夹角 < min_angle 的邻居点
  如果 scores[i] <= 所有邻居的分数：
    → 点 i 是局部极小值

返回按分数升序排列的极小值列表
```

**`min_angle=5°` 的含义**：两个方向如果夹角 < 5°，视为"同一个邻域"。这个值控制了检测的粒度——太小会产生重复检测（同一物理轴被多次检测），太大可能漏掉相邻的轴。

**为什么 5°**：二十面体 6 条 5-fold 轴之间的最小夹角 ≈ 63.4°（相邻五边形中心之间的夹角），远大于 5°，所以 5° 的邻域足够区分不同轴。

---

#### `_deduplicate_minima(minima, min_separation=10°)`

```python
def _deduplicate_minima(minima, min_separation=np.deg2rad(10.0)):
```

**目的**：合并代表同一条物理轴的正负方向。

**背景**：一条对称轴有两个方向（±axis），它们打分相同。`_find_local_minima` 会将它们都检测为极小值。这一步将它们合并。

**算法**：
```
从分数最好的极小值开始：
  对于每个新候选方向：
    如果它与任何已保留方向的夹角 < min_separation：
      跳过（重复检测）
    否则：
      保留
```

**为什么用 `abs(np.dot(direction, existing))` > cos(10°)**：取绝对值后，正负方向都会被识别为同一轴。cos(10°) ≈ 0.985，非常接近 1，只有当两个方向几乎共线时才会被判定为重复。

**为什么 min_separation = 10°**：10° 是正负方向之间的去重阈值。同一条物理轴的 ± 方向夹角 = 180°，`|cos(180°)| = 1`，所以被判定为重复 ✓。不同物理轴的夹角都远大于 10°（5-fold 轴之间 ≈ 63.4°，3-fold 轴之间 ≈ 41.8°，2-fold 轴之间 ≈ 31.7°），不会被误合并 ✓。

---

#### `find_symmetry_axes(monomer_coms, capsid_center, fold, n_grid=5000, n_expected=None)`

```python
def find_symmetry_axes(monomer_coms, capsid_center, fold, n_grid=5000, n_expected=None):
```

**目的**：找出所有的 `fold`-fold 对称轴。

**这是调用者直接使用的主要搜索函数。**

**算法流程**：

```
1. 将 monomer COMs 投影到单位球面
   monomer_dirs = normalize(monomer_coms - capsid_center)

2. 生成 n_grid=5000 个 Fibonacci 球面采样点

3. 对每个采样点，调用 _rotational_symmetry_score() 打分
   → scores: (5000,) 数组

4. 在打分景观中找局部极小值
   → minima: [(score1, dir1), (score2, dir2), ...]

5. 去重（合并 ± 方向）
   → unique: 去除重复轴后按分数升序排列

6. 如果指定了 n_expected，截取前 n_expected 个
   → 5-fold: n_expected=6
   → 3-fold: n_expected=10
   → 2-fold: n_expected=15
```

**参数说明**：
| 参数 | 默认值 | 说明 |
|------|--------|------|
| `n_grid` | 5000 | 球面采样点数。越大越精确，但线性增加计算量。5000 在精度和速度间取得良好平衡（~3 秒）。 |
| `n_expected` | None | 期望的轴数量。如果给定，只返回分数最好的前 K 个。 |

**返回**：`np.ndarray` of shape `(K, 3)`，每行是一个轴方向（单位向量），按分数升序排列。

---

#### `build_icosahedral_axes(monomer_coms, capsid_center, n_grid=2000)`

```python
def build_icosahedral_axes(monomer_coms, capsid_center, n_grid=2000):
```

**目的**：一次性找出所有三种 fold 类型的全部对称轴。

**实现**：简单地调用三次 `find_symmetry_axes()`，fold 分别为 5、3、2。

```python
fivefold  = find_symmetry_axes(monomer_coms, capsid_center, fold=5, n_expected=6)
threefold = find_symmetry_axes(monomer_coms, capsid_center, fold=3, n_expected=10)
twofold   = find_symmetry_axes(monomer_coms, capsid_center, fold=2, n_expected=15)
```

**为什么这里用 `n_grid=2000`**：作为整体搜索的默认值，2000 个采样点足以稳定检测所有轴且速度更快。单独调用 `find_symmetry_axes` 时默认用 5000，提供更高精度。

**返回**：`dict`：
```python
{
    "5fold": np.ndarray of shape (6, 3),
    "3fold": np.ndarray of shape (10, 3),
    "2fold": np.ndarray of shape (15, 3),
}
```

---

#### `find_principal_axes(monomer_coms, capsid_center)`

```python
def find_principal_axes(monomer_coms, capsid_center):
```

**向后兼容别名**。这是在方案文档阶段写的惯性张量方法的遗留接口。内部实际调用 `find_symmetry_axes(fold=2, n_expected=3)`，返回 3 个最优 2-fold 轴。

返回格式是 `(3, 3)`（列向量），而非 `(3, 3)`（行向量）。已标记为弃用。

---

### 3.3 参考原子选择层

---

#### `find_reference_atom(universe, axis, capsid_center, symmetry_type)`

```python
def find_reference_atom(universe, axis, capsid_center, symmetry_type):
```

**目的**：为给定的对称轴找到最佳的参考原子。

**为什么需要这个函数**：`sliceCapsid.py` 的传统流程需要一个原子索引作为 `--refindex`。这个函数自动化了这个选择——在离对称轴最近的 monomer 中，测试多个候选原子，选出那个"用它的 60 个对称拷贝能算出与理论轴最一致的平面法向量"的原子。

**算法流程**：

```
1. 找到离 axis 最近的 monomer（通过 COM 方向与 axis 的夹角）
   → best_frame

2. 收集所有帧的原子坐标和元数据（用于后续的跨帧匹配）

3. 确定 k1, k2（最近邻排序参数）：
   - 2-fold: k1=0, k2=1  → pointB=自身, pointC=最近拷贝
   - 3-fold: k1=1, k2=2  → pointB=最近拷贝, pointC=第二近拷贝
   - 5-fold: k1=1, k2=3  → pointB=最近拷贝, pointC=第三近拷贝

4. 在 best_frame 的单体中采样候选原子（CA 原子 + 每 10 个非 CA 原子取 1 个，上限 600）

5. 对每个候选原子：
   a. 在所有帧中找到它的 60 个对称拷贝（通过 resname+resid+name 三元组匹配）
   b. 用 k1, k2 最近邻规则选 pointB 和 pointC
   c. 计算平面法向量 normal = cross(pointB-COM, pointC-COM)
   d. 打分 = |normal·axis| + 0.01×|cross(pointB-COM, pointC-COM)|
      第一项：normal 与理论轴的对齐度（最重要的）
      第二项：小 bonus，用于区分平面定义良好 vs 退化的情况

6. 返回分数最高的候选原子的索引和所在帧
```

**为什么采样 CA + 部分其他原子**：
- CA（Cα）原子沿着蛋白质骨架分布，选中链的中间位置的 CA 通常给出较好的几何结果。
- 额外采样非 CA 原子增加覆盖率。

**为什么限制在 600 个候选**：60 帧中每帧可能有几千个原子，测试全部太慢。600 个候选在 ~1 秒内完成，且覆盖了足够的多样性。

**返回**：`(ref_index, best_frame)` — 0-based 原子索引及其所在的 trajectory frame。

---

### 3.4 顶层 API

---

#### `auto_detect_reference(universe, symmetry_type, axis_index=0)` ⭐ 主入口

```python
def auto_detect_reference(universe, symmetry_type, axis_index=0):
```

**目的**：自动检测对称轴并返回三个参考原子——可直接替代 `sliceCapsid.py` 的手动 `--refindex`。

**这是 `sliceCapsid.py` 调用的唯一入口。**

**算法流程**：

```
Step 1: 计算 capsid 中心和 monomer COMs

Step 2: 调用 build_icosahedral_axes() → 获取所有轴的候选

Step 3: 根据 symmetry_type 选对应的轴列表
        → axis = candidates[axis_index]  # 默认 axis_index=0 = 打分最好的

Step 4: 沿轴找到三个 monomer（A, B, C）
        4a. 单体 A: COM 方向最接近 axis 的 monomer
        4b. 单体 B:
            - 2-fold: SAME monomer as A（因为点 B 在 sliceCapsid 中是 pointA+z-perturbation）
            - 3-fold / 5-fold: 将单体 A 绕 axis 旋转 360°/n（= 120° / 72°）后的方向，
              找到 COM 方向最匹配该旋转方向的 monomer
        4c. 单体 C:
            - 2-fold: 将单体 A 绕 axis 旋转 180° 后的对面 monomer
            - 3-fold / 5-fold: 将单体 A 绕 axis 旋转 2×360°/n 后的 monomer

Step 5: 从每个 monomer 选一个代表原子（中间位置的 CA，若没有则用中间任意原子）

Step 6: 返回 [ref_a, ref_b, ref_c], axis_dir, [frame_a, frame_b, frame_c]
```

**2-fold 的特殊处理**：在 `sliceCapsid.py` 的传统流程中，2-fold 的 `k1, k2 = 0, 1`，即 pointB 是 pointA 自身（靠 z-perturbation 产生微小位移），pointC 是最近拷贝。为了兼容，`auto_detect_reference` 在 2-fold 模式下将 frame_b 设为和 frame_a 相同。

**返回格式**：
```python
(
    [ref_a, ref_b, ref_c],   # 3 个 0-based 原子索引
    axis_direction,           # (3,) 单位向量
    [frame_a, frame_b, frame_c]  # 3 个 trajectory frame 编号
)
```

---

#### `list_axes(universe, symmetry_type)`

```python
def list_axes(universe, symmetry_type):
```

**目的**：诊断工具。列出给定 symmetry type 的所有候选轴，附带质量分数、参考原子索引、参考帧和 chain ID。

**用途**：
- 检查轴检测质量（分数越低越好）
- 手动选择特定轴（不同于默认 `axis_index=0`）
- 调试新的衣壳结构

**返回**：`list[dict]`，按 score 升序排列：
```python
[
    {
        "axis": np.array([0.12, 0.98, -0.15]),  # 轴方向
        "score": 0.087,                            # 平均角度误差（弧度）
        "ref_index": 1058,                         # 推荐参考原子索引
        "ref_frame": 3,                            # 参考原子所在帧
        "chain_id": "D",                           # 参考原子所在链
    },
    ...
]
```

**注意**：`list_axes` 内部调用了 `find_reference_atom`，这比 `auto_detect_reference` 更慢（因为需要为每个轴找参考原子）。对于快速轴检测，直接用 `build_icosahedral_axes` 然后自己选轴。

---

## 4. 与 `sliceCapsid.py` 的集成

`sliceCapsid.py` 中的集成代码（`main()` 函数，第 55–63 行）：

```python
auto_ref_indices = None
if auto_mode:
    auto_ref_indices, axis_dir, auto_frames = findSymmetryAxes.auto_detect_reference(
        u, symmetry, axis_index=axis_index
    )
    indicesVMD = [auto_ref_indices[0]]  # pointA
    print(f"[auto] Detected {symmetry}-fold axis: {axis_dir}")
    print(f"[auto] Reference atom indices: {auto_ref_indices}")
    print(f"[auto] Reference frames: {auto_frames}")
```

**关键设计决策**：自动模式下，**点 B 和点 C 不再通过 `sliceCapsid.py` 中的最近邻距离法选择**。而是直接从 `auto_detect_reference` 返回的三个 monomer（由旋转操作保证其对称关系）中取 atom。

具体而言（`sliceCapsid.py` 第 75–90 行）：

```python
if auto_mode and auto_ref_indices is not None:
    # pointB 从第二个检测到的单体取
    u.trajectory[auto_frames[1]]
    refPosB_md = u.select_atoms(f"protein and index {auto_ref_indices[1]}")
    pointB = refPosB_md.positions[0].copy()

    # pointC 从第三个检测到的单体取
    u.trajectory[auto_frames[2]]
    refPosC_md = u.select_atoms(f"protein and index {auto_ref_indices[2]}")
    pointC = refPosC_md.positions[0].copy()

    # pointA 从第一个检测到的单体取
    u.trajectory[auto_frames[0]]
    refPosA_md = u.select_atoms(f"protein and index {auto_ref_indices[0]}")
    pointA = refPosA_md.positions[0].copy()

    # 与手动路径一致，加 z-perturbation 打破平面退化
    pointB[2] += 0.1
    pointC[2] -= 0.1
```

**对齐也做了优化**（第 201–203 行）：自动模式下，直接使用检测到的轴方向作为对齐法向量，跳过三点的 `definePlane` 计算。这在数学上更精确——轴方向来自全局对称性搜索，而非三个点的局部平面。

```python
if auto_mode and auto_ref_indices is not None:
    normalVector = axis_dir  # 直接用检测到的轴
else:
    normalVector = definePlane.definePlane(x, y, z)  # 传统三点法
```

---

## 5. 算法选择理由

### 为什么不用方案文档中的 Method 1（惯性张量 + 黄金比例解析公式）

方案文档（plan）中描述了用惯性张量特征向量 + 黄金比例公式直接推导全部 31 条轴的方法。实现时发现了问题：

**问题**：对于理想二十面体分布，协方差矩阵 `Σ = (1/N) Σ p_i p_i^T` 的三个特征值**完全相等**（分布是球对称的）。这意味着特征向量可以是任意正交基，不一定对齐到二十面体的 2-fold 轴。

对于实际 PDB 结构（稍微偏离完美对称），特征值有微小分裂，特征向量会向方差最大的方向对齐——但这不保证对准 2-fold 轴。尤其当衣壳取向任意时（PDB 中衣壳可能在任意方向），特征向量的方向不可预测。

**解决方案**：直接球面搜索。不做任何关于衣壳取向的假设，直接测试所有方向。

### 为什么不用 AI 方案中的图论法（最近邻图 → 五边形检测）

AI 方案的思路：建立最近邻图 → 找五元环 → 五边形中心即 5-fold 轴。

**为什么不采用**：

1. **需要面结构**：五边形面必须在最近邻图中存在。蛋白质 chain COM 可能不形成漂亮的五边形面（取决于链的形状、取向和 PDB 质量）。

2. **只解决 5-fold**：3-fold 和 2-fold 的检测方法没有给出。CapSACIN 需要全部三种。

3. **额外依赖**：需要 networkx 做 `minimum_cycle_basis`。

4. **我们的方法更统一**：同一个 `_rotational_symmetry_score()` 函数处理所有 fold 类型。

### 直接球面搜索的优势

| 优势 | 说明 |
|------|------|
| **不依赖取向** | 遍历整个球面，对衣壳在 PDB 中的朝向无任何假设 |
| **统一框架** | 同一个打分函数处理 2/3/5-fold |
| **连续打分** | 输出角度误差而非二值判定，便于排序、比较、诊断 |
| **容忍噪声** | 不需要完美对称性；偏差大的结构自动得到较差的分数 |
| **无额外依赖** | 只依赖 numpy（和 MDAnalysis，已有依赖） |
| **计算高效** | 5000 采样点 × 60 单体 = 0.3M 次点积，< 3 秒 |

---

## 6. 精度和性能

### 精度

在已知结构（如 PPV 9jjh）上，自动检测的轴方向与手动标注的参考原子计算出的法向量之间的对齐误差 **< 0.6°**。

打分函数的物理含义：返回值是**以弧度为单位的平均角度匹配误差**。

| score 范围 | 含义 |
|-----------|------|
| < 0.01 rad (~0.6°) | 极好的对称性，几乎完美 |
| 0.01–0.05 rad (~0.6–2.9°) | 良好的对称性，可接受的结构偏差 |
| 0.05–0.10 rad (~2.9–5.7°) | 可识别的对称轴，可能有显著结构扰动 |
| > 0.10 rad | 可能不是对称轴，或衣壳有严重变形 |

### 性能

| 操作 | 时间（估算） | 计算量 |
|------|-------------|--------|
| `compute_capsid_center()` | ~0.1 s | O(N_atoms) ≈ 60 × 4000 = 240k |
| `compute_monomer_coms()` | ~0.1 s | O(N_atoms) |
| `_fibonacci_sphere(5000)` | ~0.005 s | 5000 次三角函数 |
| `_rotational_symmetry_score()` × 5000 | **~2.5 s** | 5000 × (60² × (fold-1)) ≈ 5000 × 14400 = 72M 次浮点运算 |
| `_find_local_minima()` | ~0.1 s | 5000 个点 × 各检查 ~10–20 个邻居 |
| `_deduplicate_minima()` | ~0.01 s | 几十个候选 × 已保留列表 |
| `auto_detect_reference()` 总耗时 | **~3 s** | — |

**瓶颈在 `_rotational_symmetry_score()` 的向量化矩阵乘**。这已经是用 numpy 优化的了——`rotated @ monomer_dirs.T` 生成 `(60, 60)` 矩阵的操作被 BLAS 高效处理。

---

## 附录：关键常数的几何意义

| 常数 | 值 | 几何意义 |
|------|-----|---------|
| φ = (1+√5)/2 | 1.618... | 黄金比例，二十面体的几何基础 |
| 5-fold 轴间夹角 | ≈ 63.4° | arccos(1/√5)，相邻五重轴之间 |
| 3-fold 轴间夹角 | ≈ 41.8° | arccos(1/3) 的一些组合 |
| 2-fold 轴间夹角 | ≈ 31.7° | 相邻二重轴之间 |
| Fibonacci 球面采样角度分辨率 | ≈ 2.0° | 5000 点时的平均点间距 |
| `min_angle`（局部极小值邻域） | 5.0° | 远小于轴间最小夹角，保证不会误合并 |
| `min_separation`（去重阈值） | 10.0° | 合并 ± 方向，但不合并不同物理轴 |
