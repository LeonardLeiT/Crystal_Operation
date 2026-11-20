# CrystalOps

CrystalOps 是一个用于粒子级结构建模、分析、可视化和文件交换的 Python 工具包。所有来源的结构——代码生成或外部文件读取——都使用同一种 `CrystalStructure` 对象。

## 设计

- `crystalops.core`：共享的 `CrystalStructure`、`UnitCell` 和 `Element` 数据模型；统一坐标、模拟盒、周期条件和逐原子属性。
- `crystalops.modeling`：结构建模，包含 `builders`、`operations` 和 `lattices`；提供普通晶体、多晶、TPMS 生成器以及平移、旋转、扩胞、筛选和合并等操作。
- `crystalops.analysis`：结构及模拟结果分析，包括各向异性、周期晶界识别、Gaussian 占据界面重构和局部曲率计算；可选 PTM 晶粒恢复。
- `crystalops.visualization`：结构、曲面及分析结果的可视化，承接原 `crysplot` 绘图接口和 recipes。
- `crystalops.io`：外部文件读写；目前支持 LAMMPS data 的 `atomic` atom style。

建模、读写和分析共用 `core` 中的数据模型。分析结果由分析模块返回，可交给可视化模块显示；导入基础结构模型不需要绘图依赖。
保留的 recipes 是独立脚本，运行时需准备各自的输入数据和额外依赖。

旧的 `crystalops.builders`、`crystalops.operations`、`crystalops.lattices` 和 `crysplot` 导入路径保留为兼容入口，现有 `from crystalops import ...` 写法也继续可用。新代码建议使用按职责划分的模块：

```python
from crystalops.core import CrystalStructure, UnitCell
from crystalops.modeling import CrystalBuilder, FCC
from crystalops.visualization import StructurePlotter
```

最后一行需要安装可选绘图依赖。示例按相同职责组织，参见 [example 说明](example/README.md)。

## 安装

```bash
pip install -e .
```

需要绘图功能时：

```bash
pip install -e ".[plot]"
```

晶界曲率分析及可选的晶粒恢复依赖：

```bash
pip install -e ".[analysis,grain-segmentation,visualization]"
```

晶界分析函数、输入标签要求和结果格式见 [分析模块说明](crystalops/analysis/README.md)。

## 创建并写出晶体

```python
from crystalops.modeling import CrystalBuilder, FCC

structure = CrystalBuilder(
    unit_cell=FCC("Cu"),
    repeats=(10, 10, 10),
    lattice_constant=3.61,
).build()
structure.write("fcc.data")
```

## 读取、操作并重新写出

```python
from crystalops.core import CrystalStructure
from crystalops.modeling.operations import rotation_matrix

structure = CrystalStructure.read("fcc.data")
structure = structure.repeated((2, 1, 1))
structure = structure.rotated(rotation_matrix(z_angle=30.0))
structure = structure.translated([2.0, 0.0, 0.0])
structure.write("fcc_modified.data")
```

## 生成 Kelvin 多晶

```python
from crystalops.modeling import FCC, KelvinCrystal

structure = KelvinCrystal(
    unit_cell=FCC("Cu"),
    size=5.5,
    lattice_constant=3.61,
    preset="DSC1",
    random_seed=42,
).build()
structure.write("kelvin.data")
```

## 生成 Weaire-Phelan 多晶

按照 Yang 等人在 *Acta Materialia* 296 (2025) 121289 中的 A15 位点和
Table 1 晶体取向，生成包含两个 Dod 晶粒和三对 Tet 晶粒的周期前驱体：

```python
from crystalops.modeling import FCC, WeairePhelanCrystal

sigma3 = WeairePhelanCrystal(
    unit_cell=FCC("Cu"),
    size=14.0,          # 论文记号中的晶粒尺寸 d，盒长 a = 2d
    lattice_constant=3.61,
    twin_type="sigma3",
).build()
sigma3.write("weaire_phelan_sigma3.data")

sigma5 = WeairePhelanCrystal(
    unit_cell=FCC("Cu"),
    size=14.0,
    lattice_constant=3.61,
    twin_type="sigma5",
).build()
sigma5.write("weaire_phelan_sigma5.data")
```

`sigma3` 是论文中会向 Euler-catenoid 网络转变的构型，`sigma5` 是不含
Σ3 CTB 网络的对照构型。两者默认都使用 Table 1 的立方 Dod 取向；只有显式
设置 `dod_orientation_seed` 时才会生成实验性的随机 Dod。

构型使用加权周期 Voronoi（Laguerre）分割，使八个初始胞元等体积；
`metadata["cell_labels"]` 保存 Dod/Tet 位点名称，`grain_ids` 保存晶粒编号。
与 `KelvinCrystal` 相同，每个完整胞元先在自己的中心构建，再平移到 A15
位点；跨越周期边界的 Dod/Tet 不会被裁剪或 wrap 回主盒。
每对 Tet 使用共同的 CTB 面锚点作为晶格原点，使论文给出的 Σ3/Σ5
取向同时具有正确的平移相位；随后按照周期最小镜像距离删除小于默认
`0.5 * nearest_neighbor_distance` 的重复或过近原子。该阈值可通过
`overlap_distance_factor` 调整。

## 生成 Schwarz D-SCI 晶体

使用 `SchwarzCrystal` 构建六区域 D-SCI，并返回统一的
`CrystalStructure`：

```python
from crystalops.modeling import SchwarzCrystal

structure = SchwarzCrystal(
    size=8,
    lattice_constant=3.61,
).build()
print(structure.metadata["domain_counts"])
structure.write("schwarz.data")
```

用同一个六区域 DSC 构型可以生成不同轴比的周期模型。

```python
dsc = SchwarzCrystal(size=30, schwarz_type="DSC", z_ratio=1.4).build()
dsc.write("dsc_1-1-1.4.data")
```

`z_ratio` 设置周期盒的 `Lz/Lx`。两中心域共享一个 FCC 取向，四过渡域
保留相对于中心的四种 Σ3 孪晶取向。非立方盒按各自盒长重定位孪晶面，
使其落在实际 FCC {111} 晶面上，并按镜像关系匹配两侧原子层的堆垛。
过渡域之间的曲面用原 DSC 的 Diamond 方程判定。为使五种取向与周期盒相容，
`size` 和 `size * z_ratio` 都必须是 3 的正整数倍。原子折回主盒后，
按周期最小镜像去除间距小于半个 FCC 最近邻距离的过近原子。输出为弛豫前构型。

同一个类也可以按 Schwarz-P 方程生成 PSC：

```python
psc = SchwarzCrystal(
    size=8,
    lattice_constant=3.61,
    schwarz_type="PSC",
    random_seed=42,
).build()
psc.write("schwarz_psc.data")
```

PSC 同时包含两个晶粒：`P <= 0` 对应晶粒 1，`P > 0` 对应晶粒 2。
用户最多提供晶粒 1 的一个 `orientation`；未提供时先随机生成可公度取向。
晶粒 2 根据晶粒 1 和周期盒自动选取兼容取向。必要时盒长会在
`max_box_expansion` 范围内从目标 `size` 略微向上调整。
晶粒编号保存在 `properties["grain_ids"]`，`atom_types` 始终只表示元素类型。

这里的晶粒编号属于内存结构的逐原子属性。当前 LAMMPS `atomic` data 读写不会保存或恢复 `grain_ids`；需要分析晶粒归属时，应保留生成器返回的结构，或通过支持该属性的读入接口显式载入标签，并按原子 ID 对齐。

生成器只负责创建结构，不会自动写文件。因此同一个结果可以继续操作，也可以输出为不同文件格式。
