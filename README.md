# CrystalOps

CrystalOps 是一个用于构建、操作和交换粒子级晶体结构的 Python 工具包。所有来源的结构——代码生成或外部文件读取——都会转换为同一种 `CrystalStructure` 对象。

## 设计

- `crystalops.core`：`CrystalStructure` 和 `UnitCell` 数据模型。
- `crystalops.operations`：平移、旋转、扩胞、筛选、合并和原子类型操作。
- `crystalops.builders`：普通晶体、Voronoi/Kelvin 多晶、Math D-SCI 和 TPMS 结构生成器。
- `crystalops.io`：外部文件读写；目前支持 LAMMPS data 的 `atomic` atom style。

所有通用绘图接口都放在独立的 `crysplot` 包中。

## 安装

```bash
pip install -e .
```

需要绘图功能时：

```bash
pip install -e ".[plot]"
```

## 创建并写出晶体

```python
from crystalops import CrystalBuilder, FCC

structure = CrystalBuilder(
    unit_cell=FCC("Cu"),
    repeats=(10, 10, 10),
    lattice_constant=3.61,
).build()
structure.write("fcc.data")
```

## 读取、操作并重新写出

```python
from crystalops import CrystalStructure
from crystalops.operations import rotation_matrix

structure = CrystalStructure.read("fcc.data")
structure = structure.repeated((2, 1, 1))
structure = structure.rotated(rotation_matrix(z_angle=30.0))
structure = structure.translated([2.0, 0.0, 0.0])
structure.write("fcc_modified.data")
```

## 生成 Kelvin 多晶

```python
from crystalops import FCC, KelvinCrystal

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
from crystalops import FCC, WeairePhelanCrystal

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
from crystalops import SchwarzCrystal

structure = SchwarzCrystal(
    size=8,
    lattice_constant=3.61,
).build()
print(structure.metadata["domain_counts"])
structure.write("schwarz.data")
```

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

生成器只负责创建结构，不会自动写文件。因此同一个结果可以继续操作，也可以输出为不同文件格式。
