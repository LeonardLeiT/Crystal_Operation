# C15-Laves 晶体构造

本目录只负责构造 C15-Laves 孪晶前驱体，不包含 Kelvin、A15/W-P 或其他
TPMS 候选代码。设计约束见 [DESIGN.md](DESIGN.md)。

## 生成

在仓库根目录运行：

```powershell
python example/mywork/TPMSMiner/C15/build_c15.py `
  --lattice-constant 3.61
```

默认写入 `C15/output/`：

```text
output/
├── 01_c15_all.data             # 完整：8个Z16 + 16个Z12
├── 02_z16_cell_origin.data     # 一个代表性Z16胞元，中心在(0,0,0)
├── 03_z12_cell_origin.data     # 一个代表性Z12胞元，中心在(0,0,0)
├── 04_c15_z16_only.data        # 仅8个Z16，保留全局相对位置
├── 05_c15_z12_only.data        # 仅16个Z12，保留全局相对位置
├── cell_manifest.json          # 24个胞元的完整位置/取向信息
├── cell_manifest.csv           # 便于表格查看的位置/取向信息
└── validation.json             # 单层法向、距离、CNA和LAMMPS验证
```

两个代表性单胞元被平移到笛卡尔原点，并使用以原点为中心的模拟盒。两个
`only`构型则保留完整C15盒子中的全局相对位置。

## 坐标与取向

manifest为每个胞元记录：

- `center_reduced`：主周期盒中的分数位置；
- `construction_center_reduced_unwrapped`：构造完整跨边界胞元时使用的周期像；
- `center_cartesian_unwrapped`：对应的笛卡尔位置，单位Å；
- `orientation_matrix_crystal_to_lab`：晶体坐标到实验室坐标的3×3取向矩阵，
  使用列向量约定 `v_lab = O @ v_crystal`；
- 孪晶伙伴、公共面外法向、共享锚点、Σ值和错配角。

默认使用单层孪晶规则：同一胞元类型和同一无向公共面法向只保留一个空间层。
原完美匹配中两组相隔半周期的重复Z12层已移除，因此当前为4对Z16、6对Z12，
合计10对；未配对的4个Z12晶粒使用固定随机种子的独立普通取向，以各自
胞心为晶格锚点，并排除距离Σ3/Σ5标称角小于1°的偶然近CSL关系。所有孪晶
组合使用完全相同的4个普通晶粒取向，避免给组合对比引入额外变量。

默认使用等体积C15-Laguerre胞元。添加 `--raw-voronoi` 可生成未经等体积
加权的原始C15 Voronoi对照。添加 `--label-grains` 可让完整构型以1–24原子
类型标记完整晶体；其余四个构型始终使用单一Cu原子类型。

尺寸固定为 `size=9.25`，Cu晶格常数3.61 Å时完整构型约7.5万原子、盒长
约96.32 Å。命令行不提供尺寸扫描参数，所有孪晶组合使用同一个固定尺寸。

## 孪晶组合筛选

按照“Σ3为正候选、Σ5为负对照、共享面中心晶格相位、固定C15几何”的规则，
批量生成不同组合：

```powershell
python example/mywork/TPMSMiner/C15/build_twin_combinations.py --cna
```

结果写入 `output/twin_combinations/`。当前组合包括全Σ3、仅Z16为Σ3、仅Z12
为Σ3、Z16取向交换、Z12取向交换和全Σ5负对照。具体定义及假设见
[`twin_combinations.json`](twin_combinations.json)，生成结果及每一对的实际
立方对称约化错配角见输出目录中的 `manifest.json`。
