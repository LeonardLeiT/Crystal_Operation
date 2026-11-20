# Examples and cases

除了保留原位置的 `mywork` 和 `mov`，其他示例与 `crystalops` 的模块职责对应：

```text
example/
├── crystalops/
│   ├── core/          # 数据模型案例
│   ├── modeling/
│   │   ├── builders/   # 晶体生成器案例
│   │   └── operations/ # 结构操作案例
│   ├── analysis/      # 结构及模拟结果分析案例
│   ├── visualization/ # 绘图案例
│   └── io/            # LAMMPS 读写案例
├── mov/                    # 原项目动画脚本
└── mywork/                 # 原项目研究脚本
```

`crystalops` 下的 `test_*_case.py` 既是使用示例，也是可执行测试。
旧的 `example.crystalops.builders` 和 `example.crysplot` 可执行示例保留为兼容入口；
原 builders、operations 测试文件已迁移，只在新的规范位置收集。

运行核心案例：

```bash
python -m pytest example/crystalops -q
```

运行绘图案例：

```bash
python -m example.crystalops.visualization.plot_case
```

运行已有的各向异性分析案例：

```bash
python -m example.crystalops.analysis.analyze_anisotropy
python -m example.crystalops.analysis.analyze_anisotropy --help
```

直接运行结构生成案例：

```bash
python -m example.crystalops.modeling.builders.build_fcc --output fcc.data
python -m example.crystalops.modeling.builders.build_kelvin --seed 42 --output kelvin.data
python -m example.crystalops.modeling.builders.build_weaire_phelan \
  --twin-type sigma3 --output weaire_phelan_sigma3.data
python -m example.crystalops.modeling.builders.build_weaire_phelan \
  --twin-type sigma5 --output weaire_phelan_sigma5.data
python -m example.crystalops.modeling.builders.build_schwarz --output schwarz.data
python -m example.crystalops.modeling.builders.build_schwarz --type PSC \
  --size 8 --seed 42 --output schwarz_psc.data
python -m example.crystalops.modeling.builders.build_schwarz --series
python -m example.crystalops.modeling.builders.build_tpms --surface gyroid --output gyroid.data
python -m example.crystalops.io.modify_structure fcc.data --output fcc_modified.data
```

所有生成脚本都支持 `--help`，并且可以被导入而不会立即写文件。

Schwarz 示例的 `--series` 生成六个周期 DSC 构型，盒长比例为
1:1:1、1:1:0.7、1:1:1.2、1:1:1.3、1:1:1.4 和 1:1:1.6。默认输出目录是
`analysis_results/dsc_twin_facets`。非立方盒中的 Σ3 孪晶面按各自的
盒长重建，并匹配两侧的 FCC {111} 堆垛。每个 `.data` 配有对应的
`.grain_ids.csv`，包含六区的 `domain_id` 和五种取向对应的 `grain_id`，
可传给晶界分析示例的 `--grain-labels`。用
`--type DSC --z-ratio 1.4 --output dsc_1-1-1.4.data` 可以单独生成一个比例。

非立方 DSC 可加 `--refine-boundaries` 整理重叠的过渡晶域与相邻晶域间的孤立片段。
该选项对 1.4 和 1.7 的初始构型做过验证，需要 SciPy（`pip install ".[analysis]"`）。

```bash
python -m example.crystalops.modeling.builders.build_schwarz --type DSC --z-ratio 1.7 --refine-boundaries --output analysis_results/dsc_boundary_refined/dsc_1-1-1.7.data
```

按到 DSC 控制曲面的距离给 `Schwarz` 目录的构型逐原子标记：

```bash
python -m example.crystalops.modeling.operations.select_schwarz_atoms Schwarz --behind 2 --ahead 2
```

示例沿用 DSC 构建器的四个局部零等值面及各自的旋转、平移和截取范围。
按周期最短距离选择曲面两侧的原子，`--behind` 与 `--ahead` 均用 Å。
默认两侧各 2 Å。该选区针对建模用的曲面片，平面孪晶界及弛豫后的实际
晶界需另行提取。计算需要 SciPy 和 scikit-image（`pip install ".[analysis]"`）。
输出位于 `analysis_results/schwarz_control_band`，保留全部原子，同时写入
OVITO 可直接识别的 `Selection` 列。

Kelvin 构型始终生成 16 个完整的截角八面体晶粒，然后按 2×2×2 BCC
相对位置平移排列。每个晶粒只按自身的 6 个正方形面和 8 个六边形面
裁剪，不再用总模拟盒二次裁剪。预设只决定部分位置的固定取向，
其余位置使用由 `--seed` 控制的独立随机取向。如果需要在 OVITO 中
用 Particle Type 直接显示 16 个晶粒，可以使用：

模拟盒保持原始的 `2 × size × lattice_constant`。边界晶粒不做盒子裁剪，
因此完整晶粒中位于盒外的原子会保留，并交给 LAMMPS 周期边界处理。

```bash
python -m example.crystalops.modeling.builders.build_kelvin \
  --label-grains --output kelvin_grains.data
```

核心案例覆盖：

- 创建 FCC/BCC 晶体与超胞。
- LAMMPS data 文件写出和读回。
- 平移、旋转、周期包裹、选区和原子类型设置。
- Kelvin 多晶的随机种子复现。
- Weaire-Phelan 的论文 Σ3 构型和 Σ5 对照构型。
- Schwarz D-SCI 六个构型区域的完整性。
- Gyroid TPMS 晶体筛选。

`analysis` 提供各向异性和晶界平均曲率分析示例。晶界计算目前要求正交模拟盒、
三个方向均周期，以及完整的当前晶粒编号。生成器中的 `grain_ids` 保存在内存
结构的逐原子属性中，当前 `atomic` data 读写不会保留这些编号；文件分析可通过
`--grain-labels` 载入按原子 ID 对齐的 CSV，或显式使用 `--segment-grains` 恢复编号。

例如，对当前目录的 Cu 结构计算一个指定尺度下的晶界曲率：

```bash
python -m example.crystalops.analysis.analyze_grain_boundary_curvature PSC.data \
  --segment-grains --sigma 2.0 --spacing 0.7 --neighbor-cutoff 3.3 \
  --output-dir analysis_results/PSC_sigma2
python -m example.crystalops.analysis.analyze_grain_boundary_curvature --help
```

长度参数与输入坐标使用相同单位；所附 Cu 结构为 Å，曲率为 Å⁻¹。输出包括
原子/曲面 CSV、带曲率和有效性标记的 VTK 网格、使用的晶粒编号与 JSON 统计。
需要分别检查固定平滑尺度下的采样收敛和平滑尺度敏感性；方法和输出说明见
`crystalops/analysis/README.md`。

测试产生的临时结构文件由 pytest 自动清理。绘图案例会在当前目录生成
`structure_case.svg` 和 `tpms_case.svg`。

提取完整原子结构的 DXA defect mesh：

```bash
pip install -e ".[defect-mesh]"
python -m example.crystalops.analysis.extract_dxa_defect_mesh 1-1.data \
  --crystal-structure fcc --length-unit angstrom --output-dir analysis_results/1-1_defect_mesh
```

`--surface-smoothing` 默认是 `0`。晶体类型可选 `fcc`、`hcp`、`bcc`、
`cubic_diamond`、`hexagonal_diamond`，其中 diamond 表示原子晶格类型。
输出包括原始周期网格 NPZ、用于显示的不封口 VTK 网格及 JSON 元数据。
`--length-unit` 只标注原坐标单位，不执行换算；默认 `input_length_unit`，
上述 Cu 示例显式标注为 angstrom。

需要 PNG 和可离线旋转、缩放的 HTML 视图时：

```bash
pip install -e ".[defect-mesh,mesh-viewer]"
python -m example.crystalops.analysis.extract_dxa_defect_mesh 1-1.data \
  --crystal-structure fcc --length-unit angstrom --output-dir analysis_results/1-1_defect_mesh --plot
```

`--show-atoms` 可叠加原子并自动启用绘图，`--no-interactive` 只保存 PNG。

批量提取并对比四个构型：

```bash
python -m example.crystalops.analysis.extract_dxa_defect_mesh_batch \
  1-1.6.data 1-1.data 1.2-1.data 1.4-1.data \
  --crystal-structure fcc --surface-smoothing 0 --length-unit angstrom \
  --output-dir analysis_results/dxa_defect_mesh
```

每个构型写入以文件 stem 命名的子目录，批量入口总是保存 PNG 和离线交互 HTML。
输出根目录的 `index.html` 以 2×2 卡片嵌入四个可旋转视图，并显示网格三角数、
盒矢量长度和平滑次数；也可单独打开各 HTML 或 PNG。可加 `--show-atoms`。
