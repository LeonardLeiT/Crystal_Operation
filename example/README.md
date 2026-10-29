# Examples and cases

除了保留原位置的 `mywork` 和 `mov`，其他示例按来源包和功能分类：

```text
example/
├── crystalops/
│   ├── core/          # 数据模型案例
│   ├── operations/    # 结构操作案例
│   ├── builders/      # 晶体生成器案例
│   └── io/            # LAMMPS 读写案例
├── crysplot/               # 绘图案例
├── mov/                    # 原项目动画脚本
└── mywork/                 # 原项目研究脚本
```

`crystalops` 下的 `test_*_case.py` 既是使用示例，也是可执行测试。

运行核心案例：

```bash
python -m pytest example/crystalops -q
```

运行绘图案例：

```bash
python -m example.crysplot.plot_case
```

直接运行结构生成案例：

```bash
python -m example.crystalops.builders.build_fcc --output fcc.data
python -m example.crystalops.builders.build_kelvin --seed 42 --output kelvin.data
python -m example.crystalops.builders.build_weaire_phelan \
  --twin-type sigma3 --output weaire_phelan_sigma3.data
python -m example.crystalops.builders.build_weaire_phelan \
  --twin-type sigma5 --output weaire_phelan_sigma5.data
python -m example.crystalops.builders.build_schwarz --output schwarz.data
python -m example.crystalops.builders.build_schwarz --type PSC \
  --size 8 --seed 42 --output schwarz_psc.data
python -m example.crystalops.builders.build_tpms --surface gyroid --output gyroid.data
python -m example.crystalops.io.modify_structure fcc.data --output fcc_modified.data
```

所有生成脚本都支持 `--help`，并且可以被导入而不会立即写文件。

Kelvin 构型始终生成 16 个完整的截角八面体晶粒，然后按 2×2×2 BCC
相对位置平移排列。每个晶粒只按自身的 6 个正方形面和 8 个六边形面
裁剪，不再用总模拟盒二次裁剪。预设只决定部分位置的固定取向，
其余位置使用由 `--seed` 控制的独立随机取向。如果需要在 OVITO 中
用 Particle Type 直接显示 16 个晶粒，可以使用：

模拟盒保持原始的 `2 × size × lattice_constant`。边界晶粒不做盒子裁剪，
因此完整晶粒中位于盒外的原子会保留，并交给 LAMMPS 周期边界处理。

```bash
python -m example.crystalops.builders.build_kelvin \
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

测试产生的临时结构文件由 pytest 自动清理。绘图案例会在当前目录生成
`structure_case.svg` 和 `tpms_case.svg`。
