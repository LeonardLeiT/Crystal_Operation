# Frank–Kasper σ foam

四方 `P4₂/mnm` σ相包含30个周期胞元：`10×Z12 + 16×Z14 + 4×Z15`。
本目录采用五组Wyckoff坐标和等体积Laguerre权重构造完整胞元，再沿真实
Voronoi公共面布置共享相位的Σ3孪晶。

```powershell
python example/mywork/TPMSMiner/Sigma/build_sigma.py
```

尺寸固定为 `size=8.6`，不提供尺寸序列。Cu晶格常数3.61 Å时目标为7–8万
原子。`output/`输出7个LAMMPS构型：完整σ foam、位于原点的代表性Z12、
Z14、Z15胞元，以及仅Z12、仅Z14、仅Z15的三种构型。30个胞元的位置、
Wyckoff类型、取向和孪晶关系记录在`cell_manifest.json/csv`。

## 孪晶组合

```powershell
python example/mywork/TPMSMiner/Sigma/build_twin_combinations.py --cna
```

固定尺寸组合输出到`output/twin_combinations/`：全Σ3、仅Z12–Z12为Σ3、
仅Z14–Z14为Σ3、仅Z15–Z15为Σ3、仅Z12–Z14桥接面为Σ3，以及全Σ5
负对照。组合定义见[`twin_combinations.json`](twin_combinations.json)。
