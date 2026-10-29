# Z foam：FK Z phase / HS-I clathrate dual

本目录实现 A15 之外的基础 Frank–Kasper Z topology。采用标准 `hP7-Zr4Al3`
原型，空间群 `P6/mmm`，每个六方原胞含：

- `3 × Z12`（Wyckoff `3f`）；
- `2 × Z14`（Wyckoff `2e`）；
- `2 × Z15`（Wyckoff `2d`）。

因此严格拓扑是 `3Z12 + 2Z14 + 2Z15`，不包含 Z16。HS-I 笼形水合物是该
FK 点集的拓扑对偶。

## 生成

```powershell
python build_z_foam.py
python build_twin_combinations.py --cna
```

基础输出有 7 个 LAMMPS data：完整 Z foam、位于原点的 Z12/Z14/Z15 代表胞元，
以及仅含 Z12、Z14、Z15 的全局构型。`cell_manifest.json/csv` 记录各胞元位置、
取向和孪晶伙伴；`validation.json` 记录拓扑、距离、等体积、CNA 和 LAMMPS
读取验证。

固定 `size=13.9`，目标约 7–8 万个 Cu 原子，不进行尺寸序列扫描。7 个胞元中
设置 3 对真实公共面 CSL：每类各一对；剩余的第三个 Z12 使用独立取向。

`output/twin_combinations/` 包含全 Σ3、单类 Σ3、Z14+Z15 Σ3 和全 Σ5 负对照。
这里生成的是几何可信的退火前驱体，是否形成 TPMS 必须由后续 LAMMPS 退火与
界面拓扑分析判断。

标准 120° 六方盒在 LAMMPS restricted-triclinic 表示中满足 `|xy|=lx/2`，旧版
LAMMPS 会给出 `Triclinic box skew is large` 性能提示；本构型已通过 `read_data`
和 `run 0`，该提示不是原子重叠或盒几何错误。
