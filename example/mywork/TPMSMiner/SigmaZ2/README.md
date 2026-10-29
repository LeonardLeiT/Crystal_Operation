# SigmaZ2：独立的 Z 双周期孪晶构型库

这里把 FK σ 相沿 Z 方向复制为共同参与三维周期分区的 60 个胞元，但孪晶面不复用 `../Sigma` 的设计：

- 拓扑固定为 `20 Z12 + 32 Z14 + 8 Z15`；
- 盒长约为 `120.13 × 120.13 × 124.41 Å`，每个构型约 15 万个 Cu 原子；
- 候选只来自 60-cell 周期 Voronoi/Laguerre 图中的真实有限公共面；
- 所有与原 Sigma 15 对孪晶相同的 modulo-30 源配对均被排除，包括跨层周期像；
- 同一构型内每个晶粒最多参加一对受控孪晶；未配对晶粒保持可复现的独立随机取向和独立晶格锚点。

生成默认的三分区均衡构型（15 对 Σ3）：

```powershell
python build_sigma_z2.py
```

主输出是 `output/01_sigma_z2_balanced15.data`，并写出 60 个胞元的位置、取向和配对信息到 `cell_manifest.json/csv`。

批量生成完整构型库：

```powershell
python build_twin_library.py --cna
```

`output/twin_library/` 包含 8 套相互侧重不同的 Σ3 构型，以及用高覆盖构型的同一批面构造的 Σ5 阴性对照。完整配额、胞元类型组合、实测消取向角和 CNA 统计见其中的 `manifest.json`。

该目录仍是候选构型实验，不接入 CrystalOps 公共 builder API。早期的 sparse6 试算包含原 Sigma 的源配对，现已从输出中移除，避免误用。
