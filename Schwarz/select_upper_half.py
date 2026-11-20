from pathlib import Path
import sys

import numpy as np
from ovito.io import export_file, import_file
from ovito.modifiers import CommonNeighborAnalysisModifier

folder = Path(__file__).resolve().parent
input_path = folder / (sys.argv[1] if len(sys.argv) > 1 else "1-1.6.data")
pipeline = import_file(str(input_path))
pipeline.modifiers.append(CommonNeighborAnalysisModifier())
data = pipeline.compute()
positions = data.particles.positions
cell_vectors = data.cell[:, :3]
cell_origin = data.cell[:, 3]
fractional = np.linalg.solve(cell_vectors, (positions - cell_origin).T).T
scaled = np.pi * fractional
sin_x, sin_y, sin_z = np.sin(scaled).T
cos_x, cos_y, cos_z = np.cos(scaled).T
s_d = (
    sin_x * sin_y * sin_z
    + sin_x * cos_y * cos_z
    + cos_x * sin_y * cos_z
    + cos_x * cos_y * sin_z
)
select = (np.abs(s_d) < 0.3).astype("int32")
data.particles_.create_property("select", data=select)

other = data.particles["Structure Type"] == int(CommonNeighborAnalysisModifier.Type.OTHER)
total_other = np.count_nonzero(other)
selected_other = np.count_nonzero(other & (select == 1))
ratio = f"{selected_other / total_other:.2%}" if total_other else "N/A"
print(f"Other (select=1) / Other = {selected_other} / {total_other} = {ratio}")

export_file(
    data,
    str(folder / f"{input_path.stem}-Select.data"),
    "lammps/dump",
    columns=["Particle Identifier", "Particle Type", "Position", "select", "Structure Type"],
    restricted_triclinic=True,
    precision=17,
)
