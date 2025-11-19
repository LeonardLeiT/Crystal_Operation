from ovito.io import import_file
import numpy as np
from collections import defaultdict
from ovito.modifiers import CoordinationAnalysisModifier, PolyhedralTemplateMatchingModifier
from ovito.data import CutoffNeighborFinder
from collections import defaultdict
import math

# 1. Read Data
input_file = "Lattice-min-9664.data"
output_file = f"{input_file}-neighbor.txt"
output_wc_file = f"{input_file}-wc.txt"

# 2. Prepare pipeline
cutoff_radius = 3.5
pipeline = import_file(input_file)
pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff=cutoff_radius))
ptm = PolyhedralTemplateMatchingModifier()
ptm.output_orientation = True
pipeline.modifiers.append(ptm)
data = pipeline.compute()

# 3. Extract data
positions = data.particles.positions
types = data.particles['Particle Type']
ids = data.particles['Particle Identifier']
coordination = data.particles['Coordination']
structure_type = data.particles['Structure Type']
orientation = data.particles['Orientation']

# 4. Build neighbor type dictionary
finder = CutoffNeighborFinder(cutoff_radius, data)
neighbors_dict = {}
neighbors_type = {}
for index in range(data.particles.count):
    atom_id = ids[index]
    neighbors_dict[atom_id] = [ids[n.index] for n in finder.find(index)]
    neighbors_type[atom_id] = [types[n.index] for n in finder.find(index)]

# 5. compute WC
def compute_wc(selected_mask, label):
    unique_types = np.unique(types)
    selected_ids = ids[selected_mask]
    selected_types = types[selected_mask]
    selected_neighbors_type = {atom_id: neighbors_type[atom_id] for atom_id in selected_ids}
    total_selected = len(selected_ids)
    print("Lable_type:", label)
    print(f"Total selected: {total_selected}")
    if total_selected == 0:
        return None
    
    composition = {t: np.sum(selected_types == t) / total_selected for t in unique_types}
    print(f"Composition: {composition}")

    N_ij = {i: defaultdict(float) for i in unique_types}
    N_i = {i: 0.0 for i in unique_types}
    
    # Loop over all atoms
    for atom_id in selected_neighbors_type:
        atom_type = selected_types[np.where(selected_ids == atom_id)[0][0]]
        neigh_types = selected_neighbors_type[atom_id]
        N_i[atom_type] += len(neigh_types)
        for ntype in neigh_types:
            N_ij[atom_type][ntype] += 1
            
    # Average per atom of each type
    for i in unique_types:
        count_i = np.sum(selected_types == i)
        N_i[i] /= count_i
        for j in unique_types:
            N_ij[i][j] /= count_i

    # Compute WC
    Wc = {}
    for i in unique_types:
        Wc[i] = {}
        for j in unique_types:
            if composition[j] > 0:
                Wc[i][j] = 1-N_ij[i][j] / (N_i[i] * composition[j])
            else:
                Wc[i][j] = np.nan

    print(f"\n=== Wc results for {label} ===")
    for i in unique_types:
        for j in unique_types:
            print(f"Wc[{i}->{j}] = {Wc[i][j]:.6f}")
            
    # Compute total SRO parameter: sum of (Wc_ii)^2
    SRO_value = np.nansum([(Wc[i][i])**2 for i in unique_types if not np.isnan(Wc[i][i])])
    print(f"SRO parameter for {label}: {SRO_value:.6f}")
    
    # Symmetric Wc_ij = 0.5*(Wc_ij + Wc_ji)
    Wc_sym = {}
    for i in unique_types:
        Wc_sym[i] = {}
        for j in unique_types:
            a = Wc[i].get(j, np.nan)
            b = Wc[j].get(i, np.nan)
            if np.isnan(a) and np.isnan(b):
                Wc_sym[i][j] = np.nan
            elif np.isnan(a):
                Wc_sym[i][j] = b
            elif np.isnan(b):
                Wc_sym[i][j] = a
            else:
                Wc_sym[i][j] = 0.5 * (a + b)
                
    # Compute alpha_{i->j} = fraction of neighbors of central type i that are type j
    # alpha = N_ij[i][j] / N_i[i]
    alpha = {}
    for i in unique_types:
        alpha[i] = {}
        if np.isnan(N_i[i]) or N_i[i] == 0:
            for j in unique_types:
                alpha[i][j] = 0.0
            continue
        for j in unique_types:
            val = N_ij[i][j] / N_i[i] if not np.isnan(N_ij[i][j]) else 0.0
            # numerical safety: clip small negatives to 0
            alpha[i][j] = max(0.0, float(val))

    # Shannon entropy S_i = - sum_j alpha_ij * ln(alpha_ij) (natural log)
    shannon_per_type = {}
    chem_param_per_type = {}  # exp(-S_i)
    chem_param_all = 0.0
    for i in unique_types:
        S_i = 0.0
        for j in unique_types:
            p = alpha[i][j]
            if p > 0:
                S_i -= p * math.log(p)   # natural log
        shannon_per_type[i] = S_i
        chem_param_per_type[i] = math.exp(-S_i)  # as requested
        chem_param_all += chem_param_per_type[i] * composition[i]
    
    return Wc, composition, total_selected, SRO_value, Wc_sym, shannon_per_type, chem_param_per_type, chem_param_all

# 6. Mapping for element types
# type_to_element = {1: "Ni", 2: "Co", 3: "Cr", 4: "Fe", 5: "Ni"}
# type_to_element = {1: "Fe", 2: "Ni", 3: "Cr", 4: "Co", 5: "Cu"}
type_to_element = {1: "Cu", 2: "Ta"}
structure_labels = {0: "Other", 1: "FCC", 2: "HCP", 3: "BCC", 4: "ICO"}

# 7. Save results for each structure
unique_structs = np.unique(structure_type)
print(f"\nDetected structure types: {unique_structs}")

with open(output_wc_file, "w", encoding="utf-8") as f:
    f.write("=== Wc, symmetric Wc, SRO and chemical-order parameters ===\n")
    f.write(f"Cutoff radius = {cutoff_radius:.2f} Å\n\n")

    # whole system
    f.write("=== Whole system ===\n")
    Wc_all, comp_all, total_all, SRO_all, Wc_sym_all, S_sh_all, chem_all, chem_param_all= compute_wc(np.ones(len(types), dtype=bool), "All atoms")
    f.write(f"Total atoms: {total_all}\n")
    f.write("Composition (fraction): " + ", ".join([f"{type_to_element.get(int(k), f'T{k}')}: {v:.4f}" for k, v in comp_all.items()]) + "\n")
    f.write(f"SRO parameter (sum Wii^2): {SRO_all:.6f}\n\n")

    # write directional Wc
    f.write("Directional Wc (Wc[i->j]):\n")
    for i in sorted(Wc_all.keys()):
        for j in sorted(Wc_all[i].keys()):
            f.write(f"Wc[{type_to_element.get(int(i),'T'+str(i))}->{type_to_element.get(int(j),'T'+str(j))}] = {Wc_all[i][j]:.6f}\n")
    f.write("\n")

    # write symmetric Wc
    f.write("Symmetric Wc (0.5*(Wc[i->j]+Wc[j->i])):\n")
    for i in sorted(Wc_sym_all.keys()):
        for j in sorted(Wc_sym_all[i].keys()):
            f.write(f"Wc_sym[{type_to_element.get(int(i),'T'+str(i))}-{type_to_element.get(int(j),'T'+str(j))}] = {Wc_sym_all[i][j]:.6f}\n")
    f.write("\n")

    # write Shannon and chem parameter per central type
    f.write("Per-type Shannon S_i and chemical-order parameter exp(-S_i):\n")
    for i in sorted(S_sh_all.keys()):
        f.write(f"{type_to_element.get(int(i),'T'+str(i))}: S_shannon = {S_sh_all[i]:.6f}, exp(-S) = {chem_all[i]:.6f}\n")
    f.write(f"chemical-order parameter in whole system: {chem_param_all:.6f}\n")
    f.write("\n\n")

    # per-structure blocks
    for s in unique_structs:
        mask = (structure_type == s)
        label = structure_labels.get(int(s), f"Unknown({int(s)})")
        Wc_result, comp, total, SRO_val, Wc_sym_res, S_sh_res, chem_res, chem_param_res = compute_wc(mask, label)
        if Wc_result is None:
            f.write(f"\nStructure type {int(s)} ({label}): No atoms found.\n")
            continue

        f.write(f"\n=== Structure type {int(s)} ({label}) ===\n")
        f.write(f"Total atoms: {total}\n")
        f.write("Composition (fraction): " + ", ".join([f"{type_to_element.get(int(k), f'T{k}')}: {v:.4f}" for k, v in comp.items()]) + "\n")
        f.write(f"SRO parameter (sum Wii^2): {SRO_val:.6f}\n\n")

        f.write("Directional Wc:\n")
        for i in sorted(Wc_result.keys()):
            for j in sorted(Wc_result[i].keys()):
                f.write(f"Wc[{type_to_element.get(int(i),'T'+str(i))}->{type_to_element.get(int(j),'T'+str(j))}] = {Wc_result[i][j]:.6f}\n")
        f.write("\n")

        f.write("Symmetric Wc:\n")
        for i in sorted(Wc_sym_res.keys()):
            for j in sorted(Wc_sym_res[i].keys()):
                f.write(f"Wc_sym[{type_to_element.get(int(i),'T'+str(i))}-{type_to_element.get(int(j),'T'+str(j))}] = {Wc_sym_res[i][j]:.6f}\n")
        f.write("\n")

        f.write("Shannon and exp(-S):\n")
        for i in sorted(S_sh_res.keys()):
            f.write(f"{type_to_element.get(int(i),'T'+str(i))}: S_shannon = {S_sh_res[i]:.6f}, exp(-S) = {chem_res[i]:.6f}\n")
        f.write(f"chemical-order parameter in whole system: {chem_param_res:.6f}\n")
        f.write("\n")

print(f"All Wc, symmetric Wc, SRO and chemical-order results saved to: {output_wc_file}")

# 6. Output file
with open(output_file, "w") as f:
    f.write("ID TYPE X Y Z CN PTM_TYPE NEIGHBOR_IDS\n")
    for i in range(len(ids)):
        atom_id = ids[i]
        atom_type = types[i]
        x, y, z = positions[i]
        cn = coordination[i]
        ptm = structure_type[i]
        neighbors = neighbors_dict[atom_id]
        neighbor_types = neighbors_type[atom_id]
        neighbor_str = ",".join(map(str, neighbors))
        neighbor_type_str = ",".join(map(str, neighbor_types))
        f.write(f"{atom_id} {atom_type} {x:.4f} {y:.4f} {z:.4f} {cn} {ptm} [{neighbor_str}] [{neighbor_type_str}]\n")

print(f"Done. Output saved to: {output_file}")