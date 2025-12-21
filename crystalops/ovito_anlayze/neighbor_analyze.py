import os
import matplotlib
from ovito.io import import_file
import numpy as np
from collections import defaultdict
from ovito.modifiers import CoordinationAnalysisModifier, PolyhedralTemplateMatchingModifier
from ovito.data import CutoffNeighborFinder
from collections import defaultdict
import math
import pandas as pd

class NeighborAnalyzer:
    def __init__(self, 
                 dumpfile: str, 
                 cutoff_radius=3.5
                 ):
        """
        Automatically loads OVITO data and prepares neighbor information.
        """
        self.dumpfile = dumpfile
        self.cutoff_radius = cutoff_radius
        print(f"Initialized NeighborAnalyzer for file: {dumpfile}")

        # Initialize containers
        self.num_frames = None
        self.timesteps = []
        self.positions = []
        self.types = []
        self.ids = []
        self.structure_type = []
        self.coordination = []
        self.neighbors_dict = []
        self.neighbors_type = []

    def _build_neighbor_dict(self, data):
        """Build neighbor ID and type dictionaries for a single frame."""
        finder = CutoffNeighborFinder(self.cutoff_radius, data)
        frame_neighbors_dict = {}
        frame_neighbors_type = {}

        types = data.particles['Particle Type']
        ids = data.particles['Particle Identifier']

        for index in range(data.particles.count):
            atom_id = ids[index]
            frame_neighbors_dict[atom_id] = [ids[n.index] for n in finder.find(index)]
            frame_neighbors_type[atom_id] = [types[n.index] for n in finder.find(index)]

        # Store per-frame neighbor data
        self.neighbors_dict.append(frame_neighbors_dict)
        self.neighbors_type.append(frame_neighbors_type)

    def load_neighbor(self, 
                      outputfile: str = None
                      ):
        """
        Import file and compute structural data using OVITO modifiers.

        Args:
            outputfile: Optional filename prefix for saving.

        Returns:
            timesteps: The step of each nsnapshots
                shape: (nsnapshots)
            positions: The positions of the particles in numpy array
                shape: (nsnapshots, nparticles, 3)
            types: The atomic types of particles in numpy array
                shape: (nsnapshots, nparticles, 1)
            ids: The particle identifiers in numpy array
                shape: (nsnapshots, nparticles, 1)
            structure_type: The local structure types identified by PTM
                shape: (nsnapshots, nparticles, 1)
            coordination:
                shape: (nsnapshots, nparticles, 1)
            neighbors_dict: List of dictionaries
                shape: (nsnapshots, nparticles, variable)
            neighbors_type: List of dictionaries
                shape: (nsnapshots, nparticles, variable)
        """
        pipeline = import_file(self.dumpfile)
        pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff=self.cutoff_radius))
        ptm = PolyhedralTemplateMatchingModifier(output_orientation=True)
        pipeline.modifiers.append(ptm)

        nframes = pipeline.source.num_frames
        self.num_frames = nframes
        print(f"Loaded {self.dumpfile} with {nframes} frame(s).")

        for nframe in range(nframes):
            data = pipeline.compute(nframe)
            timestep = data.attributes['Ti  mestep']
            self.timesteps.append(timestep)  
            positions = data.particles.positions
            types = data.particles['Particle Type']  
            ids = data.particles['Particle Identifier']
            coordination = data.particles['Coordination']
            structure_type = data.particles['Structure Type']

            # Save per-frame data
            self.positions.append(positions)
            self.types.append(types)
            self.ids.append(ids)
            self.coordination.append(coordination)
            self.structure_type.append(structure_type)
            self._build_neighbor_dict(data)

            # Save to file if outputfile specified
            if outputfile:
                fname = f"{outputfile}_frame{nframe}.txt"
                with open(fname, "w") as f:
                    f.write("ID TYPE X Y Z CN PTM_TYPE NEIGHBOR_IDS NEIGHBOR_TYPES\n")
                    for i in range(len(ids)):
                        atom_id = ids[i]
                        atom_type = types[i]
                        x, y, z = positions[i]
                        cn = coordination[i]
                        ptm = structure_type[i]
                        neighbors = self.neighbors_dict[nframe][atom_id]
                        neighbor_types = self.neighbors_type[nframe][atom_id]
                        neighbor_str = ",".join(map(str, neighbors))
                        neighbor_type_str = ",".join(map(str, neighbor_types))
                        f.write(f"{atom_id} {atom_type} {x:.4f} {y:.4f} {z:.4f} {cn} {ptm} [{neighbor_str}] [{neighbor_type_str}]\n")
                print(f"Frame {nframe+1} saved to {fname}")

        return {
            "timesteps": self.timesteps,
            "positions": self.positions,
            "types": self.types,
            "ids": self.ids,
            "coordination": self.coordination,
            "structure_type": self.structure_type,
            "neighbors_dict": self.neighbors_dict,
            "neighbors_type": self.neighbors_type
        }
           

    def _computer_wc(self, ids, types, neighbors_types, composition):
        unique_types = np.unique(types)
        N_ij = {i: defaultdict(float) for i in unique_types}
        N_i = {i: 0.0 for i in unique_types}
        local_atoms = len(types)
        local_composition = {t: np.sum(types == t) / local_atoms for t in unique_types}

        # Loop over all atoms
        for atom_id in neighbors_types:
            atom_type = types[np.where(ids == atom_id)[0][0]]
            neigh_types = neighbors_types[atom_id]
            N_i[int(atom_type)] += len(neigh_types)
            for ntype in neigh_types:
                N_ij[int(atom_type)][int(ntype)] += 1
                
        # Average per atom of each type
        for i in unique_types:
            count_i = np.sum(types == i)
            N_i[i] /= count_i
            for j in unique_types:
                N_ij[i][j] /= count_i

        # Compute WC
        Wc = {i: {} for i in unique_types}
        for i in unique_types:
            for j in unique_types:
                if composition[j] > 0:
                    Wc[i][j] = 1 - N_ij[i][j] / (N_i[i] * composition[j])
                else:
                    Wc[i][j] = np.nan

        # Compute total SRO parameter: sum of (Wc_ii)^2
        SRO_WC = np.nansum([(Wc[i][i])**2 for i in unique_types if not np.isnan(Wc[i][i])])

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
                alpha[i] = {j: 0.0 for j in unique_types}
                continue
            for j in unique_types:
                val = N_ij[i][j] / N_i[i] if not np.isnan(N_ij[i][j]) else 0.0
                # numerical safety: clip small negatives to 0
                alpha[i][j] = max(0.0, float(val))
        
        # Shannon entropy S_i = - sum_j alpha_ij * ln(alpha_ij) (natural log)
        shannon_per_type = {}
        chem_param_per_type = {}  # exp(-S_i)
        SRO_shannon = 0.0
        for i in unique_types:
            S_i = -sum(p * math.log(p) for p in alpha[i].values() if p > 0)
            shannon_per_type[i] = S_i
            chem_param_per_type[i] = math.exp(-S_i)  # as requested
            SRO_shannon += chem_param_per_type[i] * composition[i]

        return {
            "Wc": Wc,
            "Wc_sym": Wc_sym,
            "SRO_WC": SRO_WC,
            "shannon_per_type": shannon_per_type,
            "chem_param_per_type": chem_param_per_type,
            "SRO_shannon": SRO_shannon,
            "local_composition": local_composition,
            "local_atoms": local_atoms
        }
    
    def _write_wc(self, outputfile, structure_labels, type_to_element, SRO, s=None):
        if s is None:
            state = 'w'
        else:
            state= 'a'
        with open(outputfile, state, encoding="utf-8") as f:          
            if s is None:
                f.write(f"\n=== Structure type all system ===\n")
            else:
                label = structure_labels.get(int(s), "All")
                f.write(f"\n=== Structure type {int(s)} ({label}) ===\n")
            if SRO is None:
                f.write(f"\nStructure type {int(s)} ({label}): No atoms found.\n")
                return

            f.write(f"Total atoms: {SRO['local_atoms']}\n")
            f.write("Composition (fraction): " + ", ".join(
            [f"{type_to_element.get(int(k), f'T{k}')}: {v:.4f}" for k, v in SRO["local_composition"].items()]) + "\n")
            f.write(f"SRO parameter (sum Wii^2): {SRO['SRO_WC']:.6f}\n\n")

            f.write("Directional Wc:\n")
            for i in sorted(SRO["Wc"].keys()):
                for j in sorted(SRO["Wc"][i].keys()):
                    f.write(f"Wc[{type_to_element.get(int(i))}->{type_to_element.get(int(j))}] = {SRO['Wc'][i][j]:.6f}\n")
            f.write("\n")

            f.write("Symmetric Wc:\n")
            for i in sorted(SRO["Wc_sym"].keys()):
                for j in sorted(SRO["Wc_sym"][i].keys()):
                    f.write(f"Wc_sym[{type_to_element.get(int(i))}-{type_to_element.get(int(j))}] = {SRO['Wc_sym'][i][j]:.6f}\n")
            f.write("\n")

            f.write("Shannon and exp(-S):\n")
            for i in sorted(SRO["shannon_per_type"].keys()):
                f.write(f"{type_to_element.get(int(i))}: S = {SRO['shannon_per_type'][i]:.6f}, exp(-S) = {SRO['chem_param_per_type'][i]:.6f}\n")
            f.write(f"Global chemical order parameter: {SRO['SRO_shannon']:.6f}\n\n")

    def Warren_Cowley(self,
                      outputfile: str = None,
                      type_to_element={1: "Ni", 2: "Co", 3: "Cr", 4: "Fe", 5: "Ni"},
                      sample: int = 1):
        """
        Compute Warren–Cowley short-range order (Wc) and related quantities
        for all frames and crystal structures.

        Args:
            outputfile: optional filename prefix for saving.
            type_to_element: mapping between type index and element symbol.
            sample: interval of frames to sample.

        Returns:
            results[frame][structure_label] = {...Wc results...}
        """
        if len(self.structure_type) == 0:
             print("No structure data found. Call load_neighbor() first.")
             self.load_neighbor()
        
        structure_labels = {0: "Other", 1: "FCC", 2: "HCP", 3: "BCC", 4: "ICO"}

        # Get unique structures and types from all frames
        unique_structures = np.unique(self.structure_type[0])
        unique_types = np.unique(self.types[0])

        total_atoms = len(self.types[0])
        composition = {t: np.sum(self.types[0] == t) / total_atoms for t in unique_types}
        print(f"Global composition: {composition}")

        results_list = []
        for frame in range(self.num_frames):
            if frame % sample != 0:
                continue

            print(f"\n=== Frame {frame} ===")
            frame_result = {'frame': frame}
            frame_result = {'step': self.timesteps[frame]} 
            frame_ids = self.ids[frame]
            frame_types = self.types[frame]
            frame_ids = np.array(frame_ids).flatten().astype(int)
            frame_neighbors_type = {atom_id: self.neighbors_type[frame][atom_id] for atom_id in frame_ids}

            sro = self._computer_wc(frame_ids, frame_types, frame_neighbors_type, composition)
            if outputfile:
                self._write_wc(f"{outputfile}_{frame}.txt",structure_labels, type_to_element, sro)

            frame_result['all_sro_wc'] = sro["SRO_WC"]
            frame_result['all_sro_shannon'] = sro["SRO_shannon"]
            for i in sorted(sro["Wc_sym"].keys()):
                for j in sorted(sro["Wc_sym"][i].keys()):
                    if j >= i:
                        frame_result[f'all_wc_[{type_to_element.get(int(i))}-{type_to_element.get(int(j))}]'] = sro['Wc_sym'][i][j]
                    else:
                        continue
            
            for s in unique_structures:
                mask = (self.structure_type[frame] == s)
                if np.sum(mask) == 0:
                    continue

                frame_ids = self.ids[frame][mask]
                frame_types = self.types[frame][mask]
                frame_neighbors_type = {
                    atom_id: self.neighbors_type[frame][atom_id] for atom_id in frame_ids
                }

                sro = self._computer_wc(frame_ids, frame_types, frame_neighbors_type, composition)
                if outputfile:
                    self._write_wc(f"{outputfile}_{frame}.txt", structure_labels, type_to_element, sro, s)
                frame_result[f'{structure_labels[s]}_sro_wc'] = sro["SRO_WC"]
                frame_result[f'{structure_labels[s]}_sro_shannon'] = sro["SRO_shannon"]
                for i in sorted(sro["Wc_sym"].keys()):
                    for j in sorted(sro["Wc_sym"][i].keys()):
                        if j >= i:
                            frame_result[f'{structure_labels[s]}_wc_[{type_to_element.get(int(i))}-{type_to_element.get(int(j))}]'] = sro['Wc_sym'][i][j]
                        else:
                            continue
            results_list.append(frame_result)
        results_df = pd.DataFrame(results_list)
        if outputfile:
            excel_file = f"{outputfile}_wc_results.xlsx"
            results_df.to_excel(excel_file, index=False)
        return results_df
    
    def structrue_fraction(self,
                           outputfile: str = None,
                           type_to_element={1: "Ni", 2: "Co", 3: "Cr", 4: "Fe", 5: "Ni"},
                      sample: int = 1):
        """
        Compute structure fraction and related elements fraction for all frames.

        Args:
            outputfile: optional filename prefix for saving.
            type_to_element: mapping between type index and element symbol.
            sample: interval of frames to sample.

        Returns:
            results_df: DataFrame containing structure fraction analysis
        """
        if len(self.structure_type) == 0:
             print("No structure data found. Call load_neighbor() first.")
             self.load_neighbor()

        structure_labels = {0: "Other", 1: "FCC", 2: "HCP", 3: "BCC", 4: "ICO"}

        # Get unique structures and types from all frames
        unique_structures = np.unique(self.structure_type[0])
        unique_types = np.unique(self.types[0])
        total_atoms = len(self.types[0])
        composition = {type_to_element.get(t, f"Type_{t}"): np.sum(self.types[0] == t) / total_atoms for t in unique_types}
        print(f"Global composition: {composition}")

        results_list = []
        for frame in range(self.num_frames):
            if frame % sample != 0:
                continue

            print(f"\n=== Frame {frame} ===")
            frame_result = {'frame': frame}
            frame_result = {'step': self.timesteps[frame]}

            # Calculate overall structure fractions
            for s in unique_structures:
                frame_result = {f'{structure_labels[s]}_all': np.sum(self.structure_type[frame] == s) / total_atoms}
                
                # Calculate element fractions within each structure type
                struct_mask = self.structure_type[frame] == s
                struct_total = np.sum(struct_mask)

                if struct_total > 0:  # Avoid division by zero
                    for t in unique_types:
                        # Count atoms of type t in structure s
                        count = np.sum((self.types[frame] == t) & struct_mask)
                        frame_result[f'{structure_labels[s]}_{type_to_element[t]}'] = count / struct_total
                else:
                    for t in unique_types:
                        frame_result[f'{structure_labels[s]}_{type_to_element[t]}'] = 0.0

            results_list.append(frame_result)
        results_df = pd.DataFrame(results_list)
        if outputfile:
            excel_file = f"{outputfile}_structure_fraction.xlsx"
            results_df.to_excel(excel_file, index=False)
            print(f"Results saved to {excel_file}")

        return results_df


import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from scipy.optimize import curve_fit

# ---- Gaussian function for curve fitting ----
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


def plot_Fe_CN_frame0_with_fit(analyzer, e_type=0, frame=0, bins=15, use_gaussian_fit=True, savefig=None):
    """
    绘制第一帧 Fe 配位数直方图 + 拟合曲线（Gaussian 或 KDE）
    """

    # --- Load data ---
    frame = 0
    types = np.array(analyzer.types[frame])
    cn = np.array(analyzer.coordination[frame])
    Fe_CN = cn[types == e_type]     # Fe type == 4

    if len(Fe_CN) == 0:
        print("No Fe atoms found.")
        return

    # --- Histogram data ---
    counts, bin_edges = np.histogram(Fe_CN, bins=bins, density=True)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    plt.figure(figsize=(6, 4))

    # --- Draw points like your image ---
    plt.plot(bin_centers, counts, "o", color="blue", label="Fe histogram")

    # =====================================================
    # Option A: Gaussian best-fit curve
    # =====================================================
    if use_gaussian_fit:
        # initial guess: peak height, mean, std
        p0 = [max(counts), np.mean(Fe_CN), np.std(Fe_CN)]
        popt, _ = curve_fit(gaussian, bin_centers, counts, p0=p0)

        x_fit = np.linspace(min(bin_edges), max(bin_edges), 200)
        y_fit = gaussian(x_fit, *popt)

        plt.plot(x_fit, y_fit, color="blue", linewidth=2,
                 label=f"Gaussian fit\nμ={popt[1]:.2f}, σ={popt[2]:.2f}")

    else:
        # =====================================================
        # Option B: KDE smoothing curve (non-parametric)
        # =====================================================
        kde = gaussian_kde(Fe_CN)
        x_fit = np.linspace(min(bin_edges), max(bin_edges), 300)
        y_kde = kde.evaluate(x_fit)
        plt.plot(x_fit, y_kde, color="blue", linewidth=2, label="KDE smoothing")

    # Formatting
    plt.xlabel("Coordination Number (CN)", fontsize=12)
    plt.ylabel("Probability Density", fontsize=12)
    plt.title("Fe Coordination Number Distribution (Frame 0)", fontsize=13)
    plt.legend()

    if savefig:
        plt.savefig(savefig, dpi=300, bbox_inches='tight')

    plt.show()

    return Fe_CN


    
# === Example usage ===
if __name__ == "__main__":
    analyzer = NeighborAnalyzer("Equli_1000_15.data", cutoff_radius=3.5)
    results = analyzer.load_neighbor()
    analyzer.Warren_Cowley('FiNiCrCoCu', sample=1, type_to_element={1: "Co", 2: "Cr", 3: "Cu", 4: "Fe", 5: "Ni"})