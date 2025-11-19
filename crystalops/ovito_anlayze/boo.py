from ovito.io import import_file, export_file
from ovito.data import CutoffNeighborFinder
import numpy as np 
from scipy.special import sph_harm
import csv

class BOO:
    """
    calculates the bond orientational order (BOO) parameters (default l=6)
    for each atom in a trajectory using OVITO Basic.
    """
    def __init__(self, 
                 dumpfile: str, 
                 cutoff: float = 3.5, 
                 l: int = 6):
        """
        Args:
            dumpfile:
            cutoff: Spherical harmonics order for Steinhard parameter.
            l: Radial cutoff parameter for neighbor finding.
        """
        self.dumpfile = dumpfile
        self.cutoff = cutoff
        self.l = l
        self.neighbors_dict = []
        self.pipeline = import_file(self.dumpfile)
        self.n_frames = self.pipeline.source.num_frames
        print(f"Loaded {self.dumpfile} with {self.n_frames} frames.")

    def _compute_qlm(self, data):
        """
        Computes complex spherical harmonic components q_lm for all atoms in one frame.

        Args:
            data (DataCollection): OVITO data collection for one frame.

        Returns:
            tuple:
                qlm (np.ndarray): Array of shape (natoms, 2l+1) containing complex q_lm values.
                neighbor_list (dict): Dictionary {atom_id: [neighbor_ids]}.
        """
        # Total number of atoms.
        natoms = data.particles.count 
        ids = data.particles['Particle Identifier']
        # Compute total number of neighbors of each atom.
        finder = CutoffNeighborFinder(cutoff=self.cutoff, data_collection=data)
        
        # Neighbor count
        N_neigh = np.zeros(natoms, dtype=int)
        # Neigbor information
        neigh_offsets = []
        neighbor_list = {}

        # Gather neighbor vectors and counts
        for iatom in range(natoms):
            neigh_iter = finder.find(iatom)  # iterator
            # Consume iterator once and store both indices and offsets
            neighbors_info = [(n.index, np.array(n.delta)) for n in neigh_iter]
            
            atom_id = ids[iatom]
            neighbor_list[atom_id] = [idx for idx, _ in neighbors_info]
            N_neigh[iatom] = len(neighbors_info)
            neigh_offsets.extend([delta for _, delta in neighbors_info])

        neigh_offsets = np.array(neigh_offsets)

        # Distance to neighbors.
        d_ij = np.linalg.norm(neigh_offsets, axis=1)

        # phi, theta for spherical harmonics
        phi = np.arctan2(neigh_offsets[:, 1], neigh_offsets[:, 0])
        theta = np.arccos(neigh_offsets[:, 2] / d_ij)

        # Compute spherical harmonics
        Y = np.zeros((2 * self.l + 1, len(d_ij)), dtype=complex)
        for m in range(-self.l, self.l + 1):
            Y[m + self.l] = sph_harm(m, self.l, phi, theta)

        # Construct Steinhard vector by summing spherical harmonics.
        qlm = np.zeros((natoms, 2 * self.l + 1), dtype=complex)
        ineigh = 0
        for iatom in range(natoms):
            qlm[iatom] = np.sum(Y[:,ineigh:ineigh+N_neigh[iatom]],axis=1)
            ineigh += N_neigh[iatom]

        # Normalize by neighbor count
        qlm /= N_neigh[:, None]

        # Normalize each q_lm vector to unit length (like reference code)
        # q_norm = np.linalg.norm(qlm, axis=1)
        # qlm = (qlm.T / q_norm).T

        return qlm, neighbor_list

    def qlm(self):
        """
        Computes q_lm for all particles across all trajectory frames.

        Returns:
            qlm (np.ndarray):
                shape (nframes, natoms, 2l+1).
            neighbors_dict (list[dict]): List of neighbor dictionaries for each frame.
        """            
        qlm = []
        for frame in range(self.n_frames):
            data = self.pipeline.compute(frame)
            qlm_frame, neighbors_frame = self._compute_qlm(data)
            qlm.append(qlm_frame)
            self.neighbors_dict.append(neighbors_frame)
        qlm = np.array(qlm)
        print(f"Computed qlm {self.l} for {self.n_frames} frames, shape = {qlm.shape}")
        return qlm

    def ql(self, 
           outputfile: str = None
           ):
        """
        Computes the rotational invariant q_l for each atom in all frames.

        Args:
            outputfile: ovito dumpfile for saving.

        Returns:
            ql (np.ndarray):
                shape (nframes, natoms).
        """
        qlm = self.qlm()  # shape: (nframes, natoms, 2l+1)
        ql = np.sqrt(4 * np.pi / (2*self.l + 1) * np.sum(np.abs(qlm)**2, axis=2))
        print(f"Computed ql {self.l} for {self.n_frames} frames, shape = {ql.shape}")
        if outputfile:
            for frame in range(self.n_frames):
                data = self.pipeline.compute(frame)
                data.particles_.create_property(f'q_{self.l}', data=ql[frame])
                export_file(data,
                            f'{outputfile}_q{self.l}_{frame}.dump',
                            'lammps/dump',
                            columns=["Particle Identifier", "Particle Type", "Position.X",
                                    "Position.Y", "Position.Z", f'q_{self.l}'])
        return ql

    def qij(self, 
            outputfile: str = None
            ):
        """
        Computes the local bond orientational order correlation (qij) for all atoms across all frames.

        Args:
            outputfile: txt for saving.

        Returns:
            qij (list of list of np.ndarray): 
                Outer list over frames,
                inner list over atoms, 
                each np.ndarray contains qij values with neighbors.
            outputfile (str): CSV file path to save results.
        """ 
        # Compute qlm for all frames
        qlm_all = self.qlm()

        all_qij = []

        # Loop over frames
        for frame_idx in range(self.n_frames):
            data = self.pipeline.compute(frame_idx)
            qlm_frame = qlm_all[frame_idx]
            neighbor_list = self.neighbors_dict[frame_idx]

            frame_qij = []
            natoms = qlm_frame.shape[0]

            ids = np.array(data.particles['Particle Identifier'])

            for i in range(natoms):
                atom_id = ids[i]                        # get the actual atom ID
                neighbors = neighbor_list[atom_id]      # use atom ID as key
                q_i = qlm_frame[i]
                norm_i = np.sqrt(np.sum(np.abs(q_i)**2))

                qij_i = []
                for j in neighbors:
                    q_j = qlm_frame[j]
                    norm_j = np.sqrt(np.sum(np.abs(q_j)**2))
                    s_ij = np.sum(q_i * np.conj(q_j)) / (norm_i * norm_j)
                    qij_i.append(np.real(s_ij))        # store real part

                frame_qij.append(np.array(qij_i))

            all_qij.append(frame_qij)

        # Save to CSV if requested
        if outputfile:
            with open(f"{outputfile}.csv", mode='w', newline='') as f:
                writer = csv.writer(f)
                # Write header
                writer.writerow([
                    "Frame", "AtomID", "AtomType", "X", "Y", "Z",
                    "NeighborIDs", "qij"
                ])

                for frame_idx in range(self.n_frames):
                    data = self.pipeline.compute(frame_idx)
                    ids = np.array(data.particles['Particle Identifier'])
                    types = np.array(data.particles['Particle Type'])
                    positions = np.array(data.particles['Position'])

                    for i, neighbors in enumerate(all_qij[frame_idx]):
                        atom_id = ids[i]
                        atom_type = types[i]
                        pos = positions[i]

                        neighbor_ids = self.neighbors_dict[frame_idx][atom_id]
                        qij_vals = neighbors  # already np.array

                        writer.writerow([
                            frame_idx, atom_id, atom_type,
                            pos[0], pos[1], pos[2],
                            str(list(neighbor_ids)),
                            str(list(qij_vals))
                        ])
            print(f"qij results saved to {outputfile}")  
        print(f"Computed qij {self.l} for {self.n_frames} frames") 
        return all_qij
    
    def alpha(self, 
              outputfile: str = None,
              qq_cut: float = 0.8
              ):
        """ 
        Computes the α parameter (fraction of connected neighbors) for all atoms across all frames based on q_ij correlation.
        
        Args: 
            outputfile: ovito dumpfile for saving.
            qq_cut (float): cutoff for product of Steinhardt vectors.
                
        Returns: 
            alpha: (np.ndarray):
                shape (nframes, natoms).
        """
        qij = self.qij(outputfile=outputfile)
        alpha = []
        for frame_idx in range(self.n_frames):
            frame_qij = qij[frame_idx]
            natoms = len(frame_qij)
            alpha_frame = np.zeros(natoms)

            for i, qij_list in enumerate(frame_qij):
                if len(qij_list) == 0:
                    alpha_frame[i] = 0.0
                    continue
                # Fraction of neighbors with q_ij > qq_cut
                alpha_frame[i] = np.sum(np.array(qij_list) > qq_cut) / len(qij_list)
            alpha.append(alpha_frame)

        alpha = np.array(alpha)

        if outputfile:
            for frame in range(self.n_frames):
                data = self.pipeline.compute(frame)
                data.particles_.create_property('alpha', data=alpha[frame])
                export_file(data,
                            f'{outputfile}_alpha_{frame}.dump',
                            'lammps/dump',
                            columns=["Particle Identifier", "Particle Type", "Position.X",
                                    "Position.Y", "Position.Z", 'alpha'])
        print(f"Computed alpha {self.l} for {self.n_frames} frames, shape = {alpha.shape}")
        return alpha
    


# Example usage
if __name__ == "__main__":
    structure = ['xtal', 'amorphous-1', 'paracrystal-1']
    structure = ['min_500_0', 'min_700_0']
    for s in structure:
        boo = BOO(dumpfile=f'{s}.data', cutoff=3)
        # qlm = boo.qlm()
        ql = boo.ql(outputfile=s)
        # qij = boo.qij(outputfile=s)
        alpha = boo.alpha(outputfile=s, qq_cut=0.2)