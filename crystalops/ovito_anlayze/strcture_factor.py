from ovito.io import import_file
from ovito.modifiers import CoordinationAnalysisModifier
import numpy as np
import matplotlib.pyplot as plt

class RDFAnalyzer:
    """
    RDFAnalyzer computes the Radial Distribution Function (RDF, g(r)) 
    and the Structure Factor (S(Q)) from atomic structure files 
    (e.g., LAMMPS dump files) using OVITO's CoordinationAnalysisModifier.
    
    Attributes:
        dumpfile (str): Path to the LAMMPS dump file.
        cutoff (float): Maximum distance for RDF calculation.
        number_of_bins (int): Number of bins to divide the distance range.
        r (np.ndarray): Array of distances corresponding to RDF.
        g_r (np.ndarray): RDF values at distances r.
        Q_values (np.ndarray): Array of Q values for structure factor.
        S_Q (np.ndarray): Computed structure factor S(Q).
    """
    def __init__(self, 
                 dumpfile: str, 
                 cutoff=20.0, 
                 number_of_bins=200
                 ):
        """
        Initialize the RDFAnalyzer.
        
        Args:
            dumpfile (str): Path to the LAMMPS dump file.
            cutoff (float, optional): Maximum distance for RDF. Default is 20.0 Å.
            number_of_bins (int, optional): Number of bins for RDF. Default is 100.
        """
        self.dumpfile = dumpfile
        self.cutoff = cutoff
        self.number_of_bins = number_of_bins
        self.r = None
        self.g_r = None
        self.Q_values = None
        self.S_Q = None

    def _save_rdf(self, 
                  outputfile: str = None
                  ):
        """Save RDF to a text file (UTF-8 encoding)."""
        with open(f'{outputfile}_RDF.txt', "w", encoding="utf-8") as f:
            np.savetxt(f, np.column_stack((self.r, self.g_r)),
                       header="r(Å)\t g(r)",
                       fmt="%.6f",
                       delimiter="\t")
        print(f"RDF results saved to: {outputfile}")

    def compute_rdf(self,
                    outputfile: str = None
                    ):
        """
        Compute the Radial Distribution Function (RDF) using OVITO.
        
        Args:
            outputfile (str, optional): If provided, save RDF to this file.

        Returns:
            list: [r, g_r] where r is the distance array and g_r is the RDF values.
        """
        pipeline = import_file(self.dumpfile)
        rdf_mod = CoordinationAnalysisModifier(cutoff=self.cutoff, number_of_bins=self.number_of_bins)
        pipeline.modifiers.append(rdf_mod)
        data = pipeline.compute()
        rdf_table = data.tables['coordination-rdf']
        self.r = rdf_table.xy()[:, 0]
        self.g_r = rdf_table.xy()[:, 1]
        print("RDF computed successfully.")
        if outputfile:
            self._save_rdf(outputfile)
        return [self.r, self.g_r]

    def _save_structure_factor(self, 
                               outputfile: str = None
                               ):
        """Save structure factor S(Q) to a text file (UTF-8 encoding)."""
        if self.Q_values is None or self.S_Q is None:
            raise ValueError("S(Q) not computed yet. Call compute_structure_factor() first.")
        with open(f'{outputfile}_SQ.txt', "w", encoding="utf-8") as f:
            np.savetxt(f, np.column_stack((self.Q_values, self.S_Q)),
                       header="Q(Å⁻¹)\t S(Q)",
                       fmt="%.6f",
                       delimiter="\t")
        print(f"Structure factor results saved to: {outputfile}")

    def _get_density_from_ovito(self):
        """Automatically extract number density ρ (atoms/Å³) from the LAMMPS dump."""
        pipeline = import_file(self.dumpfile)
        data = pipeline.compute()
        N = data.particles.count
        V = abs(data.cell.volume)
        rho = N / V
        print(f"Detected number density ρ = {rho:.5f} atoms/Å³ (N={N}, V={V:.1f} Å³)")
        return rho

    def compute_SQ(self, 
                   outputfile: str = None,
                   rho = None, 
                   Q_min=0.001, 
                   Q_max=15.0, 
                   n_Q=10000,):
        """
        Compute the Structure Factor S(Q) from the RDF using Fourier transform.
        
        Args:
            outputfile (str, optional): If provided, save S(Q) to this file.
            rho (float, optional): Atomic number density in Å^-3. Default is None.
            Q_min (float, optional): Minimum Q value in Å^-1. Default is 0.1.
            Q_max (float, optional): Maximum Q value in Å^-1. Default is 25.0.
            n_Q (int, optional): Number of Q points. Default is 1000.
        
        Returns:
            list: [Q_values, S_Q] where Q_values is the Q array and S_Q is the structure factor.
        """
        if self.r is None or self.g_r is None:
            print("RDF not computed yet. Computing RDF automatically...")
            self.compute_rdf()

        if rho is None:
            rho = self._get_density_from_ovito()

        self.Q_values = np.linspace(Q_min, Q_max, n_Q)
        self.S_Q = np.zeros_like(self.Q_values)
        for i, Q in enumerate(self.Q_values):
            integrand = (self.g_r - 1.0) * np.sin(Q * self.r) / (Q * self.r) * self.r**2
            self.S_Q[i] = 1 + 4 * np.pi * rho * np.trapz(integrand, self.r)
        print("Structure factor S(Q) computed successfully.")
        if outputfile:
            self._save_structure_factor(outputfile)
        return [self.Q_values, self.S_Q]



# === Example usage ===
if __name__ == "__main__":
    # structure = ['Equli_300_0']
    # for struc in structure:
    #     analyzer = RDFAnalyzer(f"{struc}.data", cutoff=20.0, number_of_bins=200)
    #     analyzer.compute_rdf(outputfile=struc)
    #     analyzer.compute_SQ(outputfile=struc, rho=0.05)
    structure = ['xtal', 'paracrystal-3', 'amorphous-1']
    structure = ['min_700_0']
    for struc in structure:
        analyzer = RDFAnalyzer(f"{struc}.data", cutoff=20.0, number_of_bins=200)
        analyzer.compute_rdf(outputfile=struc)
        analyzer.compute_SQ(outputfile=struc)