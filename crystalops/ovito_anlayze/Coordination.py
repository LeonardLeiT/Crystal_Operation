import warnings
warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')

from ovito.io import import_file
from ovito.modifiers import CoordinationAnalysisModifier
import numpy as np
from collections import Counter
import pandas as pd

def Coordination_Distribution(dumpfile, 
                              element_type=None, 
                              cutoff=2.5, 
                              bins=200,
                              outputfile: str = None):
    """
    Compute Coordination Distribution of a single-frame dumpfile.

    Args:
        dumpfile: Path to atom configuration file
        element_type: The target element atom type. If None, all atoms are considered.
        cutoff: Cutoff distance to count coordination atoms.
        bins: Number of bins for coordination analysis (not used in Counter method).
        outputfile: Optional filename to save Excel output.

    Returns:
        DataFrame containing coordination distribution with columns:
        ['coord_vals', 'rate', 'count']
    """
    # 1. Reading file
    pipeline = import_file(dumpfile)

    # 2. CoordinationAnalysis
    coord_mod = CoordinationAnalysisModifier(
        cutoff=cutoff,
        number_of_bins=bins
    )
    pipeline.modifiers.append(coord_mod)

    data = pipeline.compute()
    coord_numbers = data.particles['Coordination']
    types = data.particles.particle_types

    # 3. Choose target element
    if element_type:
        mask = (types == element_type)
        type_coord = coord_numbers[mask]
    else:
        type_coord = coord_numbers

    if len(type_coord) == 0:
        raise RuntimeError(f"No atoms of type {element_type} found.")

    # 4. Counts
    count_dict = Counter(type_coord)
    coord_vals = np.array(sorted(count_dict.keys()))
    counts = np.array([count_dict[k] for k in coord_vals])

    # 5. Rates
    total_atoms = len(type_coord)
    rate = counts / total_atoms

    # 6. Build DataFrame
    results_df = pd.DataFrame({
        'coord_vals': coord_vals,
        'rate': rate,
        'count': counts
    })

    # 7. Save to Excel if requested
    if outputfile:
        excel_file = f"{outputfile}_Coordination_Distribution.xlsx"
        results_df.to_excel(excel_file, index=False)

    return results_df

# === Example usage ===
if __name__ == "__main__":
    filenames = ["Cool_500_0.data", "Cool_1000_0.data"]
    target_type = 2
    cutoff_distance = 2.5

    for fname in filenames:
        # Compute coordination distribution
        Coordination_Distribution(fname, element_type=target_type, cutoff=cutoff_distance, outputfile=fname)
        print(f'{fname}')
