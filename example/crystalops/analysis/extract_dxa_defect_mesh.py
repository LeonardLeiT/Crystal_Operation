"""Extract a DXA defect mesh from a complete atomistic structure."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from crystalops.io import read_structure


def extract_dxa_defect_mesh_file(
    filename: str | Path,
    *,
    output_dir: str | Path,
    crystal_structure: str = "fcc",
    surface_smoothing: int = 0,
    trial_circuit_length: int = 14,
    circuit_stretchability: int = 9,
    input_format: str | None = None,
    length_unit: str = "input_length_unit",
    prefix: str = "defect_mesh",
    plot: bool = False,
    show_atoms: bool = False,
    interactive: bool = True,
):
    """Read a structure, extract its DXA mesh, and save optional mesh views.

    The reference ``crystal_structure`` specifies the atomic lattice used by
    DXA. ``length_unit`` labels outputs without converting coordinates.
    Lengths and cell vectors retain the input coordinate unit. Plotting
    uses the display mesh clipped at periodic boundaries without adding caps.
    """
    from crystalops.analysis import extract_dxa_defect_mesh
    from crystalops.io import write_defect_mesh

    structure = read_structure(filename, format=input_format)
    structure.metadata["length_unit"] = length_unit
    result = extract_dxa_defect_mesh(
        structure,
        crystal_structure=crystal_structure,
        surface_smoothing=surface_smoothing,
        trial_circuit_length=trial_circuit_length,
        circuit_stretchability=circuit_stretchability,
    )
    paths = write_defect_mesh(result, output_dir, prefix=prefix)
    if plot or show_atoms:
        from crystalops.visualization import plot_defect_mesh
        import matplotlib.pyplot as plt

        figure, plot_paths = plot_defect_mesh(
            result,
            output_prefix=Path(output_dir) / prefix,
            structure=structure,
            show_atoms=show_atoms,
            interactive=interactive,
        )
        try:
            paths.update(plot_paths)
        finally:
            plt.close(figure)
    return result, paths


def main() -> None:
    """Parse DXA extraction settings, export the mesh, and report saved paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="complete structure, e.g. a LAMMPS atomic .data file")
    parser.add_argument("--input-format", default=None, help="explicit reader format, e.g. lammps-data")
    parser.add_argument(
        "--length-unit", default="input_length_unit",
        help="coordinate-unit label only; performs no unit conversion",
    )
    parser.add_argument(
        "--crystal-structure",
        choices=("fcc", "hcp", "bcc", "cubic_diamond", "hexagonal_diamond"),
        default="fcc",
        help="reference atomic lattice for DXA; diamond choices specify lattice types",
    )
    parser.add_argument("--surface-smoothing", type=int, default=0, help="DXA defect-mesh smoothing iterations")
    parser.add_argument("--trial-circuit-length", type=int, default=14)
    parser.add_argument("--circuit-stretchability", type=int, default=9)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--prefix", default="defect_mesh", help="prefix of exported mesh and view files")
    parser.add_argument("--plot", action="store_true", help="save a PNG and, by default, an offline interactive HTML view")
    parser.add_argument("--show-atoms", action="store_true", help="include atoms in the mesh views; implies --plot")
    parser.add_argument("--no-interactive", dest="interactive", action="store_false", help="save only the static PNG view")
    options = parser.parse_args()
    try:
        result, paths = extract_dxa_defect_mesh_file(
            options.input,
            output_dir=options.output_dir,
            crystal_structure=options.crystal_structure,
            surface_smoothing=options.surface_smoothing,
            trial_circuit_length=options.trial_circuit_length,
            circuit_stretchability=options.circuit_stretchability,
            input_format=options.input_format,
            length_unit=options.length_unit,
            prefix=options.prefix,
            plot=options.plot,
            show_atoms=options.show_atoms,
            interactive=options.interactive,
        )
    except (ValueError, ImportError, NotImplementedError, OSError) as error:
        parser.error(str(error))
    print(f"Periodic mesh: {len(result.vertices)} vertices, {len(result.faces)} triangles")
    print(
        f"Display mesh (without caps): {len(result.clipped_vertices)} vertices, "
        f"{len(result.clipped_faces)} triangles"
    )
    for name, path in paths.items():
        print(f"{name}: {Path(path).resolve()}")


if __name__ == "__main__":
    main()
