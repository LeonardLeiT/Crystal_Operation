"""Plot saved DXA defect-mesh triangles as a PNG and an offline interactive HTML."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))


def main() -> None:
    """Read extracted geometry without rerunning DXA or changing mesh topology."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path, help="saved defect-mesh NPZ file")
    parser.add_argument("--structure", type=Path, help="optional original structure for an atom overlay")
    parser.add_argument("--show-atoms", action="store_true", help="render every atom; requires --structure")
    parser.add_argument("--output-prefix", type=Path, required=True, help="prefix for .png and optional .html outputs")
    parser.add_argument("--no-interactive", action="store_true", help="save only PNG; Plotly is not required")
    options = parser.parse_args()
    if options.show_atoms and options.structure is None:
        parser.error("--show-atoms requires --structure")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from crystalops.io import read_defect_mesh, read_structure
    from crystalops.visualization import plot_defect_mesh

    try:
        result = read_defect_mesh(options.mesh)
        structure = read_structure(options.structure) if options.structure is not None else None
        figure, paths = plot_defect_mesh(
            result, output_prefix=options.output_prefix, structure=structure,
            show_atoms=options.show_atoms, interactive=not options.no_interactive,
        )
    except (ValueError, KeyError, OSError, ImportError) as error:
        parser.error(str(error))
    plt.close(figure)
    for kind, path in paths.items():
        print(f"{kind}: {path.resolve()}")


if __name__ == "__main__":
    main()
