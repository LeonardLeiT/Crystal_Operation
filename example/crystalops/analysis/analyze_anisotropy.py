"""Run the existing anisotropy metrics through the analysis namespace."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from crystalops.analysis import aniso_alpha, aniso_lattice


def analyze_anisotropy(
    *,
    vector: tuple[float, float, float] = (1.0, 0.0, 0.0),
    lattice: tuple[float, float, float] = (1.0, 1.1, 1.2),
    reference_lattice: float = 1.0,
) -> dict[str, float]:
    """Evaluate the two legacy anisotropy functions without changing them."""
    return {
        "alpha": float(aniso_alpha(np.asarray(vector, dtype=float))),
        "lattice_deviation": float(aniso_lattice(lattice, a0=reference_lattice)),
    }


def main() -> None:
    """Read numerical inputs and print their anisotropy metrics."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vector", type=float, nargs=3, default=(1.0, 0.0, 0.0))
    parser.add_argument("--lattice", type=float, nargs=3, default=(1.0, 1.1, 1.2))
    parser.add_argument("--reference-lattice", type=float, default=1.0)
    options = parser.parse_args()

    result = analyze_anisotropy(
        vector=tuple(options.vector),
        lattice=tuple(options.lattice),
        reference_lattice=options.reference_lattice,
    )
    print(f"Alpha: {result['alpha']:.6f}")
    print(f"Lattice deviation: {result['lattice_deviation']:.6f}")


if __name__ == "__main__":
    main()
