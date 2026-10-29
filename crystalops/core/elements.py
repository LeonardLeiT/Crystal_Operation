"""Small, dependency-free periodic-table helpers."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Element:
    """Describe the element metadata required by structure file formats.

    Args:
        symbol: IUPAC element symbol.
        atomic_number: Number of protons.
        atomic_mass: Standard atomic weight in atomic mass units.
    """

    symbol: str
    atomic_number: int
    atomic_mass: float


_ELEMENT_DATA: tuple[tuple[str, int, float], ...] = (
    ("H", 1, 1.008), ("He", 2, 4.003), ("Li", 3, 6.941),
    ("Be", 4, 9.012), ("B", 5, 10.811), ("C", 6, 12.011),
    ("N", 7, 14.007), ("O", 8, 15.999), ("F", 9, 18.998),
    ("Ne", 10, 20.180), ("Na", 11, 22.990), ("Mg", 12, 24.305),
    ("Al", 13, 26.982), ("Si", 14, 28.085), ("P", 15, 30.974),
    ("S", 16, 32.065), ("Cl", 17, 35.453), ("Ar", 18, 39.948),
    ("K", 19, 39.098), ("Ca", 20, 40.078), ("Sc", 21, 44.956),
    ("Ti", 22, 47.867), ("V", 23, 50.942), ("Cr", 24, 51.996),
    ("Mn", 25, 54.938), ("Fe", 26, 55.845), ("Co", 27, 58.933),
    ("Ni", 28, 58.693), ("Cu", 29, 63.546), ("Zn", 30, 65.380),
    ("Ga", 31, 69.723), ("Ge", 32, 72.630), ("As", 33, 74.922),
    ("Se", 34, 78.971), ("Br", 35, 79.904), ("Kr", 36, 83.798),
    ("Rb", 37, 85.468), ("Sr", 38, 87.620), ("Y", 39, 88.906),
    ("Zr", 40, 91.224), ("Nb", 41, 92.906), ("Mo", 42, 95.950),
    ("Tc", 43, 98.000), ("Ru", 44, 101.070), ("Rh", 45, 102.906),
    ("Pd", 46, 106.420), ("Ag", 47, 107.868), ("Cd", 48, 112.414),
    ("In", 49, 114.818), ("Sn", 50, 118.710), ("Sb", 51, 121.760),
    ("Te", 52, 127.600), ("I", 53, 126.904), ("Xe", 54, 131.293),
    ("Cs", 55, 132.905), ("Ba", 56, 137.327), ("La", 57, 138.905),
    ("Ce", 58, 140.116), ("Pr", 59, 140.908), ("Nd", 60, 144.242),
    ("Pm", 61, 145.000), ("Sm", 62, 150.360), ("Eu", 63, 151.964),
    ("Gd", 64, 157.250), ("Tb", 65, 158.925), ("Dy", 66, 162.500),
    ("Ho", 67, 164.930), ("Er", 68, 167.259), ("Tm", 69, 168.934),
    ("Yb", 70, 173.054), ("Lu", 71, 174.967), ("Hf", 72, 178.490),
    ("Ta", 73, 180.948), ("W", 74, 183.840), ("Re", 75, 186.207),
    ("Os", 76, 190.230), ("Ir", 77, 192.217), ("Pt", 78, 195.084),
    ("Au", 79, 196.967), ("Hg", 80, 200.592), ("Tl", 81, 204.383),
    ("Pb", 82, 207.200), ("Bi", 83, 208.980), ("Po", 84, 209.000),
    ("At", 85, 210.000), ("Rn", 86, 222.000), ("Fr", 87, 223.000),
    ("Ra", 88, 226.000), ("Ac", 89, 227.000), ("Th", 90, 232.038),
    ("Pa", 91, 231.036), ("U", 92, 238.029), ("Np", 93, 237.000),
    ("Pu", 94, 244.000), ("Am", 95, 243.000), ("Cm", 96, 247.000),
    ("Bk", 97, 247.000), ("Cf", 98, 251.000), ("Es", 99, 252.000),
    ("Fm", 100, 257.000), ("Md", 101, 258.000), ("No", 102, 259.000),
    ("Lr", 103, 262.000), ("Rf", 104, 267.000), ("Db", 105, 268.000),
    ("Sg", 106, 271.000), ("Bh", 107, 272.000), ("Hs", 108, 270.000),
    ("Mt", 109, 276.000), ("Ds", 110, 281.000), ("Rg", 111, 280.000),
    ("Cn", 112, 285.000), ("Nh", 113, 284.000), ("Fl", 114, 289.000),
    ("Mc", 115, 288.000), ("Lv", 116, 293.000), ("Ts", 117, 294.000),
    ("Og", 118, 294.000),
)

_ELEMENTS = {
    symbol: Element(symbol, atomic_number, atomic_mass)
    for symbol, atomic_number, atomic_mass in _ELEMENT_DATA
}


def get_element(symbol: str) -> Element:
    """Return element metadata for ``symbol``.

    Args:
        symbol: Case-sensitive element symbol such as ``"Cu"``.

    Returns:
        Element metadata.

    Raises:
        KeyError: If no metadata is available for the requested symbol.
    """
    try:
        return _ELEMENTS[symbol]
    except KeyError as exc:
        raise KeyError(f"unknown element symbol {symbol!r}") from exc


def infer_element(atomic_mass: float, *, tolerance: float = 0.2) -> str | None:
    """Infer an element symbol from an atomic mass when the match is unique."""
    matches = [
        element.symbol
        for element in _ELEMENTS.values()
        if abs(element.atomic_mass - atomic_mass) <= tolerance
    ]
    return matches[0] if len(matches) == 1 else None
