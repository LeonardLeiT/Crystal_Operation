import numpy as np
from crystalops.config.lattice import *

lp = Lattice_plot(outfile="lattice_plot", type_fig="svg")

print("Testing Body-Centered Cubic...")
lp.BCC()

print("Testing Body-Centered Tetragonal...")
lp.BCT(c=1.3)

print("Testing Body-Centered Orthorhombic...")
lp.BCO(a=1, b=1.1, c=1.3)

print("Testing Face-Centered Cubic...")
lp.FCC()

print("Testing Face-Centered Orthorhombic...")
lp.FCO(a=1.1, b=1.2, c=1.3)

print("Testing Simple Cubic...")
lp.SC()

print("Testing Simple Tetragonal...")
lp.ST(c=1.3)

print("Testing Simple Orthorhombic...")
lp.SO(a=1.1, b=1.2, c=1.3)

print("Testing Simple Monoclinic...")
lp.SM(a=1.1, b=1.2, c=1.3, alpha_deg=80)

print("Testing Triclinic...")
lp.TRI(lengths=[1, 1.2, 1.4], angles=[80, 85, 75])

print("Testing Rhombohedral...")
lp.RHL(alpha_deg=75)

print("Testing Base-Centered Orthorhombic...")
lp.ABO(a=1.1, b=1.2, c=1.3)

print("Testing Base-Centered Monoclinic...")
lp.BCM(a=1.1, b=1.2, c=1.3, alpha_deg=100)

print("Testing Hexagonal Close Packed...")
lp.HCP(c=1.633)

print("All lattice types tested.")

