import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))
from crystalops.users.tpms import TPMS_plot

tpms = TPMS_plot(20, 1)
tpms.plot_DSchwarz('TPMS', eq = 'DSC1')
tpms.plot_Gyroid('TPMS')
tpms.plot_IWP('TPMS')