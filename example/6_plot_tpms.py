from crystalops.users.tpms import TPMS_plot

tpms = TPMS_plot(30, 2)
tpms.plot_DSchwarz('TPMS', eq = 'DSC2')
tpms.plot_Gyroid('TPMS')
tpms.plot_IWP('TPMS')