import numpy as np
import sys; sys.path.append('..')

from main import main_function
from plotting import concentration_plot

isotopes = ['922340', '922350', '922380', '952410', '952420', '952421']
conc = np.array([1, 1, 1])
flux = 1e24
years = 1
steps = 20
reactor_type = 'fast'
conc_over_time = main_function(isotopes, conc, flux, reactor_type, years, steps)
concentration_plot(isotopes, years, steps, conc_over_time)



