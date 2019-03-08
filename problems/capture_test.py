import numpy as np
import sys; sys.path.append('..')

from main import main_function
from plotting import concentration_plot

isotopes = ['922390', '932390', '942390' ]
conc = np.array([1, 1, 1])
flux = 0
years = 2/365
steps = 80
reactor_type = 'fast'
conc_over_time = main_function(isotopes, conc, flux, reactor_type, years, steps)
concentration_plot(isotopes, years, steps, conc_over_time)



