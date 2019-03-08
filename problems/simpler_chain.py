import numpy as np
import sys; sys.path.append('..')
import pandas

from main import main_function
from plotting import concentration_plot

isotopes = ['922350', '922360', '922370', '922380', '922390']
conc = np.zeros(len(isotopes))
conc[0] = 0.05
conc[3] = 0.95
flux = 0
years = 5
steps = 10000
reactor_type = 'fast'
conc_over_time = main_function(isotopes, conc, flux, reactor_type, years, steps)

data = np.vstack([np.array(isotopes), conc, conc_over_time[:,-1]])
data = data.transpose()
headers = ['Isotope', 'Conc-init', 'Conc-final']

df = pandas.DataFrame(data=data, columns=headers)
print(df)


concentration_plot(isotopes, years, steps, conc_over_time)




