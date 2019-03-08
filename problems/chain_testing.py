import numpy as np
import sys; sys.path.append('..')
import pandas

from main import main_function
from plotting import concentration_plot

# isotopes = ['922350', '922360', '922370', '922380', '922390',
#             '932390']
input_data = ['922380', 1,
              '922390', 0,
              '932390', 0]
isotopes = input_data[::2]
print(isotopes)
conc = np.array(input_data[1::2])
conc = conc/np.sum(conc)
print(conc)

flux = 0
years = 1
steps = 10000
reactor_type = 'fast'
conc_over_time = main_function(isotopes, conc, flux, reactor_type, years, 
                               steps, fixed_flux=True)

print(np.array(isotopes))
print(conc_over_time[:,-1])
print(np.array(isotopes).transpose())
print(conc_over_time[:,-1].transpose())
data = np.vstack([np.array(isotopes), conc, conc_over_time[:,-1]])
data = data.transpose()
print(data)
headers = ['Isotope', 'Conc-init', 'Conc-final']

df = pandas.DataFrame(data=data, columns=headers)
print(df)
# print('Total Concentration Sum: ',  ) 


concentration_plot(isotopes, years, steps, conc_over_time)




