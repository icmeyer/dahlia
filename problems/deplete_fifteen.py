import numpy as np
import matplotlib.pyplot as plt
import copy
import sys; sys.path.append('..')
import pandas
from pyne import nucname

from main import main_function
from plotting import concentration_plot, concentration_subplot, flux_k_plot, \
                     convergence_plot, plot_k_comparison, concentration_subplot_bypercent
from tools import print_isotopics_final, extract_PuAmCm, run_varying_enrichment, \
                  convert_to_weight_percent, convert_to_number_percent

# Flags for what to run
burnup_plot = False
calculate_post_burnup = False
find_mix_thermal = True
find_mix_fast = True

isotopes = ['922330', '922340', '922350', '922360', '922370', '922380', '922390',
                      '932340', '932350', '932360', '932370', '932380', '932390', '932400',
                                '942350', '942360', '942370', '942380', '942390', '942400', '942410', '942420', '942430',
                                                                                  '952400', '952410', '952420', '952421', '952430', '952440',
                                                                                                      '962420', '962430', '962440', '962450']
# 92 U, 93 Np, 94 Pu, 95 Am, 96 Cm
iso_names = []
for i, iso in enumerate(isotopes):
    try:
        iso_names.append(nucname.serpent(iso))
    except:
        iso_names.append(iso)
iso_names.append('Lumped-Act.')
iso_names.append('FPs')

conc = np.zeros(len(isotopes)+2)
conc[2] = 0.05
conc[5] = 0.95

nat_conc = copy.deepcopy(conc)
nat_conc[2] = 0.00711 # Weight percent
nat_conc[5] = (1-0.00711)

years = 50 / (35*1e-6*1e3) / 365.25 # MWd/kg / W/g * MW/W * g/kg / days/year
steps = int(1e3)

# Make MOX
reactor_type = 'thermal'
conc_over_time, fluxes, ks = main_function(isotopes, conc, reactor_type, years, steps)
PuAmCm_conc = extract_PuAmCm(isotopes, conc_over_time[:,-1])

weight_percent = 0.30
mox_conc = weight_percent*PuAmCm_conc + (1-weight_percent)*nat_conc

# Deplete by 15 atom percent in fast reactor
years = 100 # Large number to hit 15 percent atom depletion
steps = int(1e4)
reactor_type = 'fast'
print("Start of 15 deplete function")
conc_over_time_fast, fluxes, ks = main_function(isotopes, mox_conc, reactor_type, years, steps, depletion_percent = 15)

title = 'Pu/Am/Cm of Fast Reactor Fuel'
concentration_subplot_bypercent(iso_names, steps, conc_over_time_fast, title, PuAmCm_only = True)

# Find fast avg k
reactor_type = 'thermal'
conc_over_time_thermal, fluxes, ks = main_function(isotopes, mox_conc, reactor_type, years, steps, depletion_percent = 15)

final_concs = []
final_concs.append(conc_over_time_fast[:,-1])
final_concs.append(conc_over_time_thermal[:,-1])

final_concs = np.array(final_concs)
data = np.vstack([np.array(iso_names), final_concs*100])
data = data.transpose()
headers = ['Isotope', 'Fast MOX', 'LWR MOX']
df = pandas.DataFrame(data=data, columns=headers)
df['Fast MOX'] = pandas.to_numeric(df['Fast MOX']) 
df.at['Total', 'Fast MOX'] = df['Fast MOX'].sum()
df['LWR MOX'] = pandas.to_numeric(df['LWR MOX']) 
df.at['Total', 'LWR MOX'] = df['LWR MOX'].sum()
pandas.options.display.float_format = '{:04.2f}'.format
csv_path = './'+reactor_type+'_fifteen_percent_deplete.csv'
df.to_csv(csv_path)
csv_path =  './'+reactor_type+'_fifteen_percent_deplete_clean.csv'
df.to_csv(csv_path, float_format = '%05.2f')
print(df)
