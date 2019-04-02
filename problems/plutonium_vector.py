import numpy as np
import matplotlib.pyplot as plt
import copy
import sys; sys.path.append('..')
import pandas
from pyne import nucname

from main import main_function
from plotting import concentration_plot, concentration_subplot, flux_k_plot, \
                     convergence_plot, plot_k_comparison
from tools import print_isotopics_final, extract_PuAmCm, run_varying_enrichment, convert_to_weight_percent

burnup_plot = False
calculate_post_burnup = False
find_mix = True

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
nat_conc[2] = 0.00711
nat_conc[5] = (1-0.00711)

years = 60 / (35*1e-6*1e3) / 365.25 # MWd/kg / W/g * MW/W * g/kg / days/year
steps = int(1e3)

burnups = [2,5,10,20,30,40,50,60,70]
pu_plus_concs = []
reactor_type = 'thermal'
for burnup in burnups:
    years = burnup /  (35*1e-6*1e3) / 365.25 # MWd/kg / W/g * MW/W * g/kg / days/year
    conc_over_time, fluxes, ks = main_function(isotopes, conc, reactor_type, years, steps)
    PuAmCm_conc = extract_PuAmCm(isotopes, conc_over_time[:,-1])
    # pu_plus_concs.append(PuAmCm_conc)
    PuAmCm_conc_weight = convert_to_weight_percent(isotopes, PuAmCm_conc)
    pu_plus_concs.append(PuAmCm_conc_weight)

pu_plus_concs = np.array(pu_plus_concs)
data = np.vstack([np.array(iso_names), pu_plus_concs*100])
data = data.transpose()
headers = ['Isotope']
[headers.append(str(burnup)) for burnup in burnups]
df = pandas.DataFrame(data=data, columns=headers)
for burnup in burnups:
    df[str(burnup)] = pandas.to_numeric(df[str(burnup)]) 

pandas.options.display.float_format = '{:04.2f}'.format
for burnup in burnups:
    df.at['Total', str(burnup)] = df[str(burnup)].sum()
csv_path = './'+reactor_type+'_plutonium_vector.csv'
df.to_csv(csv_path)
csv_path =  './'+reactor_type+'_plutonium_vector_clean.csv'
df.to_csv(csv_path, float_format = '%05.2f')
print(df)

# title = str(years) + ' years after burnup'
# concentration_subplot(iso_names, times[-1], steps, conc_over_time, title)
