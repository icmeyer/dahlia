import numpy as np
import matplotlib.pyplot as plt
import copy
import sys; sys.path.append('..')
import pandas
from pyne import nucname

from main import main_function
from plotting import concentration_plot, concentration_subplot, flux_k_plot, \
                     convergence_plot, plot_k_comparison
from tools import print_isotopics_final, extract_PuAmCm, run_varying_enrichment, \
                  convert_to_weight_percent, convert_to_number_percent

# Flags for what to run
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
nat_conc[2] = 0.00711 # Weight percent
nat_conc[5] = (1-0.00711)

years = 50 / (35*1e-6*1e3) / 365.25 # MWd/kg / W/g * MW/W * g/kg / days/year
steps = int(1e3)

# Find thermal avg k
reactor_type = 'thermal'
conc_over_time, fluxes, ks = main_function(isotopes, conc, reactor_type, years, steps)
fresh_ks = ks
avg_k_fresh_thermal = np.average(ks)
PuAmCm_conc = extract_PuAmCm(isotopes, conc_over_time[:,-1])
print('Plutonium Vector')
print_isotopics_final(iso_names, conc, PuAmCm_conc)
csv_path = './'+reactor_type+'_fresh_concentrations.csv'
# flux_k_plot(years, steps, fluxes, ks)
print_isotopics_final(iso_names, conc, conc_over_time[:,-1], csv_path)
PuAmCm_thermals = [extract_PuAmCm(isotopes, conc_over_time[:,-1])]
print(reactor_type + 'Fresh Average k: ', avg_k_fresh_thermal)

if burnup_plot:
    title = '50 MWd/kg Burnup in ' + reactor_type + ' spectrum'
    concentration_subplot(iso_names, years, steps, conc_over_time, title)
    flux_k_plot(years, steps, fluxes, ks)

# Find fast avg k
reactor_type = 'fast'
conc_over_time, fluxes, ks = main_function(isotopes, conc, reactor_type, years, steps)
fresh_ks_fast = ks
avg_k_fresh_fast = np.average(ks)
csv_path = './'+reactor_type+'_fresh_concentrations.csv'
# flux_k_plot(years, steps, fluxes, ks)
print_isotopics_final(iso_names, conc, conc_over_time[:,-1], csv_path)
PuAmCm_fasts = [extract_PuAmCm(isotopes, conc_over_time[:,-1])]
print(reactor_type + 'Fresh Average k: ', avg_k_fresh_fast)

if burnup_plot:
    title = '50 MWd/kg Burnup in ' + reactor_type + ' spectrum'
    concentration_subplot(iso_names, years, steps, conc_over_time, title)
    flux_k_plot(years, steps, fluxes, ks)

### Finding the mix with the same average k 
if find_mix:
    avg_ks = []
    mix_ratios = np.linspace(0,1,15)

    # Once through
    reactor_type = 'thermal'
    avg_k_fresh = avg_k_fresh_thermal
    avg_ks = run_varying_enrichment(isotopes, iso_names, nat_conc, PuAmCm_conc,
                                    mix_ratios, reactor_type)
    once_through_mox_conc = np.interp(avg_k_fresh, avg_ks, mix_ratios)
    print('Pu+ Conc Required Once Through: ', once_through_mox_conc)
    plot_k_comparison(mix_ratios, avg_k_fresh, avg_ks, title='Fresh Ful: '+reactor_type)
    
    conc = nat_conc*(1-once_through_mox_conc) + PuAmCm_conc*(once_through_mox_conc)
    conc_over_time, fluxes, ks = main_function(isotopes, conc, reactor_type, years, steps)
    PuAmCm_conc = extract_PuAmCm(isotopes, conc_over_time[:,-1])
    PuAmCm_thermals.append(PuAmCm_conc)

    # Twice through
    reactor_type = 'thermal'
    avg_ks = run_varying_enrichment(isotopes, iso_names, nat_conc, PuAmCm_conc,
                                    mix_ratios, reactor_type)
    twice_through_mox_conc = np.interp(avg_k_fresh, avg_ks, mix_ratios)
    print('Pu+ Conc Required Twice Through: ', twice_through_mox_conc)
    plot_k_comparison(mix_ratios, avg_k_fresh, avg_ks, title='Once Through MOX: '+reactor_type)
    
    # conc = nat_conc*(1-once_through_mox_conc) + PuAmCm_conc*(once_through_mox_conc)
    # conc_over_time, fluxes, ks = main_function(isotopes, conc, reactor_type, years, steps)
    # PuAmCm_conc = extract_PuAmCm(isotopes, conc_over_time[:,-1])
    
    # Make data frame and write PuAmCm Concentrations to file
    PuAmCm_thermals = np.array(PuAmCm_thermals)
    data = np.vstack([np.array(iso_names), PuAmCm_thermals*100])
    data = data.transpose()
    headers = ['Isotope']
    stages = ['5% Enriched', 'Once Burned MOX']
    [headers.append(stage) for stage in stages]
    df = pandas.DataFrame(data=data, columns=headers)
    for stage in stages:
        df[stage] = pandas.to_numeric(df[stage]) 
    pandas.options.display.float_format = '{:04.2f}'.format
    for stage in stages:
        df.at['Total', stage] = df[stage].sum()
    df.at['Avg k', 'Isotope'] = avg_k_fresh
    df.at['Weight percent required', '5% Enriched'] = once_through_mox_conc
    df.at['Weight percent required', 'Once Burned MOX'] = twice_through_mox_conc
    csv_path = './'+reactor_type+'_plutonium_vectors_thermal_stages.csv'
    df.to_csv(csv_path)
    csv_path =  './'+reactor_type+'_plutonium_vectors_thermal_stages_clean.csv'
    df.to_csv(csv_path, float_format = '%05.2f')
    print(df)


