import numpy as np
import sys; sys.path.append('..')
import pandas
from pyne import nucname

from main import main_function
from plotting import concentration_plot, concentration_subplot, flux_k_plot, \
                     convergence_plot

burnup_plot = True
estimate_convergence = False
calculate_post_burnup = True


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
flux = 0
years = 60 / (35*1e-6*1e3) / 365.25 # MWd/kg / W/g * MW/W * g/kg / days/year
steps = int(1e3)
reactor_type = 'fast'
conc_over_time, fluxes, ks = main_function(isotopes, conc, flux, reactor_type, years, steps)

data = np.vstack([np.array(iso_names), conc*100, conc_over_time[:,-1]*100])
data = data.transpose()
headers = ['Isotope', 'Conc-init', 'Conc-final']

pandas.options.display.float_format = '{:2.4f}'.format
df = pandas.DataFrame(data=data, columns=headers)
df['Conc-init'] = pandas.to_numeric(df['Conc-init'])
df['Conc-final'] = pandas.to_numeric(df['Conc-final'])
df.loc['Total'] = pandas.Series(df['Conc-final'].sum(), index = ['Conc-final'])
df.at['Total', 'Conc-init'] = df['Conc-init'].sum()
csv_path = './'+reactor_type+'_concentrations.csv'
df.to_csv(csv_path)
print(df)
print(df['Conc-init'].sum())

if burnup_plot:
    title = '60 MWd/kg Burnup in ' + reactor_type + ' spectrum'
    concentration_subplot(iso_names, years, steps, conc_over_time, title)
    flux_k_plot(years, steps, fluxes, ks)


if estimate_convergence:
    # Check temporal convergence
    final_concs = []
    time_steps = np.floor(np.logspace(0,3,10))
    print(time_steps)
    for steps in time_steps:
        conc_over_time, fluxes, ks = main_function(isotopes, conc, flux, reactor_type, years, int(steps))
        final_concs.append(conc_over_time[:,-1])
    
    convergence_plot(iso_names, time_steps, final_concs)

if calculate_post_burnup:
    burnup_conc = conc_over_time[:,-1]
    times = np.logspace(-1, 4, 6)
    final_concs = np.empty([len(times), len(isotopes)+2])
    for i, time in enumerate(times):
        # Calculate decay of isotopes after removal
        flux = 0
        years = time # MWd/kg / W/g * MW/W * g/kg / days/year
        steps = int(1e2)
        conc_over_time, fluxes, ks = main_function(isotopes, burnup_conc, flux, reactor_type, years, steps, fixed_flux=True)
        final_concs[i,:]  = conc_over_time[:,-1]
    
    print(burnup_conc)
    print(final_concs)
    data = np.vstack([np.array(iso_names), burnup_conc*100, final_concs*100])
    data = data.transpose()
    headers = ['Isotope', 'Burnup-Conc']
    [headers.append(str(time)) for time in times]
    
    df = pandas.DataFrame(data=data, columns=headers)
    df['Burnup-Conc'] = pandas.to_numeric(df['Burnup-Conc'])
    for time in times:
        df[str(time)] = pandas.to_numeric(df[str(time)]) 

    pandas.options.display.float_format = '{:2.4f}'.format
    df.loc['Total'] = pandas.Series(df['Burnup-Conc'].sum(), index = ['Burnup-Conc'])
    for time in times:
        df.at['Total', str(time)] = df[str(time)].sum()
    csv_path = './'+reactor_type+'_concentrations_post_burnup.csv'
    df.to_csv(csv_path)
    print(df)
    
    title = str(years) + ' years after burnup'
    concentration_subplot(iso_names, times[-1], steps, conc_over_time, title)
