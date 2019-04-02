import numpy as np
import pandas
import copy

import main

def run_varying_enrichment(isotopes, iso_names, nat_conc, PuAmCm_conc, mix_ratios, reactor_type):
    avg_ks = []
    for ratio_mox in mix_ratios:
        conc = nat_conc*(1-ratio_mox) + PuAmCm_conc*(ratio_mox)

        years = 50 / (35*1e-6*1e3) / 365.25 # MWd/kg / W/g * MW/W * g/kg / days/year
        steps = int(1e3)
        conc_over_time, fluxes, ks = main.main_function(isotopes, conc, reactor_type, years, steps)

        print('mox isotopics')
        print_isotopics_final(iso_names, conc, conc_over_time[:,-1])
        avg_ks.append(np.average(ks))
    return avg_ks


def extract_PuAmCm(isotopes, final_conc):
    index_of_Pu_start = isotopes.index('942350')
    index_of_Cm_end   = isotopes.index('962450')
    PuAmCm_conc = np.zeros_like(final_conc)
    PuAmCm_only_conc = final_conc[index_of_Pu_start:index_of_Cm_end+1]
    PuAmCm_only_conc_norm = PuAmCm_only_conc/np.sum(PuAmCm_only_conc)
    PuAmCm_conc[index_of_Pu_start:index_of_Cm_end+1] = PuAmCm_only_conc_norm
    return PuAmCm_conc

def print_isotopics_final(iso_names, conc_initial, conc_final, csv_path=None):
    data = np.vstack([np.array(iso_names), conc_initial*100, conc_final*100])
    data = data.transpose()
    headers = ['Isotope', 'Conc-init', 'Conc-final']
    
    pandas.options.display.float_format = '{:2.4f}'.format
    df = pandas.DataFrame(data=data, columns=headers)
    df['Conc-init'] = pandas.to_numeric(df['Conc-init'])
    df['Conc-final'] = pandas.to_numeric(df['Conc-final'])
    df.loc['Total'] = pandas.Series(df['Conc-final'].sum(), index = ['Conc-final'])
    df.at['Total', 'Conc-init'] = df['Conc-init'].sum()
    if csv_path:
        df.to_csv(csv_path)
    print(df)

def convert_to_weight_percent(isotopes, conc):
    mass_concs = np.zeros(len(isotopes))
    for i, iso in enumerate(isotopes):
        mass = np.float(iso[2:5])
        print('mass', mass, i, iso)
        mass_concs[i] = conc[i]*mass
    mass_concs_norm = mass_concs/np.sum(mass_concs)
    conc_return = copy.deepcopy(conc)
    conc_return[0:len(isotopes)] = mass_concs_norm
    conc_return_norm =  conc_return/np.sum(conc_return)
    return conc_return_norm

def convert_to_number_percent(isotopes, conc):
    number_concs = np.zeros(len(isotopes))
    for i, iso in enumerate(isotopes):
        mass = np.float(iso[2:5])
        number_concs[i] = conc[i]/mass
    number_concs_norm = number_concs/np.sum(number_concs)
    conc_return = copy.deepcopy(conc)
    conc_return[0:len(isotopes)] = number_concs_norm
    conc_return_norm =  conc_return/np.sum(conc_return)
    return conc_return_norm
