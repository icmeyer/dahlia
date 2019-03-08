import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.constants import physical_constants

# from solver import solve_matrix
from plotting import concentration_plot

from data import FAST_DATA, THERMAL_DATA, YR2S, INV_YR2S, MOL, EV2JOULE


use_pyne = True
if use_pyne:
    from pyne import data, nucname

# Note: Isotopes must be listed by ascending Z and then ascending A
def main_function(isotopes, conc, flux, reactor_type, years, steps,
                  fixed_flux=False):
    if reactor_type=='fast':
        DATA = FAST_DATA
    elif reactor_type=='thermal':
        DATA = THERMAL_DATA

    fluxes = []
    flux = evaluate_flux(isotopes, conc, DATA, fixed_flux)
    ks = []
    ks.append(evaluate_kinf(isotopes, conc, DATA, reactor_type))
    depletion_matrix = np.zeros([len(isotopes)+2, len(isotopes)+2])
    
    flux_depletion = calc_flux_depletion(isotopes, flux, DATA)
    print('Flux depletion matrix')
    print(flux_depletion)
    # plt.spy(flux_depletion)
    # plt.title('flux depletion')
    # plt.show()
    
    decay_depletion = calc_decay_depletion(isotopes, DATA)
    print('Decay depletion matrix')
    print(decay_depletion)
    # plt.spy(decay_depletion)
    # plt.title('decay depletion')
    # plt.show()

    # print('By half-life')
    # print(np.log(2)/decay_depletion/3600)

    depletion_matrix = flux_depletion + decay_depletion
    # plt.spy(depletion_matrix)
    # plt.show()

    time = YR2S*years # Convert to seconds
    dt = time/steps
    conc_over_time = np.zeros([len(isotopes)+2, steps])
    conc_over_time[:, 0] = conc
    for i in range(steps-1):
        conc = matrix_exp_solve(depletion_matrix, dt, conc)
        conc_over_time[:, i+1] = conc
        # print('flux evaluation')
        flux = evaluate_flux(isotopes, conc, DATA, fixed_flux)
        fluxes.append(flux)
        ks.append(evaluate_kinf(isotopes, conc, DATA, reactor_type))
        depletion_matrix = calc_flux_depletion(isotopes, flux, DATA) + decay_depletion

    return conc_over_time, fluxes, ks

def calc_decay_depletion(isotopes, DATA):
    # Create dictionary to keep track of daughters
    daughter_flags = dict()
    for iso in isotopes:
        daughter_flags[iso] = False

    dec_depletion = np.zeros([len(isotopes)+2, len(isotopes)+2])
    # Decay from j to i
    for i, iso_i in enumerate(isotopes):
        for j, iso_j in enumerate(isotopes):
            # Off-diagonal terms
            if j < i or j > i:
                ratio = data.branch_ratio(iso_j, iso_i)
                # print('print(iso_j, iso_i, ratio)')
                # print(iso_j, iso_i, ratio)
                # print(iso_j, iso_i, ratio)
                if ratio > 0:
                    decay_const = data.decay_const(iso_j)
                    dec_depletion[i,j] -= ratio*decay_const 
                    daughter_flags[iso_j] = True

                if iso_j=='932340':
                    print(iso_j, iso_i, data.branch_ratio(iso_j,iso_i))

            # Diagonal terms
            elif i == j:
                # Add decay constant here
                # print(iso_j, data.half_life(iso_j))
                dec_depletion[i,j] +=  data.decay_const(iso_j)

    for i, iso in enumerate(isotopes):
        if not daughter_flags[iso] and dec_depletion[i,i] != 0:
            print("Warning: no decay daughter for ", iso)
            dec_depletion[-2,i] += -dec_depletion[i,i]

    return dec_depletion


def calc_flux_depletion(isotopes, flux, DATA):
    flux_depletion = np.zeros([len(isotopes)+2, len(isotopes)+2])
    # Capture from j to i
    for i, iso_i in enumerate(isotopes):
        for j, iso_j in enumerate(isotopes):
            if j <= i:
                # Neutron Capture
                diff = float(iso_i) - float(iso_j)
                one_A_diff = [9, 10, 11]
                # print(diff)
                # print(iso_i, iso_j, i, j)
                if diff in one_A_diff:
                    try:
                        cap_xs =  DATA[iso_j]['cap']*1e-24
                    except:
                        cap_xs = 0
                    # Americium meta-stable cases
                    # print(iso_i, iso_j, cap_xs)
                    if iso_i == '952420' and iso_j == '952410':
                        cap_xs = 0.11*cap_xs
                    elif iso_i == '952421' and iso_j == '952410':
                        cap_xs = 0.89*cap_xs
                    # print("omg here!!!")
                    # print('cap from data',DATA[iso_j]['cap'])
                    # print('cap', cap_xs)
                    # print(flux*cap_xs)
                    flux_depletion[i,j] -= flux*cap_xs

                if i == j:
                    try:
                        fis_xs = DATA[iso_i]['fis']*1e-24
                        abs_xs = fis_xs + DATA[iso_i]['cap']*1e-24
                    except:
                        abs_xs = 0
                        fis_xs = 0
                    flux_depletion[i,j] += flux*abs_xs
                    flux_depletion[-1,j] += -flux*fis_xs
    return flux_depletion


def evaluate_flux(isotopes, conc, DATA, fixed=False):
    if fixed:
        flux = 0
        return flux

    # Assume power/mass
    power_per_mass = 35 # W/g or J/s/g
    fission_E = 200*1e6*EV2JOULE # J/fission
    # macro_per_density = 0 # cm^2/g
    macro = 0 # cm^-1
    for i, iso in enumerate(isotopes):
        try:
            fis_xs = DATA[iso]['fis']*1e-24
            # macro_per_density += fis_xs*MOL/nucname.anum(iso)
            macro += fis_xs*conc[i]/nucname.anum(iso)*MOL
        except:
            fis_xs = 0
    if macro > 0:
        flux = power_per_mass/(fission_E * macro)
    else: 
        flux = 0

    # print('Flux!')
    # print(flux)
    # print(type(flux))
    # print("{:.2e}".format(flux))
    return flux

def evaluate_kinf(isotopes, conc, DATA, reactor_type):
    numerator = 0
    denom = 0
    for i, iso in enumerate(isotopes):
        try:
            fis_xs = DATA[iso]['fis']*1e-24
            abs_xs = DATA[iso]['cap']*1e-24 + fis_xs
            macro_abs = abs_xs*conc[i]/nucname.anum(iso)*MOL
            denom += macro_abs
        except: 
            continue

        if reactor_type == 'thermal':
            if iso == '922350':
                fis_xs = DATA[iso]['fis']*1e-24
                abs_xs = DATA[iso]['cap']*1e-24 + fis_xs
                macro_fis = fis_xs*conc[i]/nucname.anum(iso)*MOL
                macro_abs = abs_xs*conc[i]/nucname.anum(iso)*MOL

                nubar = 2.4
                numerator += nubar*macro_fis
                denom += macro_abs

            elif iso == '942390':
                fis_xs = DATA[iso]['fis']*1e-24
                abs_xs = DATA[iso]['cap']*1e-24 + fis_xs
                macro_fis = fis_xs*conc[i]/nucname.anum(iso)*MOL
                macro_abs = abs_xs*conc[i]/nucname.anum(iso)*MOL

                nubar = 2.9
                numerator += nubar*macro_fis
                denom += macro_abs

        elif reactor_type == 'fast':
            if iso == '922350' or iso =='922380':
                fis_xs = DATA[iso]['fis']*1e-24
                abs_xs = DATA[iso]['cap']*1e-24 + fis_xs
                macro_fis = fis_xs*conc[i]/nucname.anum(iso)*MOL
                macro_abs = abs_xs*conc[i]/nucname.anum(iso)*MOL

                nubar = 2.6
                numerator += nubar*macro_fis
                denom += macro_abs

            elif iso == '942390':
                fis_xs = DATA[iso]['fis']*1e-24
                abs_xs = DATA[iso]['cap']*1e-24 + fis_xs
                macro_fis = fis_xs*conc[i]/nucname.anum(iso)*MOL
                macro_abs = abs_xs*conc[i]/nucname.anum(iso)*MOL

                nubar = 3.1
                numerator += nubar*macro_fis
                denom += macro_abs

    kinf = numerator/denom
    return kinf

def matrix_exp_solve(matrix, dt, conc):
    # print('exponent evaluation')
    # print(matrix)
    # print(np.exp(-matrix*dt))
    # print(conc)
    return np.matmul(expm(-matrix*dt),conc)

