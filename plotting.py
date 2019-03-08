import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

font = {'weight' : 'bold',
        'size'   : 13}
rc('font', **font)
rc('text', usetex=True)


def concentration_plot(isotopes, time, steps, conc_over_time, title):
        times = np.linspace(0, time, steps)
        styles = ['-','-.','.',':','v','^','<','>','s','p','h']
        n_styles = len(styles)
        
        for i, iso in enumerate(isotopes):
            plt.semilogy(times, conc_over_time[i, :],styles[i%n_styles])
            # plt.plot(times, conc_over_time[i, :])
        plt.legend(isotopes)
        plt.grid(True)
        plt.xlabel('Time (y)')
        plt.title(title)
        plt.show()

def concentration_subplot(isotopes, time, steps, conc_over_time, title):
        times = np.linspace(0, time, steps)
        # n_plots = U, Np, Pu, Am, Cm
        # 92 U, 93 Np, 94 Pu, 95 Am, 96 Cm
        iso_one_letter = ['U', 'N', 'P', 'A', 'C', 'O'] # O is for other
        n_plots = len(iso_one_letter)
        plot_axes = {}
        plot_legends = {}
        count = 0
        # fig, axarr = plt.subplots(int(n_plots/2), 2)
        for letter in iso_one_letter:
            count += 1
            plot_axes[letter] = plt.subplot(int(n_plots/2), 2, count)
            plot_legends[letter] = []
        
        for i, iso in enumerate(isotopes):
            try: 
                ax = plot_axes[iso[0]]
                ax.semilogy(times, conc_over_time[i, :])
                plot_legends[iso[0]].append(iso)
                # plt.plot(times, conc_over_time[i, :])
            except:
                ax = plot_axes['O']
                ax.semilogy(times, conc_over_time[i, :])
                plot_legends['O'].append(iso)
                continue

        for letter in iso_one_letter:
            ax = plot_axes[letter]
            ax.legend(plot_legends[letter])
            ax.grid(True)

        plt.xlabel('Time (y)')
        plt.show()
        
def flux_k_plot(time, steps, fluxes, ks):
    plt.subplot(2,1,1)
    times = np.linspace(0, time, steps)
    plt.plot(times[1:], fluxes)
    plt.grid(True)
    plt.xlabel('Time (y)')
    plt.ylabel('$\phi$ $(n/cm^2/s)$')

    plt.subplot(2,1,2)
    plt.plot(times[:], ks)
    plt.grid(True)
    plt.xlabel('Time (y)')
    plt.ylabel('$k_{\inf}$')
    plt.show()

def convergence_plot(isotopes, time_steps, final_concs):
    n_times = len(time_steps)
    # Plot only 235-U 239-Pu and FPs convergence
    rel_errors = np.zeros([3, n_times])
    rel_errors[:,0] = 1
    for i in range(n_times-1):
        for j, iso in enumerate(isotopes):
            if iso == 'U-235':
                prev_conc = final_concs[i][j]
                new_conc = final_concs[i+1][j]
                rel_errors[0,i+1] = np.abs((new_conc - prev_conc)/prev_conc)
            elif iso == 'Pu-239':
                prev_conc = final_concs[i][j]
                new_conc = final_concs[i+1][j]
                rel_errors[1,i+1] = np.abs((new_conc - prev_conc)/prev_conc)
            elif iso == 'FPs':
                prev_conc = final_concs[i][j]
                new_conc = final_concs[i+1][j]
                rel_errors[2,i+1] = np.abs((new_conc - prev_conc)/prev_conc)
    # plt.semilogy(time_steps, rel_errors[0,:],'.',markersize=12)
    # plt.semilogy(time_steps, rel_errors[1,:],'.',markersize=12)
    # plt.semilogy(time_steps, rel_errors[2,:],'.',markersize=12)
    plt.semilogy(time_steps, rel_errors[0,:])
    plt.semilogy(time_steps, rel_errors[1,:])
    plt.semilogy(time_steps, rel_errors[2,:])
    plt.legend(['U-235','Pu-239','Fission Products'])
    plt.title('Temporal Convergence of Concentrations')
    plt.ylabel('Absolute Relative Errror Between Final Concentrations')
    plt.xlabel('Total Number of Time Steps')
    plt.show()

        
