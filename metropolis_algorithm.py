#
#   Metropolis Algorithm Implementation
#   1. Generate n particles with energy 1<E<E_Max
#   2. For N trials, on each trial, change the energy of 1 randomly chosen 
#       particle by +1 or -1, also chosen randomly. Then, calculate 
#       the TOTAL energy of the system and compare it to the energy BEFORE 
#       the change was made. If the change is less than 0, proceed. If the 
#       energy change dE is greater than 0, proceed only if the value of some 
#       random number 0<r<1 is less than/equal to the value of exp(-beta*dE). 
#       Otherwise, scrap the change and continue to the next trial.
#   3. Put the results in plottable, sorted form and plot them.
#
# Code Standards: https://www.python.org/dev/peps/pep-0008
#

#import matplotlib.mlab as mlab
#from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import argparse
import collections
import numpy as np
import math
import random
from time import strftime

def log(message):
    print message

def gen_mcs_dist(temp, num_particles, num_trials):
    n_particles = num_particles
    n_trials = num_trials

    temp = float(temp)
    beta = 1.0/temp

    e_0 = 0.0
    e_i = 0.0
    e_j = 0.0

    e_list = []
    e_avg_list = []

    def randi():
        return random.randint(0,1)

    def rand():
        return random.random()

    def get_p(dE):
        return math.exp(-1.0*beta*dE)

    def get_p_list(e_vals):
        p_list = []
        for e_val in e_vals:
            p_list.append(get_p(e_val))
        return p_list

    def mult(a, b):
        list_out = []
        if len(a) != len(b):
            print "got lists of different length"
        else:
            for i in range(0, len(a)):
                list_out.append(a[i] * b[i])
        return list_out

    trials = range(0, n_trials)

    for x in trials:
        if x is 0:
            e_i = e_0                               #Initialize
        e_j = e_i + (1 if randi() is 1 else -1)     #Chagne by +/- 1
        if e_j < 0:
            e_j = e_i + (1 if randi() is 1 else -1) #Retry for negative 
                                                    #energies
        else:
            dE = e_j - e_i                          #calculate dE
            if dE <= 0:                             
                e_i = e_j                           #Accept if -1
            else:
                w = get_p(dE)                       #Accept +1
                if rand() <= w:                     # only if w
                    e_i = e_j                       # is <= random
                else:                               # between 0,1
                    e_i = e_j*w                     #Otherwise accept +1*w (<1)

        e_list.append(e_i)
        p_list = get_p_list(e_list)                 #calculate <E> for each
        Z = sum(p_list)                             # trial
        e_avg = n_particles * sum(mult(e_list, p_list)) / Z
        e_avg_list.append(e_avg)
    return e_avg_list

def gen_complete_distributions(num_particles, 
                                plots_temperatures, 
                                num_trials,
                                run_number):
    e_single_average_list = []
    e_lists = []
    for z in range(0, len(plots_temperatures)):
        log('Starting: Monte Carlo run number ' + str(
                run_number) + ', temperature is ' + str(plots_temperatures[z]))
        e_list = gen_mcs_dist(plots_temperatures[z], 
                num_particles, 
                num_trials)
        e_avg_val = e_list[len(e_list) - 1]
        e_single_average_list.append(e_avg_val) #Try last E_avg
        e_lists.append(e_list)
    return e_lists, e_single_average_list

def main(args_dict):
    #--------------------------------------------------------------------------
    num_particles = int(args_dict['num_particles'])
    plots_temperatures = [1,5,10,50,100,500,1000,5000,10000,100000]
    num_trials = int(args_dict['num_trials'])
    trials = range(0, num_trials)
    runs = int(args_dict['num_runs'])
    log(strftime("%Y-%m-%d %H:%M:%S"))
    log('Doing ' + str(runs) + ' total simulations with ' + str(
            len(plots_temperatures)) + ' temperature variations per run')
    #--------------------------------------------------------------------------
    average_energies_lists = []                    #Later we will plot <E>/T
    e_lists_lists = []                             #ex:e_lists_lists[run][temp]
    #----------------------------------------------#Do multiple simulations
    for xx in range(0, runs):                      # and store the results
        e_lists, e_single_average_list = gen_complete_distributions(
                                            num_particles, 
                                            plots_temperatures, 
                                            num_trials,
                                            xx)
        average_energies_lists.append(e_single_average_list)
        e_lists_lists.append(e_lists)
    #-----------------------------------------------------------Figure1 config
    fig = plt.figure(figsize=(8.5,11))
    plot_index = 0
    for z in range(0, len(plots_temperatures)):
        plot_index += 1
        #--------------------------------------------------------Subplot config
        plt.subplot(5, 2, plot_index)
        plt.title(r"$T=" + str(plots_temperatures[z]) + r"$", fontsize=12)
        if plot_index == 9:
            plt.ylabel(r"$\mathrm{\langle E \rangle}$", 
                        fontsize=12)
        if plot_index == 9:
            plt.xlabel(r"$\mathrm{Monte\,Carlo\,Steps}$", fontsize=12)
        plt.tick_params(axis='y', which='major', labelsize=5)
        plt.tick_params(axis='y', which='minor', labelsize=5)
        plt.tick_params(axis='x', which='major', labelsize=5)
        plt.tick_params(axis='x', which='minor', labelsize=5)
        plt.grid(axis="both", alpha=0.10, linestyle="-")
        for xx in range(0, runs):                               #Plot All Runs 
            plt.scatter(trials,                                 #at Once!
                        e_lists_lists[xx][z], 
                        s=0.25, 
                        alpha=0.10)
        #-------------------------------------------------------Subplot Config
    fig.suptitle(r"$\mathrm{Trials}=" + str(num_trials) + r",\,n=" +
                    str(num_particles) + r"\mathrm{\,Particles,\,}" + 
                    str(runs) + r"\mathrm{\, Runs}$", fontsize=18)
    plt.savefig("mcs_" + str(num_particles) + "_particles_" + str(num_trials) + "_trials.png")
    #plt.savefig("mcs_" + str(num_particles) + "_particles_run" + 
    #                str(run_number) + ".svg")
    #-----------------------------------------------------------Figure1 End
    #
    #-----------------------------------------------------------Figure2 Begin
    fig = plt.figure(figsize=(8.5,11))
    for average_energy_list in average_energies_lists:
        plt.loglog(plots_temperatures, average_energy_list)
    analytical_energy_averages = []
    for temp in plots_temperatures:
        exp_val = math.exp(-1.0/temp)
        analytical_energy_value = exp_val / (1.0 - exp_val)
        analytical_energy_averages.append(analytical_energy_value)
    plt.loglog(plots_temperatures, analytical_energy_averages)
    #-----------------------------------------------------------Figure2 config
    fig.suptitle(r"$\langle E\left(T\right)\rangle \mathrm{\,values:\,}n=" + 
                    str(num_particles) + r"\mathrm{,\,Trials}=" + 
                    str(num_trials) + r"$", fontsize=18)
    plt.ylabel(r"$\langle E\left( T\right) \rangle \mathrm{\, values}" + 
                    r"$", fontsize=16)
    plt.xlabel(r"$\mathrm{Temperature\,(K)}$", fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.grid(axis="both", which='major', alpha=0.25, linestyle="-")
    plt.grid(axis="both", which='minor', alpha=0.10, linestyle="-")
    plt.savefig("e_averages_" + str(num_particles) + "_particles.svg")
    plt.savefig("e_averages_" + str(num_particles) + "_particles.png")
    log("Finished")
    log(strftime("%Y-%m-%d %H:%M:%S"))
    #-----------------------------------------------------------Figure2 config

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Does the Metropolis algorithm and exports plots.')
    parser.add_argument("-n", "--num-trials", default=5000, 
        help="Number of times to randomly alter energy")
    parser.add_argument("-p", "--num-particles", default=20, 
        help="Number of particles in the macrostate (system)")
    parser.add_argument("-R", "--num-runs", default=5, 
        help="Number of Monte Carlo simulations to do (recommend 3-5+")
    args_dict = vars(parser.parse_args())
    main(args_dict)
