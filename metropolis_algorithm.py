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

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import argparse
import collections
import numpy as np
import math
import random
from scipy.optimize import curve_fit

def log(message):
    print message

def gen_mcs_dist(temp, num_particles, num_trials, energy_max):
    beta = float(1/float(temp)) #Thermodynamic Beta = 1/kT
    
    def get_random_energy_change_value():
        try:
            if random.randint(0,1) == 1:
                return 1
            else:
                return -1
        except Exception as e:
            log("Failed to get random energy change value:" + str(e))
            return 0

    particle_collection = []
    for i in range(0, num_particles):
        particle_collection.append(random.randint(1, energy_max))

    #Step 2, N times: Change energy of random particle, compute total 
    # energy of system, then accept or do not accept.
    steps = []
    energy_steps = []
    num_rejects = 0
    for i in range(0, num_trials):
        particle_collection_unaltered = []
        for particle in particle_collection:
            particle_collection_unaltered.append(particle)
        #Calculate total energy of system
        initial_macrostate_energy = sum(particle_collection)
        microstate_energy_change = get_random_energy_change_value()
        #Alter the energy of a random particle
        random_particle_index = random.randint(0, len(particle_collection) - 1)
        if particle_collection[random_particle_index] > 1:
            particle_collection[random_particle_index] += microstate_energy_change
        altered_macrostate_energy = sum(particle_collection)
        dE = altered_macrostate_energy - initial_macrostate_energy
        if dE > 0:
            random_num_r = random.random()
            energy_condition_value = math.exp(-(beta) * dE)
            if random_num_r > energy_condition_value:
                num_rejects += 1
                particle_collection = particle_collection_unaltered
        total_macrostate_energy = sum(particle_collection)
        steps.append(i)
        energy_steps.append(total_macrostate_energy)

    probabilities_list = []
    partition_fct = float(0)
    #First get the partition function by adding all particles' exponential fcts
    for particle in particle_collection:
        partition_fct += math.exp(-beta * particle)
    for particle in particle_collection:
        probabilities_list.append(1/partition_fct * math.exp(-beta * particle))
    return particle_collection, probabilities_list, steps, energy_steps

def gen_complete_distributions(num_particles, 
                plots_temperatures, 
                maximum_energy_all_plots, 
                num_trials,
                run_number):
    complete_data_set = []
    monte_carlo_steps_array = []
    fig = plt.figure(figsize=(8.5,11))
    plot_index = 0
    for z in range(0,10):
        plot_index += 1
        plt.subplot(5, 2, plot_index)
        plt.title(r"$T=" + str(plots_temperatures[z])+r"$", fontsize=12)
        #Decorative purposes; only show y-axis label on left, x label on bottom
        if plot_index % 2 != 0:
            plt.ylabel(r"$\mathrm{Probability}$", fontsize=12)
        if plot_index >= 9:
            plt.xlabel(r"$\mathrm{Energy}$", fontsize=12)
        x_vals, y_vals, x_steps, y_steps = gen_mcs_dist(
                plots_temperatures[z], 
                num_particles, 
                num_trials, 
                maximum_energy_all_plots)
        monte_carlo_steps_array.append([x_steps,y_steps])
        #plt.xticks([0,max(x_vals)])
        #plt.yticks([0,max(0.11)])
        plt.tick_params(axis='y', which='major', labelsize=5)
        plt.tick_params(axis='y', which='minor', labelsize=5)
        plt.tick_params(axis='x', which='major', labelsize=3)
        plt.tick_params(axis='x', which='minor', labelsize=3)
        plt.grid(axis="both", alpha=0.10, linestyle="-")
        plt.scatter(x_vals, y_vals, s=5, c="r")
        complete_data_set.append({
            'E': x_vals, 
            'P': y_vals, 
            'T': plots_temperatures[z],
            'N': num_particles,
            'E_avg': 0})
    fig.suptitle(r"$\mathrm{Trials}=" + str(num_trials) + r",\,n=" +
                    str(num_particles) + r"\mathrm{\,Particles\,-\,Run\," + str(run_number) + r"}$", fontsize=18)
    plt.savefig("distribution_" + str(num_particles) + "_particles" + 
                    str(run_number) + ".png")
    plt.savefig("distribution_" + str(num_particles) + "_particles" + 
                    str(run_number) + ".svg")
    #Render an entirely different set of plots for Monte Carlo Steps
    fig = plt.figure(figsize=(8.5,11))
    plot_index = 0
    for z in range(0,10):
        plot_index += 1
        plt.subplot(5, 2, plot_index)
        plt.title(r"$T=" + str(plots_temperatures[z])+r"$", fontsize=12)
        #Decorative purposes; only show y-axis label on left, x label on bottom
        if plot_index % 2 != 0:
            plt.ylabel(r"$\mathrm{Total\,Energy}$", fontsize=12)
        if plot_index >= 9:
            plt.xlabel(r"$\mathrm{Monte\,Carlo\,Steps}$", fontsize=12)
        plt.tick_params(axis='y', which='major', labelsize=5)
        plt.tick_params(axis='y', which='minor', labelsize=5)
        plt.tick_params(axis='x', which='major', labelsize=5)
        plt.tick_params(axis='x', which='minor', labelsize=5)
        plt.grid(axis="both", alpha=0.10, linestyle="-")
        plt.scatter(monte_carlo_steps_array[z][0], 
                    monte_carlo_steps_array[z][1], 
                    s=1, 
                    c="r")
    fig.suptitle(r"$\mathrm{Trials}=" + str(num_trials) + r", n=" +
                    str(num_particles) + r"\mathrm{\,Particles}$", fontsize=18)
    plt.savefig("mcs_" + str(num_particles) + "_particles" + 
                    str(run_number) + ".png")
    #plt.savefig("mcs_" + str(num_particles) + "_particles" + 
    #                str(run_number) + ".svg")


    #Now compute the average energy <E> for each data set
    E_averages = []
    for data_set in complete_data_set:
        #<E>=Sum(E_i * P_i)
        avg_E = float(0)
        for x in range(0, num_particles):
            try:
                avg_E += data_set['E'][x] * math.exp(-(1/float(data_set['T']) *
                                                    data_set['E'][x]))
            except Exception as e:
                log('Failed to get avg E, x='+str(x)+',len(data_set[E])=' + 
                        str(len(data_set['E']))+ ': ' + str(e))
        data_set['E_avg'] = avg_E
        E_averages.append(avg_E)
    return complete_data_set, E_averages

def main(args_dict):
    num_particles = int(args_dict['num_particles'])
    plots_temperatures = [1,5,10,50,100,500,1000,5000,10000,100000]
    maximum_energy_all_plots = int(args_dict['energy_max'])
    num_trials = int(args_dict['num_trials'])
    runs = 1
    average_energies = []
    for xx in range(0, runs):
        complete_data_set, E_averages = gen_complete_distributions(
                                            num_particles, 
                                            plots_temperatures, 
                                            maximum_energy_all_plots, 
                                            num_trials,
                                            xx)
        average_energies.append(E_averages)
    fig = plt.figure(figsize=(8.5,11))
    for average_energy_list in average_energies:
        plt.loglog(plots_temperatures, average_energy_list)
    #calculate analytical <E> values
    analytical_energy_averages = []
    for temp in plots_temperatures:
        exp_val = math.exp(-(float(1)/float(temp)))
        analytical_energy_value = exp_val / (1 - exp_val)
        analytical_energy_averages.append(analytical_energy_value)
    plt.loglog(plots_temperatures, analytical_energy_averages)
    fig.suptitle(r"$\langle E\left(T\right)\rangle \mathrm{\,values:\,}n=" + 
                    str(num_particles) + r"\mathrm{,\,Trials}=" + 
                    str(num_trials) + "$", fontsize=18)
    plt.ylabel(r"$\langle E\left( T\right) \rangle \mathrm{\, values}" + 
                    r"$", fontsize=16)
    plt.xlabel(r"$\mathrm{Temperature\,(K)}$", fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.grid(axis="both", which='major', alpha=0.25, linestyle="-")
    plt.grid(axis="both", which='minor', alpha=0.10, linestyle="-")
    plt.savefig("e_averages_" + str(num_particles) + "_particles.svg")
    plt.savefig("e_averages_" + str(num_particles) + "_particles.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Does the Metropolis algorithm and exports plots.')
    parser.add_argument("-n", "--num-trials", default=10000, 
        help="Number of times to randomly alter energy")
    parser.add_argument("-E", "--energy-max", default=1000, 
        help="Maximum energy value, with 1 being the lowest energy")
    parser.add_argument("-p", "--num-particles", default=20, 
        help="Number of particles in the macrostate (system)")
    #parser.add_argument("-T", "--temperature", default=300, 
    #    help="Desired equilibrium temperature")
    args_dict = vars(parser.parse_args())
    main(args_dict)
