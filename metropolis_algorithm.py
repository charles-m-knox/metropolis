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

def log(message):
    print message

def get_complete_plot(temp, num_particles, num_trials, energy_max):
    beta = (1/temp) #Thermodynamic Beta = 1/kT
    
    energy_progression = [] #For each trial, we will record the total energy

    def get_highest_energy(particle_collection):
        highest_energy_val = 0
        try:
            for particle in particle_collection:
                if particle['E'] > highest_energy_val:
                    highest_energy_val = particle['E']
            return int(highest_energy_val)
        except Exception as e:
            log("Failed to get highest energy: " + str(e))
            return int(energy_max)

    def get_random_particle_index(num_particles):
        try:
            return random.randint(0,num_particles)
        except Exception as e:
            log('Failed to generate random particle index given ' 
                + str(num_particles) + ' as input:' + str(e))
            return 0

    def alter_energy_of_particle(particle, energy_change):
        try:
            particle['E'] += int(energy_change)
            particle['P'] = math.exp(-(1/temp) * float(particle['E']))
        except Exception as e:
            log('Failed to alter energy of particle:' + str(e))

    def alter_random_particle(particle_collection, energy_change):
        try:
            random_particle_index = get_random_particle_index(
                                        len(particle_collection)-1) 
            alter_energy_of_particle(
                particle_collection[random_particle_index],
                energy_change)
        except Exception as e:
            log('Failed to alter random particle: ' + str(e))

    def generate_particle(energy_max):
        return_dict = {'P': 0, 'E': int(1)}
        try:
            return_dict['E'] = int(random.randint(1, energy_max))
            return_dict['P'] = math.exp(
                                    -(1/temp)
                                    * float(return_dict['E']))
        except Exception as e:
            log('Failed to generate particle: ' + str(e))
        return return_dict

    def generate_particle_list(num_particles, energy_max):
        return_list = []
        try:
            for x in range(0,num_particles):
                particle = generate_particle(energy_max)
                return_list.append(particle)
        except Exception as e:
            log('Failed to generate particle list:' + str(e))
        return return_list

    def get_total_energy(particle_collection):
        total_energy = 0
        try:
            for particle in particle_collection:
                total_energy += particle['E']
        except Exception as e:
            log('Failed to get total energy of macrostate:' + str(e))
        return total_energy

    def get_random_energy_change_value():
        try:
            if random.randint(0,1) == 1:
                return 1
            else:
                return -1
        except Exception as e:
            log("Failed to get random energy change value:" + str(e))
            return 0

    #Step 1: Initialize an dictionary list of N length, with each item
    # containing { P_i, E_i } - these are particles.
    particle_collection = generate_particle_list(num_particles, energy_max)
    
    #Step 2, N times: Change energy of random particle, compute total 
    # energy of system, then accept or do not accept.
    for i in range(0, num_trials):
        particle_collection_unaltered = particle_collection
        unaltered_macrostate_energy = get_total_energy(particle_collection)
        microstate_energy_change = get_random_energy_change_value()
        alter_random_particle(particle_collection, microstate_energy_change)
        altered_macrostate_energy = get_total_energy(particle_collection)
        dE = altered_macrostate_energy - unaltered_macrostate_energy
        # If the energy change is DE<0, accept the change
        #if dE < 0:
        #    log("Accepting change")
        #If the energy change is DE>0, accept the change only if 
        # r<exp(-b*DE);  where r is a randomly generated number 0<r<1:
        if dE > 0:
            random_num_r = random.random()
            energy_condition_value = math.exp(-(beta) * 
                dE)
            if random_num_r > energy_condition_value:
                particle_collection = particle_collection_unaltered
        elif dE == 0:
            log("Got zero for energy change. Doing nothing.")
        total_macrostate_energy = get_total_energy(particle_collection)
        energy_progression.append(total_macrostate_energy)

    #Count everything up.
    count_collection = []
    highest_energy = get_highest_energy(particle_collection)
    for x in range(0, highest_energy):
        count_collection.append(0)
    for particle in particle_collection:
        try:
            count_collection[int(particle['E']) - 1] += 1 
        except Exception as e:
            log("Failed to count energy of " + str(particle['E']) 
                + ": " + str(e))
    count_collection.sort()
    count_collection.reverse()
    return range(0, highest_energy), count_collection, 

def main(args_dict):
    fig = plt.figure(figsize=(8.5,11))
    plot_index = 0
    plots_particles = [20,20,100,100,1000,1000,10000,10000,100000,100000]
    plots_temperatures = [100,5500,100,5500,100,5500,100,5500,100,5500]
    maximum_energy_all_plots = int(args_dict['energy_max'])
    trials = int(args_dict['num_trials'])
    fig.suptitle(r"$\mathrm{Trials}="+str(trials)+r"$", fontsize=18)
    for z in range(0,10):
        plot_index += 1
        plt.subplot(5, 2, plot_index)
        plt.title(r"$n=" + str(plots_particles[z]) +
                    r",\,T=" + str(plots_temperatures[z])+r"$", fontsize=12)
        #Decorative purposes; only show y-axis label on left, x label on bottom
        if plot_index % 2 != 0:
            plt.ylabel(r"$\mathrm{Count}$", fontsize=12)
        if plot_index >= 9:
            plt.xlabel(r"$\mathrm{Energy}$", fontsize=12)
        x_vals, y_vals = get_complete_plot(
                plots_temperatures[z], 
                plots_particles[z], 
                trials, 
                maximum_energy_all_plots)
        plt.xticks([0,max(x_vals)])
        plt.yticks([0,max(y_vals)])
        plt.tick_params(axis='both', which='major', labelsize=6)
        plt.grid(axis="y", alpha=0.25, linestyle="-")
        plt.scatter(x_vals, y_vals, s=2, c="r")
    plt.savefig("figure_out.png")
    plt.savefig("figure_out.svg")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Does the Metropolis algorithm and exports plots.')
    parser.add_argument("-n", "--num-trials", default=100, 
        help="Number of times to randomly alter energy")
    parser.add_argument("-E", "--energy-max", default=100, 
        help="Maximum energy value, with 1 being the lowest energy")
    #parser.add_argument("-p", "--num-particles", default=1, 
    #    help="Number of particles in the macrostate (system)")
    #parser.add_argument("-T", "--temperature", default=300, 
    #    help="Desired equilibrium temperature")
    args_dict = vars(parser.parse_args())
    main(args_dict)
