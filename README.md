# Metropolis Algorithm for Energy, Temperature of Particles

## Usage

I'm currently running:

    python metropolis_algorithm.py -R 8 -n 2500 -p 20


## Introduction

Uses the metropolis algorithm (i.e. Monte Carlo) to simulate the progression of a system of n particles at various Temperature values.

## What it does

The program generates a random distribution of n particles, then does the Metropolis progression algorithm N number of times (i.e. 1000), then it plots the resulting probability vs energy plots for 10 different Temperatures.

It does this process 5 times and calculates the average energy, <E>, at each temperature, and then plots all 5 sets of <E> values.

## Output

![20 particles, progression] (https://raw.githubusercontent.com/chuck-knox/metropolis/master/mcs_20_particles_2500_trials.png)

![20 particles, <E>] (https://raw.githubusercontent.com/chuck-knox/metropolis/master/e_averages_20_particles.png)
