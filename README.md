# Metropolis Algorithm for Energy, Temperature of Particles

## Usage

I'm currently running:

    python metropolis_algorithm.py -R 8 -n 2500 -p 20


## Introduction

Uses the metropolis algorithm (i.e. Monte Carlo) to simulate the progression of a system of n particles at various Temperature values.

Requires matplotlib, numpy, python2

## What it does

The program generates a random distribution of n particles, then does the Metropolis progression algorithm n number of times (i.e. 1000), then it plots the resulting probability vs energy plots for 10 different Temperatures.

It does this process R times and calculates the average energy, at each temperature, and then plots all R sets of  values along with an analytical calculation.

## Output

![200 particles] (https://raw.githubusercontent.com/chuck-knox/metropolis/master/output/png/mcs_200_particles_2500_trials.png)

![200 particles, <E>] (https://raw.githubusercontent.com/chuck-knox/metropolis/master/output/png/e_averages_200_particles.png)


## Ideas

To make the process more streamlined, each run could be sent on its own thread. I'd do it if I had had more time before having to turn in this project.
