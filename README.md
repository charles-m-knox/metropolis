# Metropolis Algorithm for Energy, Temperature of Particles

## Usage

I'm currently running:

 python metropolis_algorithm.py -R 3 -n 20000 -p 20
 
and we will see the output when it's done.

## Introduction

Uses the metropolis algorithm (i.e. Monte Carlo) to simulate the progression of a system of n particles at various Temperature values.

## What it does

The program generates a random distribution of n particles, then does the Metropolis progression algorithm N number of times (i.e. 1000), then it plots the resulting probability vs energy plots for 10 different Temperatures.

It does this process 5 times and calculates the average energy, <E>, at each temperature, and then plots all 5 sets of <E> values.

## Output

![20 particles, distribution] (https://raw.githubusercontent.com/chuck-knox/metropolis/master/distribution_20_particles0.png)

![20 particles, <E>] (https://raw.githubusercontent.com/chuck-knox/metropolis/master/e_averages_20_particles.png)

![500 particles, <E>] (https://raw.githubusercontent.com/chuck-knox/metropolis/master/e_averages_500_particles.png)
