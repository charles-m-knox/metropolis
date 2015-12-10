import numpy as np
import math
import random
import matplotlib.pyplot as plt

#Start with some random microstate E_0

n_particles = 20
n_trials = 5000

temp = 2.0
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
        e_i = e_0               #Initialize
    e_j = e_i + (1 if randi() is 1 else -1)
    if e_j < 0:
        e_j = e_i + (1 if randi() is 1 else -1)
    else:
        #calculate dE
        dE = e_j - e_i
        if dE <= 0:
            e_i = e_j
        else:
            w = get_p(dE)
            if rand() <= w:
                e_i = e_j
            else:
                e_i = e_j*w

    e_list.append(e_i)

    #calculate <E>
    p_list = get_p_list(e_list)
    Z = sum(p_list)
    e_avg = n_particles * sum(mult(e_list, p_list)) / Z
    e_avg_list.append(e_avg)

plt.plot(trials, e_avg_list)
plt.savefig('algorithm_test.png')
