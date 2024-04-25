"""
This script simulates the decay of 211Pb nuclei into different isotopes over time. 
It calculates the decay sequences for a large number of lead-211 (211Pb) nuclei and 
tracks the number of nuclei for each isotope in the decay chain: 
- 211Pb -> 211Bi -> 211Po -> 207Pb.
The decay constants for each decay step are provided (tau1, tau2, tau3). 
The simulation progresses in discrete time steps (dt) over a specified time interval (Nt). 

The script initializes arrays to store the decay times for each nucleus, 
as well as arrays to track the number of nuclei for each isotope at each time step. 

Decay State Check:
    - If the current time t[i] is less than the decay time t1[j] of nucleus j, it remains as 211Pb.
    - If the time falls between t1[j] and t1[j] + t2[j], nucleus j decays into 211Bi.
    - If the time falls between t1[j] + t2[j] and t1[j] + t2[j] + t3[j], nucleus j decays into 211Po.
    - If the time exceeds t1[j] + t2[j] + t3[j], nucleus j decays into stable 207Pb.

Update Counts:
    - The counts tn1, tn2, tn3, and tn4 are incremented based on the decay state of nucleus j at time t[i].

It then iterates through each time step, updating the number of nuclei for each isotope 
based on the decay times of individual nuclei.

The resulting decay curves for each isotope are plotted using Matplotlib.
"""

import numpy as np
import matplotlib.pyplot as plt
import random as rnd
from math import *

# Constants
N0 = 10000     # Initial number of 211Pb nuclei
tau1 = 3125    # Decay constant for 211Pb -> 211Bi
tau2 = 184.5   # Decay constant for 211Bi -> 211Po
tau3 = 0.746   # Decay constant for 211Po -> 207Pb

dt = 100       # Time interval for each step
Nt = 1000      # Number of time steps

# Initialize arrays to store decay times and number of nuclei for each isotope
t1 = np.zeros(N0)   # Decay time for 211Pb -> 211Bi
t2 = np.zeros(N0)   # Decay time for 211Bi -> 211Po
t3 = np.zeros(N0)   # Decay time for 211Po -> 207Pb

n1 = np.zeros(Nt+1)   # Number of 211Pb nuclei over time
n2 = np.zeros(Nt+1)   # Number of 211Bi nuclei over time
n3 = np.zeros(Nt+1)   # Number of 211Po nuclei over time
n4 = np.zeros(Nt+1)   # Number of 207Pb nuclei over time

# Decay times inverse transform sampling from an exp distribution
for i in range(0, N0):
    t1[i] = -tau1 * log(rnd.random())
    t2[i] = -tau2 * log(rnd.random())
    t3[i] = -tau3 * log(rnd.random())

# Decay simulation
t = np.linspace(0, dt * Nt, Nt+1)  # Time array
n1[0] = N0   # Initial number of 211Pb nuclei

# Loop over time steps
for i in range(1, Nt+1):
    tn1, tn2, tn3, tn4 = 0, 0, 0, 0
    
    # Loop over each nucleus
    # for each nucleus check the time against the decay time
    # to determine the state of decay for that nucleus at time t[i]
    for j in range(1, N0):
        # Count nuclei for each isotope based on decay times
        if t[i] < t1[j]:
            # Nucleus j remains as 211Pb
            tn1 += 1
        if t1[j] < t[i] < t1[j] + t2[j]:
            # Nucleus j decays into 211Bi
            tn2 += 1
        if t1[j] + t2[j] < t[i] < t1[j] + t2[j] + t3[j]:
            # Nucleus j decays into 211Po
            tn3 += 1
        if t[i] > t1[j] + t2[j] + t3[j]:
            # Nucleus j decays into stable 207Pb
            tn4 += 1
    
    # Update number of nuclei arrays at time t[i]
    n1[i] = tn1
    n2[i] = tn2
    n3[i] = tn3
    n4[i] = tn4

# Plot decay curves for each isotope
plt.plot(t, n1, 'k', label='211Pb')
plt.plot(t, n2, 'c', label='211Bi')
plt.plot(t, n3, 'g', label='211Po')
plt.plot(t, n4, 'r', label='207Pb')
plt.legend()
plt.title('211Pb Decay')
plt.xlabel('Time [s]')
plt.ylabel('Number of Nuclei')
plt.show()
