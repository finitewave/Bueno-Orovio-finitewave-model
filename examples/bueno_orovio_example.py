"""
Example script to run a 0D model simulation and plot the results.

This script sets up a simple stimulation protocol, runs the simulation,
and plots the membrane potential over time.
"""

import numpy as np
import matplotlib.pyplot as plt

from implementation.bueno_orovio_0d import BuenoOrovio0D, Stimulation


stimulations = [Stimulation(t_start=0.1, duration=0.2, amplitude=5.0)]
t_max = 300.0

model = BuenoOrovio0D(dt=0.01, stimulations=stimulations)
model.run(t_max=t_max)

time = np.arange(0, t_max, model.dt)
plt.plot(time, model.history['u'])
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (u)')
plt.title('0D Model Simulation')
plt.grid()
plt.show()

