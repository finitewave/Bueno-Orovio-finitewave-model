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

# Convert to mV using model scaling
# V = np.array(model.history['u']) * 85.7 - 84.0

plt.plot(model.times, model.history['u'])
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Potential (mV)')
plt.title('0D Model Simulation')
plt.grid()
plt.show()

