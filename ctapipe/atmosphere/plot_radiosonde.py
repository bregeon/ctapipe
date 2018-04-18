#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot radiosonde data

@author: Johan Bregeon
"""

# @TODO
# * plot all samples for P and T on one graph
# * plot Dewpoint, Wind Dir and Wind Speed

from ctapipe.utils import get_dataset
from matplotlib import pyplot as plt

from ctapipe.atmosphere.radiosonde_reader import pRadioSondeReader

# Read input data
reader=pRadioSondeReader()
reader.readFile(get_dataset('radiosonde_68110_2009_2014.fsl'))
prof=reader.getSounding('2009-08-23')
prof.dump()

fig = plt.figure(figsize=(6, 8))
ax  = fig.add_subplot(111)

pressure, = ax.plot(prof.Pressure, prof.Height, color='blue')
temp,     = ax.plot(prof.Temperature, prof.Height, color='green')
legend_handles = [pressure, temp]
legend_labels = ['Pressure', 'Temperature']


plt.title('Atmospheric (P, T) Profiles')
ax.set_xlabel("Pressure, Temperature")
ax.set_xbound(0, 1100)
ax.set_ylabel("Altitude")
ax.set_ybound(0, 40000)

ax.legend(legend_handles, legend_labels,
          bbox_to_anchor=(1.02, 1.), loc=1,
          borderaxespad=0., fontsize=10)
plt.show()
