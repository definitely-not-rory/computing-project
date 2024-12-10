#IMPORTS
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.rcParams.update({'font.size': 14})

import matplotlib.font_manager as font_manager
plt.rcParams["axes.formatter.use_mathtext"]=True
plt.rcParams['font.family']='serif'
cmfont = font_manager.FontProperties(fname=mpl.get_data_path() + '/fonts/ttf/cmr10.ttf')
plt.rcParams['font.serif']=cmfont.get_name()
plt.rcParams['mathtext.fontset']='cm'
plt.rcParams['axes.unicode_minus']=False

pos_x_values=np.linspace(0,100,11)

vel_x_values=pos_x_values-5

pos_y_vals=5*np.ones(len(pos_x_values))

vel_y_vals=5.5*np.ones(len(vel_x_values))

fig,ax=plt.subplots()
'''
ax.set_ylim(4.5,6)
fig.set_figheight(3)
for i in range(-1,20):
    ax.axvline(i*5,c='k',linestyle='dashed',alpha=0.6)
ax.scatter(pos_x_values,pos_y_vals,c='r',s=200)
ax.scatter(vel_x_values,vel_y_vals,c='b',s=200)
ax.set_yticklabels([])
ax.text(-26,5,'Positions')
ax.text(-27,5.5,'Velocities')
ax.set_yticks([])
ax.set_xlabel('Time')

plt.show()'''
fig,ax=plt.subplots()
x_values=[100,50,25,20,10,5,2,1]
energy_values=[20.377916509189397,10.309808564070153,4.675415385187427,3.653676748687907,1.7494092466366389,0.8651408783329116,0.3097225504129889,0.12068003588746343]
rad_values=[2182173.75643804,4.6329776102827624,2.068899089766374,1.6522387093050206,0.8885825370285261,0.5126203155491615,0.2875939737229105,0.21205992585962008]
plt.scatter(x_values,energy_values)
ax.set_xscale('symlog')
ax.set_xlabel('log10(dt)')
ax.set_ylabel('Absolute Percentage Change in Energy')
plt.show()

fig,ax=plt.subplots()
x_values=[100,50,25,20,10,5,2,1]
energy_values=[20.377916509189397,10.309808564070153,4.675415385187427,3.653676748687907,1.7494092466366389,0.8651408783329116,0.3097225504129889,0.12068003588746343]
rad_values=[2182173.75643804,4.6329776102827624,2.068899089766374,1.6522387093050206,0.8885825370285261,0.5126203155491615,0.2875939737229105,0.21205992585962008]
plt.scatter(x_values,rad_values)
ax.set_xscale('symlog')
ax.set_yscale('symlog')
ax.set_xlabel('log10(dt)')
ax.set_ylabel('log10(Absolute Percentage Change from Circular Orbit)')
plt.show()

'''
PLOTS:
%energy as a func of time

intro text:

explain newtonian gravity

explain 2 problem solution

explain how >2 is unsolvable

explain numerical methods

method:

explain leapfrog

explain why higher order methods are less effective (store more variables, impractical for N particles)

explain energy-time assumptions (from project plan)

explain numerical issues w cons of energy

explain cons of momentum

explain CoM coord system

explain why its useful

Results:

plot of orbits

plot of energies

plot of earth deviation from circular


conclusions:

extension:

expand to an 'N' particle system

within a solar system context, incorporate further planets/bodies (moons, long period asteroids etc)

model more complex gravitational system - movement around black holes, galaxy dynamics, cluster collapse, virial theorem, formation of universe

use in non-astro fields, such as plasma, fluids and semiconductors

'''