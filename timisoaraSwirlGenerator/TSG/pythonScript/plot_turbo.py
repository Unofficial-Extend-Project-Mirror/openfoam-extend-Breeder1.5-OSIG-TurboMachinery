#!/usr/bin/env python
#Made by Olivier Petit, 2010. Any problem, contact Turbomachinery group.

description = """ Plot the different information available from the turboPerformance functionObject.\n
Those information are: the head of the turbine, the torque, the efficiency and the forces plotted against the time. """

import os
import math
from pylab import *
from numpy import loadtxt

# Function definitions

def addToPlots( timeName ):
	fileName =  "./turboPerformance/" + timeName + '/turboPerformance.dat'
	time, head, TOmega, eff, Fx, Fy, Fz = loadtxt( fileName ,skiprows=1,usecols=(0,1,2,3,4,5,6),unpack=True )
	
	figure(1);
	plot( time, head )  # fig1 active
	ylabel(' Head (m) ');   xlabel(' Time (s) '); title(' Head vs. Time ')
	grid()
	savefig("HeadVsTime.eps",format='eps')
	
	figure(2); 
	plot( time, TOmega ) # fig2 active
	ylabel(' TOmega (W) '); xlabel(' Time (s) '); title(' Axial Power vs. Time ')
	grid()
	savefig("TOmegaVsTime.eps",format='eps')
	
	figure(3); 
	plot( time, eff )   # fig3 active	
	ylabel(' Eff (%) '); xlabel(' Time (s) '); title(' Efficiency vs. Time ')
	grid()
	savefig("EffVsTime.eps",format='eps')
	
	figure(4); 
	plot( time, Fx, time, Fy ) # fig4 active     
	ylabel(' Forces (N) '); xlabel(' Time (s) '); title(' Axial Forces vs. Time ')
	grid()
	savefig("ForcesVsTime.eps",format='eps')
	

# - - - Main program - - - - -

for dirStr in os.listdir("./turboPerformance/"):
	addToPlots( dirStr )
show()
