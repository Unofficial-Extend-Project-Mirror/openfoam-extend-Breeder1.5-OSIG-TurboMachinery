#!/usr/bin/env python
#Made by Olivier Petit, 2010. Any problem, contact Turbomachinery group.

description = """ Plot the pressure information that comes from the probes, defined in system/controlDict. In this case, there are 30 different probes."""

import os, sys
import math
from pylab import *
from numpy import loadtxt

# Function definitions 

def addToPlots( timeName ):
	fileName =  "./probes1/" + timeName + '/p'
	probename = [ "p1", "p2", "p3", "p4", "p5", "p6", "p7"]
	probenumber = { "p1": '1', "p2":'2', "p3":'3', "p4":'4', "p5":'5', "p6":'6', "p7":'7'}
	
	for i in probename:
		x = probenumber [i]
		x = float(x)
		i=[]
		time=[]
		abc =loadtxt(fileName, skiprows=4)
		for z in abc:
			time.append(z[0])
			i.append(z[x])
		figure(x);
		leg = str(x)
		legend = "probe" + leg + "VsTime.eps"
		plot(time,i)
		ylabel(' p/rho ');   xlabel(' Time (s) '); title(' probe ' + leg + ' vs. Time ')
		grid()
		savefig(legend,format='eps')
	
	
	
	
for dirStr in os.listdir("./probes1/"):
	addToPlots( dirStr )
	
#Activate this line if show is needed
show()
