#!/usr/bin/env python
#Made by Olivier Petit, 2010. Any problem, contact Turbomachinery group.

description = """ Compare the experimental and numerical radial/tangential velocity at the two windows located in the U9 spiral """

import matplotlib
import matplotlib.pyplot as plt
from pylab import *
from numpy import loadtxt
import math, os, getopt



# Definition of the different path where the data can be found

#sampletime = str(raw_input("Please enter the time step you want to sample: "))
sampletime = sys.argv[1]
Mean = sys.argv[2]
datadir = '../measurements/LDVdata/'
sampledir = './sets/'+sampletime+"/"
names = ['W_0', 'W_1', 'W_2']
legexp = {'W_0': "experimental",'W_1': "experimental", 'W_2': "experimental"}
legsim = {'W_0': "numerical",'W_1': "numerical", 'W_2': "numerical"}


# Definition of some constant

R = 0.050  # minimal radius of the test section
Q = 0.030  # discharge in m^3/s (30 L/s)
VThroat = Q / (math.pi*R**2) # mean velocity at the throat
Conv = pi / 180. #conversion degree to radian

# Definition of the functions

def mydraw():
    x, Ux, Uy, Uz = loadtxt(samplefile, unpack=True)
    length1, Um_Mean, Um_RMS, length2, Ut_Mean, Ut_RMS =loadtxt(datafile, unpack=True)
    
# Reduce to dimensionless dimensions, and transform from cartesian to cylindrical coordinates
   
    x = array(x)
    Ux = array(Ux)
    Uy = array(Uy)
    Uz = array(Uz)
    length1 = array(length1)
    length2 = array(length2)
    Um_Mean = array(Um_Mean)
    Um_RMS = array (Um_RMS)
    Ut_Mean = array(Ut_Mean)
    Um = (Uz * math.cos(beta * Conv) + Uy * math.sin(beta * Conv))
    Utang = -Ux * math.cos (angle * Conv)
    xprim = x / R
    Um_prim = Um / VThroat # dimensionless radial velocity
    Utang_prim = Utang / VThroat
        
# Plot the results
  
    leg_sim = ["Measured meridional velocity", "Measured tangential velocity"]
    leg_exp = ["Numerical meridional velocity", "Numerical tangential velocity"]

    plot(xprim,Um_prim, 'k',linewidth=3.0, label= leg_exp[0]) # Computed results from the sets file
    legend (loc='lower right')
 
    dashes = [5,2,10,5]
    l, = plot(xprim, Utang_prim, 'k--',linewidth=3.0, label= leg_exp[1]) # Computed results from the sets file
    l.set_dashes(dashes)
    legend (loc='lower right')
 
    plot (length1, Um_Mean, 'k_', linewidth = 4.0, label= leg_sim[0]) # Experimental results
    plt.errorbar(length1,Um_Mean, Um_RMS, fmt = 'k_', ecolor = 'k', elinewidth = '0.5')
    legend (loc='lower right')

    plot (length2, Ut_Mean, 'k:', linewidth = 4.0, label= leg_sim[1]) # Experimental results
    plt.errorbar(length1,Ut_Mean, Ut_RMS, fmt = 'k:', ecolor = 'k', elinewidth = '0.5')
    legend (loc='lower right')

# Definition of the plot parameters
   
    grid()
    ylabel(r"Dimensionless velocity [-]",fontsize=14)
    xlabel(r"Dimensionless survey axis [-]",fontsize=14)
    ax = list(axis())
    ax[2] = -2
    ax[0] = 0
    axis(ax)
    savefig("window_"+n+".png",format='png')
    #show()
    clf()
    pass


# Creation of the plot using the function mydraw()

for n in names:
    if n == 'W_0':
        beta = 25
        angle = 0
    elif n == 'W_1':
        beta = 8.5
        angle = 0
    elif n=='W_2' :
        beta = 8.5
        angle = 180
    betastr = str(beta)
    print "plot for window " + n + ",the angle of the window is: " + betastr
    
# Definition of the path of the file to load
   
    xy_files = "_" + Mean + ".xy"
    samplefile = sampledir + n + xy_files
    datafile = datadir + n + "_tot"
    mydraw() # Call the function draw defined further up
    
