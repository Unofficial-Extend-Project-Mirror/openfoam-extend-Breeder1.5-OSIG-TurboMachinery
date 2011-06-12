#!/usr/bin/env python
#Made by Olivier Petit, 2010. Any problem, contact Turbomachinery group.

description = """ Extract different informations from the log file of the simulation. The different informations, plotted against time are: max Courant number, mean Courant number, and Velocity magnitude """

import sys, os, math
from pylab import *
from numpy import loadtxt

# Definition of different files

infilename = sys.argv[1] # Argument that should be given in runScript.py, when this sub-script is called.
outfilename = './logInfo_tmp'
logfile = './logInfo'

def extractinfo():
#Extract the different informations described in the description, from the log file

    ifile = open(infilename, 'r')
    ofile = open(outfilename, 'w')
    
    lines = ifile.readlines()

    for line in lines:
        words = line.split()
        for word in words:
            index =words.index(word) 
            if word[0:] == 'Time':
                if words[index+1]=='=':
                    time = words[index+2]
                    ofile.write('%s ' % time)
                else:
                    print "the following lines will not be written (they are part of the definition of the different dictionnaries:\n" + line 
            if word[0:] == 'Courant':
                courantMean = words[index+3]
                courantMax = words[index+5]
                VMag = words[index+8]
                courant = courantMean + ' ' + courantMax + ' ' + VMag
                #print courant
                ofile.write('%s\n' % courant)
    ofile.close()
    
    ofile = open(outfilename, 'r')
    lfile = open(logfile, 'w')
    
    newlines = ofile.readlines()
    for x in range(1,len(newlines)-1):
        lines = newlines[x]
        lfile.write('%s' % lines)
        
    ofile.close()
    lfile.close()
    os.remove(outfilename)
    pass

def plotloginfo():

#Plot the different informations picked up against the time

    t, CoMean, CoMax, VMag = loadtxt("./logInfo", unpack = True)
    
    figure(1);
    plot(t, VMag)
    title("Velocity Magnitude vs time")
    xlabel("Time")
    ylabel("Velocity Magnitude")
    grid()
    savefig("VMAgVsT.eps", format = 'eps')
    
    figure(2);
    plot(t, CoMean)
    title("Mean courant number vs time")
    xlabel("Time")
    ylabel("Mean Courant Number")
    grid()
    savefig("CoMeanVsT.eps", format = 'eps')
    
    figure(3);
    plot(t, CoMax)
    title("Max courant number vs time")
    xlabel("Time")
    ylabel("Velocity Magnitude")
    grid()
    savefig("MaxCoVsT.eps", format = 'eps')
    
    pass


#Main code

extractinfo()
plotloginfo()
#show()
os.remove(logfile)
