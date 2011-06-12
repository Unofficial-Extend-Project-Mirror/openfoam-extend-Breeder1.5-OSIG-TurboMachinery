#!/usr/bin/env python

import shutil, sys, os, glob
from subprocess import call
import getopt


# Definition of the different functions

def usage():
	print "Usage: you need three arguments to run this script. The first one should be the time step, the second one the log file, and the last one should be U or UMean depending on whether UMean is sampled or not."
	print """Example: runScript.py arg1 arg2 arg3. arg1 should be a timestep folder, arg2 should be the name of the log file and arg3 U or UMean."""
def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:v", ["help", "output="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    output = None
    verbose = False
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-o", "--output"):
            output = a
        else:
            assert False, "unhandled option"

if __name__ == "__main__":
	main()

# Those are the 2 arguments that the script requires. The first one should be a time step, while the second one should be a log file.
sampletime = sys.argv[1]
logFile = sys.argv[2]
Mean = sys.argv[3]

# set the current directory as working directory
origdir = os.getcwd()


# Create figures directory for holding all plots
if os.path.isdir('figures'):    	# does figures directory exist?
	shutil.rmtree('figures')	# yes, remove old directory
os.mkdir('figures')                     # create figures dir


# Construct the plots for the spiral windows
print "creating the plots for the draft tube windows"
cmd='python '+origdir+'/pythonScript/plot_draftTube_window.py ' + " " + sampletime + " " + Mean
retcode = call(cmd,shell=True)

#Extract the information from the log file
print "Creating the plots from the log file"
cmd='python '+origdir+'/pythonScript/plot_logInfo.py ' + " " + logFile
retcode = call(cmd,shell=True)

# Construct the plot for the probes
print "creating the plots for the different probes"
cmd='python '+origdir+'/pythonScript/plot_probes.py'
retcode = call(cmd,shell=True)

# Construct the plot for the turboPerformances
print "creating the plots for the turbo performances"
cmd='python '+origdir+'/pythonScript/plot_turbo.py'
retcode = call(cmd,shell=True)

    
# move.eps image files to figures directory
filelist=glob.glob('*.eps')
for file in filelist:
    shutil.move(file, origdir+"/figures/.")
filelist=glob.glob('*.png')
for file in filelist:
    shutil.move(file, origdir+"/figures/.")    
os.chdir("..")

