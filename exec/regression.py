#!/usr/bin/python
import os
import sys


# This script verifies that MPI and Serial builds produce identical results
# for a specified set of rundecks.
# Usage:
#    From the decks subdirectory issue the command:
#      ../exec/regression.py  <runsource1> [<rundeck2> ...]
#
#    Note that you must issue the command from a batch process so that MPI can be
#    used.
#
# Options:
#    * If the environment variable VERBOSE is set to True, then the script
#      will display all commands as they execute.
#    * If the environment variable DEBUG is set to True, then the script will
#      display all commands, but not actually execute them.
#
# Issues:
#   1) The script currently leaves a fair bit of detritus lying around in the decks subdirectory.
#      cleanup command should be added in the future.
#   2) Currently the details of failures are sent to /dev/null.   An extra log file (or set of files) should
#      eventually be managed to contain the details and leave STDOUT to handle the big picture items.


# Issue a shell command and raise an exception if the result is not 0.
# Environment variables:
#     VERBOSE=True  displays commands in the log
#     DEBUG=True    displays commands in the log, but does not actually
#                   submit to system.
def systemCommand(commandString):
    if os.environ.has_key('DEBUG'):
        debug = os.environ['DEBUG']
    else:
        debug = False
        
        if debug:
            os.system("echo "+commandString)
        else:
            if os.environ.has_key('VERBOSE') and os.environ['VERBOSE']:
                print "   Command: " + commandString

            status = os.system(commandString + " &> /dev/null")
            if (status != 0):
                raise Exception('unix', commandString)

# This procedure returns a suitable rundeck name for
# a given run source and mode (mpi,serial, openmp, hybrid).
# Must be consistent across other procedures, but otherwise
# is arbitrary.
# For now it assumed to be specified by the user
def rundeckName(runSource, mode):
    return runSource

def checkpointFileName(runCase, duration):
    if runCase['mode'] == 'serial':
        return runCase['rundeck'] + "_" + runCase['mode'] + "_" + duration
    else:
        return runCase['rundeck'] + "_" + runCase['mode'] + "_" + duration + "_NPES=" + runCase['NPES']

# Return a dict that specifies the configuration under test
def newConfiguration(runSource, mode):
    return { 'runSource':runSource,
             'mode':mode,
             'rundeck':rundeckName(runSource,mode),
             'RUN':"RUN="+rundeckName(runSource,mode)
             }

# Return a dict that specifies the configuration under test along
# with the number of processors being executed.
def newRunCase(configuration, npes):
    configuration['NPES'] = str(npes)
    return configuration

# Return additional arguments to "make" needed for building non serial configurations.
def getBuildOptions(configuration):
    mode = configuration['mode']
    if mode == 'serial':
        return ""
    elif mode == 'mpi':
        return "ESMF=YES"
    elif mode == 'openmp':
        return "MP=YES"
    
# Return additional arguments to "make setup" needed for running non serial configurations.
def getRunOptions(runCase):
    mode = runCase['mode']
    if mode == 'serial':
        return ""
    elif mode == 'mpi':
        npes = runCase['NPES']
        return "ESMF=YES NPES=" + npes
    elif mode == 'openmp':
        npes = runCase['NPES']
        return "MP=YES NPROC=" + npes

# For serial case we need to build "aux" to create the CMPE002 for comparing results.
# Otherwise we just need gcm.
def getBuildTargets(configuration):
    mode = configuration['mode']
    if mode == 'serial':
        return "gcm aux"
    else:
        return "gcm"
    
# Build a configuration.
def build(configuration):
    run = configuration['RUN']
    options = getBuildOptions(configuration)
    targets = getBuildTargets(configuration)
    try:
        systemCommand("make vclean")
        systemCommand("make " + targets + " " + run + " " + options)
    except:
        print "   Failed to build " + configuration['rundeck']
        raise
    
def run1hr(runCase):
    rundeck = runCase['rundeck']
    options = getBuildOptions(runCase)
    try:
        systemCommand("make setup_nocomp " + runCase['RUN'] + " SETUP_FLAGS=-wait " + getRunOptions(runCase))
        systemCommand("cd " + rundeck + "; cp fort.2 " + checkpointFileName(runCase, "1hr"))
    except:
        message = "   Failed to run 1 hour test for " + rundeck
        if runCase.has_key('NPES'):
            message += " on " + runcase['NPES'] + " processors."
        print message
        raise
    
def run1dy(runCase):
    rundeck = runCase['rundeck']
    options = getBuildOptions(runCase)
    expectedRC = 13; # modelE convention
    restart = "./" + rundeck
    if runCase['mode'] == 'mpi':
        restart += " -np " + runCase['NPES']
    restart += " -r"
    try:
        systemCommand("cd " + rundeck + "; cp " + checkpointFileName(runCase, "1hr") + " fort.2")
        systemCommand("cd " + rundeck + "; touch I; " + restart + "; test `head -1 run_status` -eq " + str(expectedRC))
        systemCommand("cd " + rundeck + "; cp fort.2 " + checkpointFileName(runCase, "1dy"))
    except:
        message = "  Failed to run 1 day continuation for " + rundeck
        if runCase.has_key('NPES'):
            message += " on " + runcase['NPES'] + " processors."
        print message
        raise
    

def compare(runA, runB, duration):
    try:
        rundeck = runA['rundeck']
        numLinesExpected = "4"
        serial = newConfiguration(runA['runSource'],'serial')
        file1 = rundeck + "/" + checkpointFileName(runA, duration)
        file2 = rundeck + "/" + checkpointFileName(runB, duration)
        try:
            cmp = serial['rundeck']+"_bin/CMPE002P " + file1 + " " + file2
            systemCommand("numLines=`" + cmp + "| wc -l`; test " + numLinesExpected + " -eq $numLines")
        except:
            print "   Prognostics differ in " + file1 + " and " + file2
            raise
        try:
            cmp = serial['rundeck']+"_bin/CMPE002 " + file1 + " " + file2
            systemCommand("numLines=`" + cmp + "| wc -l`; test " + numLinesExpected + " -eq $numLines")
        except:
            print "   Diagnostics differ in " + file1 + " and " + file2 + " but prognostics agree. Continuing ..."
    except:
        print "   Comparison between " + file1 + " and " + file2 + " failed."
        raise
    
for rundeck in sys.argv[1:]:
    print "Verifying rundeck " + rundeck;
    serialConfiguration = newConfiguration(rundeck, "serial")
    mpiConfiguration = newConfiguration(rundeck,"mpi")

    if os.environ.has_key('NP_LIST'):
        npList = os.environ['NP_LIST']
    else:
        npList = [1, 4]

    try:
        build(serialConfiguration)
        run1hr(serialConfiguration)
        run1dy(serialConfiguration)

        try:
            build(mpiConfiguration)
            for npes in npList:
                runCase = newRunCase(mpiConfiguration, npes)
                run1hr(runCase)
                compare(serialConfiguration, runCase, '1hr')
        except:
            print "  ... abandoning mpi configuration."

        try:
            for npes in npList:
                runCase = newRunCase(mpiConfiguration, npes)
                run1dy(runCase)
                compare(serialConfiguration, runCase, '1dy')
        except:
            print "  ... abandoning 24 hr mpi runs."
            raise

        print "  ... verification complete."

    except:
        print "  ... abandoning rundeck."

            
print "Done"
                    
                    
