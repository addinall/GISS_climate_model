LATEST GISS CLIMATE MODEL
=========================

### Cloned - Mark Addinall

Instead of just dis-beleiving the climate models from afar, I have been studying them.  This is the 'latest' offering from GISS.  Most of the code is VERY VERY OLD and a few parameters have been tarted up.

I am going to re-engineer this nonsense into a modern language and then all can see and point fingers in mirth.

What I want to know is where all the R&D money went?

----
  
NOTE: The full documentation on how to run the model is in doc/HOWTO.html. 
This summary is not complete, and not necessarily up-to-date either.

PLEASE READ THE FULL DOCUMENTATION - IT REALLY WILL MAKE YOUR LIFE EASIER!

The directory tree of the modelE has the following structure:

  modelE
        |
        |-/model   (the source code for GCM model)
        |
        |-/aux     (auxiliary programs such as pre- and post-processing)
        |
        |-/exec    (various scripts needed to compile and setup the model)
        |
        |-/doc     (directory for documentation)
        |
        |-/decks   (directory for rundecks)
                |
                |-<run_name_1>.R     (rundeck for the run <run_name_1>)
                |
                |-/<run_name_1>_bin  (directory for binaries for <run_name_1>)
                |
                |-/<run_name_1>      (link to directory where you setup
                |                     and run <run_name_1>)
                |-<run_name_2>.R
                ................

           Configuring the model on your local system

   Intended working directory is directory modelE/decks. In this
directory type "gmake config". This will create a default ~/.modelErc
file in your home directory. This should be edited so that run output,
rundeck libraries etc. can be properly directed, and so that the
compile options (multiple processing, compiler name , NetCDF libraries
etc.)  can be set appropriately.

               Compiling and running the model.

   All rundecks should be created inside this directory and all "make"
commands should be run from there. The following is a typical example
of how to compile and setup a run with the name "my_run":

      cd decks                      # go to directory decks
      gmake rundeck RUN=my_run      # create rundeck for "my_run"

You will need to edit the rundeck in order to choose a configuration
that is appropriate. Once that is done...

      gmake gcm RUN=my_run          # compile the model for "my_run"
      gmake setup RUN=my_run        # run setup script for "my_run"

Make sure that you create the rundeck with "gmake rundeck ..." before
running any other commands for this run, otherwise the Makefile will
not understand you. You can skip "gmake gcm ..." and just do "gmake setup..."
in which case gcm will be compiled automatically.
Another command you want to run (after creating the rundeck) is:

      gmake aux RUN=my_run

This will compile auxiliary programs in /aux. All the binaries (both for
model and for auxiliary programs) are put into /decks/my_run.bin .

The following is a list of targets currently supported by Makefile:

 config  - copy the default .modelErc setup to your home directory.
 rundeck - create new rundeck
 depend  - create dependencies for specified rundeck
 gcm     - compile object files and build executable for specified rundeck
 aux     - compile standard auxiliary programs
 auxqflux- compile auxiliary programs for computing qflux
 auxdeep - compile auxiliary programs for setting deep ocean
 setup   - do setup for specified rundeck
 clean   - remove listings and error files
 vclean  - remove object files, .mod files and dependencies
 newstart- remove all files in the run directory
 exe     - compile gcm and put executable into RUN directory
 htmldoc - create web-based documentation for this RUN

If you run "gmake" without arguments it will print a short help.


