#!/bin/bash
#PBS -l select=1:ncpus=4:proc=wood
#PBS -l walltime=4:00:00
#PBS -W group_list=a940a
#PBS -N modelE_test
# #PBS -m "abe"
# #PBS -M "Thomas.L.Clune@nasa.gov"

#BSUB -n 4
#BSUB -W 06:00
#BSUB -q general_lng
#BSUB -P k3002
#BSUB -J modelE_test
#BSUB -B
# #BSUB -u "Thomas.L.Clune@nasa.gov"

# Script to build standard modelE rundecks under OpenMP and MPI for the purposes of
# detecting defects in paralellization.

RUNDECKS="E1M20 E1oM20 E1F20 E001tr"
case `hostname` in
    discover* | borg* )
        export BASELIBDIR=/usr/local/other/baselibs/ESMF222rp3_NetCDF362b6_10.1.017_intelmpi/Linux
        export MODES="serial other MPI OpenMP"
        export NETCDFHOME=
        export MPIDISTR=intel
	export USE_COMPILER=intel
	export COMPILER_VERSION=10
	export GCMSEARCHPATH=/discover/nobackup/projects/giss/prod_input_files;;
    *) 
        export NOBACKUP=$HOME
	export BASELIBDIR=/home/trayanov/baselibs/v2_2rp2_nb2/LinuxIA64
        export MODES="serial MPI OpenMP";;
esac


export MPI_NPROCS="1 4"
export OMP_NPROCS="1 4"
export SCRATCH=$NOBACKUP/modelE_scratch
export RESULTS=$NOBACKUP/modelE_regression
export BASELINE=$NOBACKUP/modelE_baseline
export MODELE_DATA=$NOBACKUP/modelE_data
export SAVEDISK=$SCRATCH/savedisk
export STATUS=0

function setNumLines {
    local rundeck=$1
    case $rundeck in
	E1M20 | E1F20)  NUM_LINES=110;;
        E1oM20) NUM_LINES=139;;
	E001tr) NUM_LINES=123;;
    esac
}

# build a rundeck in a particular configuration
function buildRundeck {
    local rundeck=$1
    local mode=$2
    local run=`buildName $rundeck $mode`
    local logfile=$3

    logMessage "    ... building $run."

    gmake rundeck RUN=$run RUNSRC=$rundeck          >> $logfile 2>&1

    case $mode in
	    MPI)    EXTRAS="ESMF=YES NPES=1" ;;
	    OpenMP) EXTRAS="MP=YES NPROC=1 EXTRA_FFLAGS=-mp" ;;
	    other) EXTRAS="EXTRA_FFLAGS=-mp" ;;
    esac

    gmake gcm RUN=$run $EXTRAS SETUP_FLAGS=-wait  >> $logfile 2>&1

    if [ "$mode" == "serial" ]; then 
	gmake aux RUN=$run     >> $logfile 2>&1
    fi
    logMessage "    ... done building $run."
}

function environmentModules {
    source ${MODULESHOME}/init/bash
    module purge
    case `hostname` in
	discover* | borg* )  module load comp/intel-10.1.023 lib/mkl-10.1.2.024 mpi/impi-3.2.1.009
    esac
    echo "****************"
    echo " module list: "
    module list | echo
    echo $USE_COMPILER
    echo "****************
    export
"
}

function setEnvironment {
    environmentModules
}


function checkoutCode {
    local PRIVATE_DIRECTORY=$1
    local logfile=$2
    export CVS_RSH=ssh
    export CVSROOT=simplex.giss.nasa.gov:/giss/cvsroot
    cvs co -d $PRIVATE_DIRECTORY modelE >> $logfile
}

createModelErc() {
    cat > ~/.modelErc <<EOF
DECKS_REPOSITORY=$MODELE_DATA/decksRepository
CMRUNDIR=$MODELE_DATA/cmrun
EXECDIR=$MODELE_DATA/exec
OVERWRITE=YES
OUTPUT_TO_FILES=YES
VERBOSE_OUTPUT=YES
MP=NO
UMASK=002
SAVEDISK=$SAVEDISK
GCMSEARCHPATH=$GCMSEARCHPATH
NETCDFHOME=
COMPILER=$USE_COMPILER
EOF

    source ~/.modelErc
}

makeOutputDirectories() {
    echo "SAVEDISK is $SAVEDISK"
    mkdir -p $MODELE_DATA $DECKS_REPOSITORY $CMRUNDIR $EXECDIR $SAVEDISK >> $MAIN_LOG 2>&1
}

function setUp {
    mkdir -p $RESULTS
#
#   Remove previous log files if any
#
    rm -f $RESULTS/*.log

    export MAIN_LOG=$RESULTS/main.log
    rm -f $MAIN_LOG
    touch $MAIN_LOG
    logMessage "setUp"
    
    setEnvironment
    createModelErc
    makeOutputDirectories
}


function cleanUp {
    logMessage "cleanUp"
    # delete source directories
    rm -rf $MODELE_DATA
    rm -rf ~/.modelErc
    rm -rf $SCRATCH
}

function runTest {
    local rundeck=$1
    local mode=$2
    local logfile=$3
    local run=`buildName $rundeck $mode`

    case $mode in
	serial) runSerial $rundeck $logfile;;
        MPI)    runMPI $rundeck $logfile;;
	OpenMP) runOpenMP $rundeck $logfile;;
	other) runOther $rundeck $logfile;;
    esac
}

function compare {
    local rundeck=$1
    local mode=$2
    local file1=$3
    local file2=$4
    local logfile=$5

    local run=`buildName $rundeck $mode`
    local CMP=$SCRATCH/$rundeck.serial/decks/$rundeck.serial_bin/CMPE002 
    
    ln -s $file1 $run.f1
    ln -s $file2 $run.f2

    local numLinesFound=`$CMP $run.f1 $run.f2 | wc -l`
    local minLinesExpected=`$CMP $run.f1 $run.f1 | wc -l`
    local maxLinesExpected

    local rc

    if [ $mode == OpenMP ]; then 
	let "maxLinesExpected=$minLinesExpected+22"
    else
	maxLinesExpected=$minLinesExpected
    fi

    if [[ ($numLinesFound -le $maxLinesExpected) && ($minLinesExpected -le $numLinesFound) ]]; then
	rc=0
    else
	local numArgs=$#
	logError "Comparison failed for files $file1 and $file2."
	if [ $numArgs == 5 ]; then
	    echo "Numline discrepancy. Expected $maxLinesExpected but found $numLinesFound." >> $logfile
	    $CMP $run.f1 $run.f2 >> $logfile
	else
	    logMessage "Numline discrepancy. Expected $maxLinesExpected but found $numLinesFound."
	    $CMP $run.f1 $run.f2 >> $MAIN_LOG
	fi
	rc=1
    fi
    rm $run.f1 $run.f2
    return $rc
}

function runSerial {
    local rundeck=$1
    local logfile=$2
    local mode=serial
    local run=`buildName $rundeck $mode`

    # serial
    logMessage "   ... serial 1 hour run ..."
    gmake setup_nocomp RUN=$run SETUP_FLAGS=-wait  >> $logfile 2>&1
    if [ $? == 0 ]; then
	logMessage '   ... 1 hour serial run has completed.'
	cp $run/fort.2 $RESULTS/$run.1hr
    else
	logError $rundeck $mode "     ... 1 hour serial run has failed.  Aborting"
	return
    fi;

    logMessage "   ... serial 1 day continuation."
    pushd $CMRUNDIR/$run
    local NORMAL_TERMINATION=13
    ./$run -r >> $logfile 2>&1
    if [ $? == $NORMAL_TERMINATION ]; then
	echo '... 1 day serial continuation has completed.' >> $MAIN_LOG
	cp fort.2 $RESULTS/$run.1dy
    else
	logError $rundeck $mode "     ... 1 day serial run has failed.  Aborting"
	popd
	return
    fi;
    popd
}

#
# The following function is a kludge to support the fact that on the Altix,
# OpenMP results only match serial results if both are compiled with "-mp"
#
function runOther {
    local rundeck=$1
    local logfile=$2
    local mode=other
    local run=`buildName $rundeck $mode`

    # serial
    logMessage "   ... other 1 hour run ..."
    gmake setup_nocomp RUN=$run SETUP_FLAGS=-wait  >> $logfile 2>&1
    if [ $? == 0 ]; then
	logMessage '   ... 1 hour other run has completed.'
	cp $run/fort.2 $RESULTS/$run.1hr
    else
	logError $rundeck $mode "     ... 1 hour other run has failed.  Aborting"
	return
    fi;

    logMessage "   ... other 1 day continuation."
    pushd $CMRUNDIR/$run
    local NORMAL_TERMINATION=13
    ./$run -r >> $logfile 2>&1
    if [ $? == $NORMAL_TERMINATION ]; then
	echo '... 1 day other continuation has completed.' >> $MAIN_LOG
	cp fort.2 $RESULTS/$run.1dy
    else
	logError $rundeck $mode "     ... 1 day other run has failed.  Aborting"
	popd
	return
    fi;
    popd
}

function runMPI {
    local rundeck=$1
    local mode=MPI
    local logfile=$2
    local run=`buildName $rundeck $mode`

    logMessage "   ... MPI tests ..."
    local CMP=$SCRATCH/$rundeck.serial/decks/$rundeck.serial_bin/CMPE002
    local NORMAL_TERMINATION=13

    for np in $MPI_NPROCS; do
        gmake setup_nocomp RUN=$run ESMF=YES NPES=$np SETUP_FLAGS=-wait >> $logfile 2>&1
	compare $rundeck $mode $RESULTS/$rundeck.serial.1hr $run/fort.2 $logfile
	if [ $? == 0 ]; then
	    logMessage "   ... 1 hour $mode run on $np processes matches serial results."
	else
	    logError $rundeck $mode "MPI 1 hour run on $np processes does not match serial results."
	    STATUS=1
	    return 1
	fi
	pushd $CMRUNDIR/$run
	./$run -np $np -r >> $logfile 2>&1
	if [ $? == $NORMAL_TERMINATION ]; then
	    logMessage "   ... 1 day $mode continuation on $np processes has completed on $np processes."
	else
	    logError $rundeck $mode "1 day continuation on $np processes failed."
	    STATUS=1
	    popd
	    return 1
	fi
	popd
	compare $rundeck $mode $RESULTS/$rundeck.serial.1dy $run/fort.2 $logfile
	if [ $? == 0 ]; then
	    logMessage "   ... 1 day $mode continuation on $np processes matches serial results."
	else
	    logError $rundeck $mode "1 day run on $np processes does not match serial results."
	    STATUS=1
	    return 1
	fi
    done
}

function runOpenMP {
    local rundeck=$1
    local mode=OpenMP
    local logfile=$2
    local run=`buildName $rundeck $mode`

    logMessage "   ... OpenMP tests ..."
    local CMP=$SCRATCH/$rundeck.serial/decks/$rundeck.serial_bin/CMPE002
    local NORMAL_TERMINATION=13

    for np in $OMP_NPROCS; do
        gmake setup_nocomp RUN=$run MP=YES NPROC=$np SETUP_FLAGS=-wait >> $logfile 2>&1
	local referenceMode
	case `hostname` in
	    palm | explore* ) referenceMode=other;;
            *) referenceMode=serial;;
	esac

	compare $rundeck $mode $RESULTS/$rundeck.$referenceMode.1hr $run/fort.2 $logfile
	if [ $? == 0 ]; then
	    logMessage "   ... 1 hour $mode run on $np threads matches serial* results.."
	else
	    logError $logfile $mode "1 hour run on $np threads does not match serial* results."
	    STATUS=1
	    return 1
	fi
	pushd $CMRUNDIR/$run
	./$run -np $np -r >> $logfile 2>&1
	if [ $? == $NORMAL_TERMINATION ]; then
	    logMessage "   ... 1 day $mode run on $np threads has completed."
	else
	    logError $rundeck $mode '1 day continuation on $np threads failed.' >> $MAIN_LOG
	    popd
	    return 1
	fi
	popd
	compare $rundeck $mode $RESULTS/$rundeck.$referenceMode.1dy $run/fort.2 $logfile
	if [ $? == 0 ]; then
	    logMessage "   ... 1 day $mode continuation on $np threads matches serial results.."
	else
	    logError $rundeck $mode "1 day continuation on $np threads does not match serial results."
	    STATUS=1
	    return 1
	fi
    done
}


function checkBaseline {
    local rundeck=$1
    local CMP=$SCRATCH/$rundeck.serial/decks/$rundeck.serial_bin/CMPE002

    # if _all_ runs complete _and_ _all_ checks pass
    for duration in 1hr 1dy; do
	compare $rundeck serial $RESULTS/$rundeck.serial.$duration $BASELINE/$rundeck.serial.$duration $MAIN_LOG
	if [ $? == 0 ]; then
	    logMessage "   ... baseline $rundeck $duration run matches serial results."
	else
	    logMessage "   ... baseline $duration run does not match serial results." >> $MAIN_LOG
	    if [ $STATUS == 0 ]; then
		logMessage "   ... All consistency checks passed => Assuming intentional change to results!"
		cp $RESULTS/$rundeck.serial.$duration $BASELINE/$rundeck.serial.$duration
	    else
		logMessage "   ... Not updating baseline results due to failed consistency check."
	    fi
	fi
    done
}

function process {
    local rundeck=$1
    local mode=$2
    local logfile=$3
    local run=`buildName $rundeck $mode`

    
      echo "    ... cvs co:" $rundeck $mode $logfile >> $logfile
    checkoutCode $run $logfile
      echo "    ... build:" $rundeck $mode $logfile  >> $logfile
    pushd $run/decks
    buildRundeck $rundeck $mode $logfile
      echo "    ... run and test:" $rundeck $mode    >> $logfile
    runTest $rundeck $mode $logfile
    popd
      echo "    ... done:" $rundeck $mode            >> $logfile
}

function logMessage {
    local message=$1
    echo  $message >> $MAIN_LOG
}

function logError {
    local rundeck=$1
    local mode=$2
    local message=$3

    echo " "                                            >> $MAIN_LOG
    echo "********************************************" >> $MAIN_LOG
    echo "Error in rundeck: " $rundeck "  mode: " $mode >> $MAIN_LOG
    echo "     " $message                               >> $MAIN_LOG
    echo "********************************************" >> $MAIN_LOG
    echo " "                                            >> $MAIN_LOG
}

function clearStatus {
    STATUS=0
}

# Build a string to be used as a unique identifier for
# a given build from rundeck and mode.
#
function buildName {
    local rundeck=$1
    local mode=$2
    echo $1.$2
}    

setUp

echo "Hello"

mkdir $SCRATCH
pushd $SCRATCH

logMessage "Regression testing on platform $HOSTNAME."
logMessage "   Processing rundecks $RUNDECKS."
logMessage "   Processing modes $MODES."
logMessage " "

for rundeck in $RUNDECKS; do
    logMessage "Beginning processing for rundeck $rundeck ..."
    clearStatus

    for mode in $MODES; do
	run=`buildName $rundeck $mode`
	logfile=$RESULTS/$run.log
	checkoutCode $run $logfile
	pushd $run/decks
	buildRundeck $rundeck $mode $logfile &
	popd
    done
    
    wait

    setNumLines $rundeck
    
    for mode in $MODES; do
	logMessage "    ... running tests for $mode mode."
	run=`buildName $rundeck $mode`
	logfile=$RESULTS/$run.log
	pushd $run/decks
	runTest $rundeck $mode $logfile
	popd
	logMessage "    ... done running tests for $mode mode."
    done

    logMessage "Comparing $rundeck results against archive."
    checkBaseline $rundeck

    logMessage "Completed processing for rundeck $rundeck."
    logMessage " "
done

popd
#cleanUp
logMessage "All processing completed."

# mail log to standard mailing list
export MAILING_ADDRESS=giss-modelE-regression@sourcemotel.gsfc.nasa.gov
export SUBJECT="modelE regression tests on $HOSTNAME"
case `uname` in
    OSF1)  export MAIL=Mail;;
    Linux) export MAIL=mail;;
esac
$MAIL -s "$SUBJECT" $MAILING_ADDRESS < $MAIN_LOG

