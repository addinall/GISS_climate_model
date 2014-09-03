#!/usr/bin/perl
## setup_e - set up GCM modelE run
## This script is not supposed to be executed by users but rather
## is a part of Makefile functionality. It is run as a part of
##           gmake setup RUN=RunID
## Inside the Makefile it is called (from .../decks): setup_e RunID
## It looks for $HOME/.modelErc file and extracts all necessary options from
## there. If such file is not present the default options which are
## specified below will be used. Those default options are adjusted for
## working environment of ra.giss.nasa.gov.
##
## The options setup_e is looking for are:
## CMRUNDIR - directory to which all run directories will be linked.
## EXECDIR - path to directory with modelE scripts and with some executables.
## SAVEDISK - a directory (big) where all run directories will be created.
## GCMSEARCHPATH - directory to search for gcm input files.
## MAILTO - email address of the user. (default `whoami`)
## UMASK - the value of 'umask' to be used for model runs.

## default settings
$CMRUNDIR="/u/cmrun";
$EXECDIR="/u/exec";
$SAVEDISK="/raid1";
$GCMSEARCHPATH="/u/cmrun";
$NETCDFHOME="/usr/bin";
$MAILTO="";
$UMASK=002;
$MIN_STACK=0;
$omp=0;
$mpi=0;
$nproc=1;
$flag_wait=0;
$MPIDISTR="";
$debugger="";
$QSUB_STRING="";
$MPIRUN_COMMAND="";

## if $HOME/.modelErc is present get settings from there

if (exists $ENV{MODELERC}) {
  $modelerc = $ENV{MODELERC}
}
else {
  $modelerc = (getpwuid($>))[7]."/.modelErc";
}

if ( -f $modelerc ) {
    print "Using settings from $modelrc \n" ;
    open MODELERC, $modelerc or die "can't open $modelerc";
    while(<MODELERC>) {
	$CMRUNDIR = $1 if /^ *CMRUNDIR *= *(\S+)/;
	$EXECDIR = $1 if /^ *EXECDIR *= *(\S+)/;
	$SAVEDISK = $1 if /^ *SAVEDISK *= *(\S+)/;
	$GCMSEARCHPATH = $1 if /^ *GCMSEARCHPATH *= *(\S+)/;
	$MAILTO = $1 if /^ *MAILTO *= *(\S+)/;
	$UMASK = oct "$1" if /^ *UMASK *= *(\S+)/;
	$NETCDFHOME = $1 if /^ *NETCDFHOME *= *(\S+)/;
	$QSUB_STRING = $1 if /^ *QSUB_STRING *= *\"?([^ ^#][^#^"]*).*\n/;
        $MPIRUN_COMMAND = $1 if /^ *MPIRUN_COMMAND *= *\"?([^ ^#][^#^"]*).*\n/;
    }
    close MODELERC;
} else {
    print "$HOME/.modelErc is not present. Using default settings.\n";
}

if ( $MAILTO =~ /^\s*$/ ) { $MAILTO = `whoami`; chop $MAILTO; }
print "CMRUNDIR = $CMRUNDIR\n";
print "EXECDIR = $EXECDIR\n";
print "GCMSEARCHPATH = $GCMSEARCHPATH\n";
print "SAVEDISK = $SAVEDISK\n";
print "MAILTO = $MAILTO\n";
printf "UMASK = %03lo\n", $UMASK;
print "QSUB_STRING = $QSUB_STRING\n";
print "MPIRUN_COMMAND = $MPIRUN_COMMAND\n";

while ($_ = $ARGV[0], /^-/) {
    shift;
    last if /^--$/;
    if (/^-omp\b/) { $omp = 1; $nproc = shift; next;}
    if (/^-mpi\b/) { $mpi = 1; $nproc = shift; next;}
    if (/^-wait\b/) { $flag_wait = 1; next;}
    if (/^-mpidistr\b/) { $MPIDISTR = shift; next;}
    if (/^-debug\b/) { $debugger = shift; $flag_wait = 1; next;}
    print "setup: unknown option $_ \n"; exit 1;
}

if ( $omp && $mpi ) {
    print "you can't use -omp and -mpi at the same time\n";
    exit 1;
}

if ( $nproc !~ /^\d+$/ || $nproc < 1 ) {
    print "number of processors should be >= 1\n";
    print "you have: $nproc\n";
    exit 1;
}

if ( $#ARGV != 0 ) { 
    print "This script is not supposed to be run outside of Makefile\n";
    print "Inside the Makefile it is called:";
    print "setup_e.pl {-omp|-mpi} num_proc runID\n";
    print "$#ARGV\n";
    exit 1; 
}

$runID = shift;
$rfile = "$runID.R";
$RunDir = "$SAVEDISK/$runID";
$lockfile = "$RunDir/lock";
$CMRUN = $CMRUNDIR;
$DeckDir = `pwd`; chop $DeckDir;
$umask_inv = $UMASK ^ 0777;
$umask_str = sprintf "%lo", $UMASK;
$NETCDFBIN = "$NETCDFHOME/bin";
$netcdf_template_file = "$runID.nctemp";

## check if this run is already running
if ( -f $lockfile ) {
    print "            **********************                \n";
    print "$1 seems to be already running in $SAVEDISK/$runID\n";
    print "If you think it is an error, then most probably this\n";
    print "task was interrupted in an unusual way. Please check.\n";
    print "Then remove the lock file:\n";
    print "$lockfile\n";
    print "and restart the setup.\n";
    print "            **********************                \n";
    exit 1;
}

## check if the rundeck is present
if ( ! -s $rfile ) {
   print "File does not exist: $rfile\n";
   exit 1;
}

print "setting up run $runID\n";
print "output files will be saved in $SAVEDISK/$runID\n";

## Create the run directory if necessary
if ( ! -d $RunDir ) {
    mkdir $RunDir, 0777 & $umask_inv 
	or die "Can't create $RunDir. Aborting setup.\n";
}

## Check that link is not already correct (added by gavin)
if ( -e $runID ) {
    if ( `ls -l $runID` !~ /-> *$RunDir$/ ) {
	print "./$runID exists and is pointing to something else:\n";
	print "Please check. Aborting setup.\n";
	exit 1;
    }
} else {
    symlink $RunDir, $runID or die "Can't create link $runID in local dir";
}

## Also link to $CMRUN if different from $SAVEDISK (added by gavin)
if ( $CMRUN ne $SAVEDISK ) {
    if ( -e "$CMRUN/$runID" ) {
	die "Can't create a link $CMRUN/$runID" if ! -l "$CMRUN/$runID";
	if (`ls -l $CMRUN/$runID` !~ /-> *$RunDir$/ ) {
	    unlink "$CMRUN/$runID" or die "Can't rm old link $CMRUN/$runID";
	    symlink "$RunDir", "$CMRUN/$runID" or die "Can't create link";
	}
    } else {
        symlink "$RunDir", "$CMRUN/$runID" or die "Can't create link from $RunDir to $CMRUN/$runID";
    }
}

## Make sure that we really can write to dir/link $runID
if ( ! ( -d $runID && -r $runID && -w $runID && -x $runID ) ) {
    print "Couldn't create the link $runID\n";
    print "or run directory has wrong permissions. Aborting setup.\n";
    exit 1;
}

## If NetCDF template file is present copy it to run directory
if ( -s $netcdf_template_file ) {
    `cp $netcdf_template_file $runID/nctemp`;
}

## Switching to Run directory
chdir "$runID" or die "Can't chdir to $runID : $!\n";

open PRT, ">$runID".".PRT" or die "can't open ${runID}.PRT for writing\n";
print PRT "0Run $runID\n";

## Move executable from "$runID"_bin directory to $runID
if ( ! -f "$DeckDir/$runID"."_bin/$runID.exe" ) {
    print "$DeckDir/$runID"."_bin/$runID.exe","not found\n";
    exit 1;
}
`mv -f "$DeckDir/$runID"_bin/$runID.exe $runID.exe`;
chmod 0777 & $umask_inv, "$runID.exe";

open RFILE, "$DeckDir/$rfile" or die "can't open $DeckDir/$rfile";

## Check signature
$_ = <RFILE>;
if ( ! /^\s*$runID/ ) {
    print "inconsistent Naming: $rfile is not $_\n" ;
    exit 1;
}
print PRT " $_";

## Read rfile until "Data input files:" is encountered
$run_options = 0;
while (<RFILE>) {
    print PRT " $_";
    if ( /^\s*Run *Options/ ) { $run_options = 1; }
    if ( /^\s*STACKSIZE *= *(\d+)/ && $run_options ) { $MIN_STACK = $1; }
    if ( /Data *input *files/ ) { last; }
}

print "Requested Stack  = $MIN_STACK\n" if ( $MIN_STACK ) ;

## Use the subsequent information to create the link and unlink files
while (<RFILE>) {
    if ( /Label *and *Namelist/ ) { last; }
    print PRT " $_";
    s/!.*//;
    push @data_files, /([\w.+_-]+\s*=\s*[\w.\/+_-]+)/g;
}

open RUNIDLN, ">$runID"."ln" or die "can't open ${runID}ln for writing\n";
open RUNIDULN, ">$runID"."uln" or die "can't open ${runID}ln for writing\n";

$flag_missing_data_files = 0;

foreach $_ ( @data_files ) {
    ($name, $dest) = split /\s*=\s*/;
    if ( $dest !~ /^\// ) { $dest = "$GCMSEARCHPATH/$dest"; }
    if ( ! -e "$dest" ) {
	print "$dest not found\n";
	$flag_missing_data_files = 1;
	#exit 1;
    }
    if ( $name !~ /^(AIC|OIC|GIC)$/ ) {
	print RUNIDLN "ln -fs $dest $name\n";
	print RUNIDULN "rm $name\n";
    } else {
	`ln -sf $dest $name` ; 
	print "using $dest for IC only\n";
    }
}

print RUNIDULN "chmod -R ugo+r . 2> /dev/null\n";

if ( $flag_missing_data_files ) { exit 1; }


close RUNIDLN;
close RUNIDULN;

open RUNTIMEOPTS, ">runtime_opts" or die "can't create runtime_opts\n";
## Architecture-dependent settings
## Linux setting is only for Intel 7 compilers, but harmless otherwise
$uname = `uname`;
if ( $uname =~ /IRIX64/ ) {
    print RUNTIMEOPTS "export PAGESIZE_DATA=64 PAGESIZE_STACK=64\n";
} elsif ( $uname =~ /AIX/ ) {
    print RUNTIMEOPTS "export XLFRTEOPTS=NAMELIST=OLD\n";
} elsif ( $uname =~ /Linux/ ) {
    print RUNTIMEOPTS "export F_UFMTENDIAN=big\n";
    print RUNTIMEOPTS "export G95_ENDIAN=BIG\n";   #needed for g95
}
print RUNTIMEOPTS <<EOF;
    if [ `ulimit -s` != 'unlimited' ] ; then
      if [ `ulimit -s` -lt $MIN_STACK ] ; then
        ulimit -s $MIN_STACK || \\
          echo "!!! Could not set required stack size !!!" ;
        echo "current stack: `ulimit -s`"
      fi
    fi
EOF
close RUNTIMEOPTS;

## Architecture-dependent commands to start openMP/MPI
$omp_run = '';
$mpi_run = 'echo no support for MPI, will not run ';
if ( $uname =~ /IRIX64/ ) {
    $omp_run = "export MP_SET_NUMTHREADS=\$NP; ";
    $mpi_run = "mpirun -np \$NP ";
} elsif ( $uname =~ /OSF1/ ) {
    $omp_run = "set OMP_SET_NUMTHREADS=\$NP; ";
    $mpi_run = "prun -s -n \$NP ";
} elsif ( $uname =~ /AIX/ ) {
    $omp_run = "export MP_SET_NUMTHREADS=\$NP; "; #?? check
} elsif ( $uname =~ /Linux/ ) {
    $omp_run = "export OMP_NUM_THREADS=\$NP; ";
    if ( $MPIDISTR =~ /mvapich2/ ) {
        $mpi_run = "mpirun_rsh -np \$NP -hostfile \\\$PBS_NODEFILE ";
    } else {
        $mpi_run = "mpirun \$MPI_FLAGS -np \$NP ";
    }
    if ( $MPIDISTR =~ /SCALI/ ) {
	$mpi_run .= "-inherit_limits ";
    }
    if ( $MPIDISTR =~ /openmpi/ ) {
        $mpi_run .= "--mca btl_openib_warn_no_hca_params_found 0 ";
    }
} elsif ( $uname =~ /Darwin/ ) {
    $omp_run = "export OMP_NUM_THREADS=\$NP; ";
    $mpi_run = "mpirun \$MPI_FLAGS -np \$NP ";
}

# overwrite the command with the one from modelErc if present
if ( $MPIRUN_COMMAND ) {
    $mpi_run = $MPIRUN_COMMAND." ";
}


chmod 0777 & $umask_inv, "${runID}ln", "${runID}uln", "runtime_opts";

open I, ">I" or die "can't open 'I' for writing\n";
## Check signature
$_ = <RFILE>;
if ( ! /^\s*$runID[ (]/ ) {
    print "inconsistent Naming: $rfile is not start of $_\n" ;
    exit 1;
}
print I;
$_ = <RFILE>;
print I;
while (<RFILE>) {
    chop;
    s/!.*//g;
    s/^\s*/ /;
    s/\s*$//;
    next if ( ! $_ ); 
    print I "$_\n";
}
close I;
close PRT;
chmod 0666 & $umask_inv, "I", "${runID}.PRT";

## Run the 1st hour then clean up
print "starting the execution \n";
print "current dir is ", `pwd`;

## Switching to background
if ( (!$flag_wait) && ($pid = fork) ) {
    print "Starting 1st hour in the background.\n\n";
    exit 0;
}

if ( ! defined $pid ) {
    print "Fork failed. Continuing in foreground.\n";
}

if ( $mpi ) {
    $run_command = $mpi_run;
    $qsub_command = $QSUB_STRING;
} else {
    $run_command = $omp_run; # serial also fits here
    $qsub_command = "";
}

## a hack to run setup in debugger

if ( $debugger ) {
print <<`EOC`;
    umask $umask_str
    touch lock
    . ./runtime_opts
    echo "-99" > run_status
    echo "INPUT not yet completed" >> run_status
    ./"$runID"ln
EOC
exec "$debugger ./$runID.exe; rm -f AIC GIC OIC; ./${runID}uln; rm -f lock";
}

## If on some machines MPI can't be used interactively, a hack
## can be intruduced here to run 1st hour in serial mode


## hack to set "ulimit" for openmpi jobs
# seems it is not needed any more (and sometimes causes problems)
#$stack_in_wrapper = ($MIN_STACK>1024000) ? $MIN_STACK : 1024000;
#if ( $MPIDISTR =~ /openmpi/ ) {
#print <<`EOC`;
#    echo "#!/bin/sh" > ${runID}.wrapper
#    echo "ulimit -s $stack_in_wrapper" >> ${runID}.wrapper
#    echo `pwd`/"$runID.exe \\\$\@" >> ${runID}.wrapper
#EOC
#} else {
#    `ln -sf $runID.exe $runID.wrapper`;
#}
#`chmod 755 $runID.wrapper`;

## Running the model
print <<`EOC`;
    umask $umask_str
    touch lock
    . ./runtime_opts
    echo "-99" > run_status
    echo "INPUT not yet completed" >> run_status
    if [ -s nctemp ] ; then $NETCDFBIN/ncgen -b -o nctemp.nc nctemp ; fi
    ./"$runID"ln
    NP=$nproc
    echo "#!/bin/sh" > setup_command
    echo cd `pwd` >> setup_command
    echo "$run_command ./"$runID".exe -i I >> ${runID}.PRT" >> setup_command
    chmod 755 setup_command
    $qsub_command ./setup_command
    touch run_status  # don't trust clocks on batch nodes
    rc=`head -1 run_status`
    rm -f AIC GIC OIC
    ./"$runID"uln
    rm -f lock
    exit \$rc
EOC

# print warning if any
open(PRTFILE, "$runID.PRT") or die "can't open $runID.PRT\n";
while(<PRTFILE>) {
    if ( /^ ?WARNING/ ) { print ; }
}
	close PRTFILE;

$rcode = $? >> 8;
# mpirun returns 0 on success...
if ( $rcode != 13 && $rcode != 12 ) {
    print " Problem encountered while running hour 1 :\n";
    if ( $rcode != 1 ) {
	$error_message = `tail -1 run_status`; chop $error_message;
    } else {
	$error_message = "Unknown reason (Segmentation fault?)";
    }
    print " >>> $error_message <<<\n";
    exit 4 ;
} else {
    print "1st hour completed successfully\n";
}

## HACK needed for proper setup of mvapich2 runs
#  (need to remove extra '\')
$run_command =~ s/\\(?=\$PBS_NODEFILE)//g;

## Create executable script file RUNID
open RUNID, ">$runID" or die "can't open $runID for writing\n";
print RUNID <<EOF;
\#!/bin/sh
    PRTFILE=${runID}.PRT
    IFILE="I"
    RESTART=0
    NP="\$MP_SET_NUMTHREADS"
    if [ "\$NP"x = x ] ; then NP=1; fi
    while [ \$\# -ge 1 ] ; do
      OPT=\$1 ; shift
      case \$OPT in
        -q)
            PRTFILE='/dev/null'
            ;;
        -l)
            PRTFILE="\$1" ; shift
            ;;
        -i)
            IFILE="\$1" ; shift
            ;;
        -r)
            RESTART=1
            ;;
        -np)
            NP="\$1" ; shift
            ;;
         *)
            echo "Warning: wrong option ignored: \$OPT"
            ;;
      esac
    done
    umask $umask_str
    if [ `head -1 run_status` -eq 13 ] ; then
      if [ `find run_status -newer I` ] ; then
        echo 'run seems to have finished'
      exit 1; fi; fi
    if [ -f lock ] ; then
      echo 'lock file present - aborting' ; exit 1 ; fi
    touch lock
    echo '-99' > run_status
    echo 'INPUT not yet completed' >> run_status
    . ./runtime_opts
    if [ -s nctemp ] ; then $NETCDFBIN/ncgen -b -o nctemp.nc nctemp ; fi
    ./${runID}ln
    $run_command ./${runID}.exe -i ./\$IFILE > \$PRTFILE
    rc=`head -1 run_status`
    ./${runID}uln
    rm -f lock
    if [ \$RESTART -ge 1 ] ; then exit \$rc ; fi
    exec $EXECDIR/runpmE $runID \$rc $MAILTO
EOF
close RUNID;
chmod 0777 & $umask_inv, $runID;
## end of RUNID script
`ln -sf $runID E`;

## overwrite RUNID.I omitting ISTART=.. line.
`umask $umask_str; grep -v "ISTART=" I > I.tmp; mv -f I.tmp I`;

## setup finished normally
exit 0 ;

