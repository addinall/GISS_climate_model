#package RegressionTools;
use CommandEntry;

my $standardFlags = "SETUP_FLAGS=-wait";
my $HYCOM = "E1fzhyc";

my $extraFlags;
$extraFlags{SERIAL}   = "";
$extraFlags{MPI}      = "ESMF=YES NPES=\$npes";
$extraFlags{OPENMP}   = "EXTRA_FFLAGS=-mp MP=YES NPROCS=\$npes";
$extraFlags{SERIALMP} = "EXTRA_FFLAGS=-mp"; # Intel kludge for matching serial and OpenMP.

$extraFlags{E1M20} ="";
$extraFlags{E1F20} ="";
$extraFlags{E1oM20} ="";
$extraFlags{E001tr} ="";
$extraFlags{$HYCOM} ="EXTRA_FFLAGS+=-DCHECK_OCEAN_HYCOM";

my $numLinesCMP;
$numLinesCMP{E1M20}  = 120;
$numLinesCMP{E1F20}  = 120;
$numLinesCMP{E1oM20} = 110;
$numLinesCMP{E001tr} = 123;
$numLinesCMP{$HYCOM} = 138;

sub setModuleEnvironment {
#    require perlModule;
    require "$ENV{MODULESHOME}/init/perl";
    module (purge);
    module (load, "comp/intel-10.1.023", "lib/mkl-10.1.2.024", "mpi/impi-3.2.1.009");
}

sub createTemporaryCopy {
  my $env = shift;
  my $tempDir = shift;

  my $referenceDir = $env -> {REFERENCE_DIRECTORY};
  my $commandString = "mkdir -p $tempDir;  cp -r -u $referenceDir/* $tempDir;";
  return (CommandEntry -> new({COMMAND => $commandString}));

  }

# Build rundeck and compile
# Do a vclean just in case receiving a dirty copy.
# Creat CMPE002P if SERIAL configuration - used for verification
sub compileRundeck {
  my $env = shift;
  my $installDir = shift;
  my $compiler = shift;

  my $rundeck    = $env -> {RUNDECK};
  my $configuration = $env -> {CONFIGURATION};
  my $resultsDir = $env -> {RESULTS_DIRECTORY};
  $resultsDir .="/$compiler";
  
  my $flags = "$standardFlags $extraFlags{$configuration} $extraFlags{$rundeck}";
  $flags =~ s/(\$npes)/1/eeg;

  my $expName;
  if (@_) {$expName = shift}
  else {$expName = "$rundeck.$configuration.$compiler"};

  my $logFile = "$resultsDir/$expName.buildlog";
  unlink($logFile); # delete it

  my $commandString = <<EOF;
  export MODELERC;
  echo "MODELERC is $MODELERC";
  echo "CONTENTS: ";
  cat $MODELERC;
  cd $installDir/decks;
  make rundeck RUN=$expName RUNSRC=$rundeck;
  make vclean RUN=$expName;
  make -j gcm RUN=$expName $flags COMPILER=$compiler;
EOF

  my $binDir = $expName . "_bin";
  if ($configuration eq SERIAL) {
    $commandString .= "make aux RUN=$expName $flags COMPILER=$compiler;\n";
    $commandString .= "cp $binDir/CMPE002P $resultsDir/CMPE002P.$rundeck;\n";
    $commandString .= "cp $binDir/CMPE002 $resultsDir/CMPE002.$rundeck;\n";
  }

  return (CommandEntry -> new({COMMAND => $commandString, QUEUE => "", STDOUT_LOG_FILE => "$logFile"}));
}

sub runConfiguration {
    my $env = shift;
    my $installDir = shift;
    my $npes = shift;
    my $compiler = shift;

    my $rundeck    = $env -> {RUNDECK};
    my $configuration = $env -> {CONFIGURATION};
    my $resultsDir = $env -> {RESULTS_DIRECTORY};
    $resultsDir .="/$compiler";

    my $flags = "$standardFlags $extraFlags{$configuration}";
    $flags =~ s/(\$npes)/$npes/eeg;
    my $expName = "$rundeck.$configuration.$compiler";

    my $suffix;
    my $logFile = "$resultsDir/$expName.runlog";

    if ($configuration eq "SERIAL" or $configuration eq "SERIALMP") {
	$suffix = "";
    }
    else {
	$suffix = ".np=$npes";
    }
    
    my $continue1Day = "./$expName -r ";
    if ($configuration eq "MPI" or $configuration eq "OPENMP") {$continue1Day = "./$expName -np $npes -r ";}

    my $commandString = <<EOF;
    echo "Using flags: $flags "
    export MODELERC;
    cd $installDir/decks;
    rm -f $expName/fort.2
    make setup_nocomp RUN=$expName $flags COMPILER=$compiler;
    cd $expName;
    # Remove old result so that we notice if things fail really badly
    rm -f $resultsDir/$expName.1hr$suffix;
    cp fort.2 $resultsDir/$expName.1hr$suffix;
    $continue1Day;
    # Remove old result so that we notice if things fail really badly
    rm -f $resultsDir/$expName.1dy$suffix;
    cp fort.2 $resultsDir/$expName.1dy$suffix;
EOF
  return (CommandEntry -> new({COMMAND => $commandString, QUEUE => "", STDOUT_LOG_FILE => "$logFile", NUM_PROCS => $npes}));
}

sub runSetup {
  my $env = shift;

  my $installDir    = $env -> {INSTALL_DIR};
  my $rundeck       = $env -> {RUNDECK};
  my $configuration = $env -> {CONFIGURATION};
  my $npes          = $env -> {NUM_PROCS};
  my $resultsDir    = $env -> {RESULTS_DIR};
  $resultsDir .="/$compiler";

  my $expName;
  if (@_) {$expName = shift}
  else {$expName = "$rundeck.$configuration.$compiler"};

  my $commandString = <<EOF;
cd $dir/decks;
make setup_nocompile RUN=$expName $standardFlags;
cp $expName/fort.2 $resultsDir/$expName.1hr;
pushd $expName;
./$expName -r;
cp $expName/fort.2 $resultsDir/$expName.1dy;
EOF

  return (CommandEntry -> new({COMMAND => $commandString}))
}

sub checkResults {
  my $env = shift;
  my $tempDir = shift;
  my $npes = shift;
  my $compiler = shift;

  my $installDir    = $env -> {INSTALL_DIR};
  my $rundeck       = $env -> {RUNDECK};
  my $configuration = $env -> {CONFIGURATION};
  my $resultsDir    = $env -> {RESULTS_DIR};
  $resultsDir .="/$compiler";
  my $cmp           = "$resultsDir/CMPE002P";
  my $expName;

  if (@_) {$expName = shift}
  else {$expName = "$rundeck.$configuration.$compiler"};

  my $suffix;
  if ($configuration eq "SERIAL" or $configuration eq "SERIALMP") {
      $suffix = "";
  }
  else {
      $suffix = ".np=$npes";
  }

  my $commandString = <<EOF;
cd $resultsDir;
$cmp $rundeck.SERIAL.1hr $expName.1hr$suffix;
$cmp $rundeck.SERIAL.1dy $expName.1dy$suffix;
EOF

  return (CommandEntry -> new({COMMAND => $commandString, STDOUT_LOG_FILE => "$expName.log"}))
}

sub writeModelErcFile {
    my $env = shift;
    my $modelerc = $ENV{MODELERC};
    my $commandString = "rm $modelerc\n";

    while (my ($var, $value) = each(%$env) ) {
	$commandString .= "echo $var=$value >> $modelerc\n";
    }
    $commandString .= "mkdir -p $env{DECKS_REPOSITORY} $env{CMRUNDIR} $env{SAVEDISK} $env{EXECDIR} \n";
    return (CommandEntry -> new({COMMAND => $commandString}))
}

sub cvsCheckout {
#    my $what = shift;
  my $env = shift;

  my $scratchDirectory = $env -> {SCRATCH_DIRECTORY};
  my $cvsroot          = $env -> {CVSROOT};

  my $commandString = <<EOF;
(pushd $scratchDirectory
cvs -d $cvsroot co modelE
popd)
EOF

  return (CommandEntry -> new({COMMAND => $commandString}))
}

1;
