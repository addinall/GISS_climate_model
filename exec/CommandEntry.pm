#!/usr/bin/perl
package CommandEntry;

sub new {
  my $proto = shift;
  my ($args) = @_;

  my $class = ref($proto) || $proto;
  my $self  = {};
  
  $self -> {COMMAND} = " ";
  $self -> {NUM_PROCS} = 1;
  $self -> {QUEUE} = LOCAL; # or BATCH
  $self -> {DEPENDENCY} = ""; # place holder for list of references to other commands
  $self -> {STATUS} = PENDING;
  $self -> {STDOUT_LOG_FILE} = "LOG";
  $self -> {STDERR_LOG_FILE} = "ERR";

  # override defaults with arguments
  while ( my ($key, $value) = each(%$args) ) {
    $self -> {$key} = $value;
  }

  bless ($self, $class);
  return $self;
}

sub dependency {
  my $self       = shift;

  return ($self -> {DEPENDENCY});
}

sub status {
  my $self = shift;
  return ($self -> {STATUS});
}

sub isComplete {
  my $self       = shift;

  return ( $self -> {STATUS} eq COMPLETE );
}

sub isRunning {
  my $self       = shift;

  return ( $self -> {STATUS} eq RUNNING );
}

sub isPending {
  my $self       = shift;

  return ( $self -> {STATUS} eq PENDING );
}

sub isReady {
  my $self       = shift;
  
  my $pending = ($self -> isPending());
  my $dependency = ($self -> dependency());

  if ($dependency eq "") { return $pending }
  else {return ($pending and ($self -> dependency() -> isComplete()) )};
}

sub setStatus {
  my $self = shift;
  my $newStatus = shift;
  $self -> {STATUS} = $newStatus;
}

sub setDependency {
  my $self = shift;
  my $newDependency = shift;
  $self -> {DEPENDENCY} = $newDependency;
}

sub runInBatch {
  my $self = shift;
  my $commandString = shift;
  my $queue = shift;
  my $queueString = " ";
  if ($queue) {$queueString = "-q $queue\n";};

  my $NODE_SIZE = 8;

  my $nCPUS = $self -> {NUM_PROCS};
  my $nodes = 1 + int(($nCPUS-1)/$NODE_SIZE);
  my $ncpus = $nCPUS;
  $ncpus = $NODE_SIZE if ($ncpus > $NODE_SIZE);

  $walltime = "3:00:00\n";
#  if ($queue == "datamove") {$walltime = "0:30:00\n"};

  my $script = <<EOF;
#!/bin/bash
#PBS -l select=$nodes:ncpus=$ncpus:proc=neha
#PBS -l walltime=$walltime
#PBS -W group_list=a940a
#PBS -N regression
#PBS -j oe
##PBS -o foo
#PBS $queueString

cd \$PBS_O_WORKDIR
. /usr/share/modules/init/bash
module purge
module load comp/intel-10.1.023 lib/mkl-10.1.2.024 mpi/impi-3.2.1.009

$commandString
EOF

`echo '$script' > pbstmp; qsub -V pbstmp; rm pbstmp`;

}

sub launch {
  my $self = shift;
  my $semaphore = shift;
  my $logFile = $self -> {STDOUT_LOG_FILE};
  my $logErr = $self -> {STDERR_LOG_FILE};

  my $commandString   = "(" . $self -> {COMMAND} . ") 2>&1 >> $logFile; /usr/bin/touch $semaphore;";

  my $mode = $self -> {QUEUE};

  if ($self -> {QUEUE} eq LOCAL) {
    `($commandString) &`; # run in background
  }
  else {
    $self -> runInBatch("$commandString", $self -> {QUEUE});
  }
}

sub launchMultiple {
  my $cmds       = shift;
  my $semaphores = shift;

  my $commandString;

  foreach (@$cmds) {
    $commandString .= "(";
    $commandString .= ($_) -> {COMMAND};
    $commandString .= "; touch " . shift(@$semaphores) . ") &;";
  }
  $commandString .= "wait";

  if ($self -> {QUEUE} eq LOCAL) {
    `$commandString`;
  }
  else {
    $self -> runInBatch($commandString, $self -> {QUEUE});
  }
}

1;
