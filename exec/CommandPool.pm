#!/usr/bin/perl
package CommandPool;
use CommandEntry;

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my $self  = {};
  
  $self -> {TABLE} = [];
  $self -> {SEMAPHORE_PREFIX} = ".command";
  $self -> {FINAL_COMMAND} = 0;
  
  bless ($self, $class);
  makePristine($self -> {SEMAPHORE_PREFIX});
  return $self;
}

sub makePristine{
    my $semaphore = shift;

    system("rm -f $semaphore*");
    return 0;
}

sub numCommands {
  my $self = shift;
  
  my $table = $self -> {TABLE};
  my $n = @$table;
  return $n;
}

sub get {
  my $self = shift;
  my $id = shift;

  return ($self -> {TABLE} -> [$id]);
}

sub add {
  my $self = shift;
  my $newCommand = shift;

  my $n = ($self -> numCommands());
  $self -> {TABLE} -> [$n] = $newCommand;

  if (@_) {$newCommand -> setDependency(shift);}

  return $n;
}

sub setFinal {
  my $self = shift;
  my $newCommand = shift;

  $self -> {FINAL_COMMAND} = $newCommand;

}

sub allComplete {
  my $self       = shift;

  my $id;  
  for ($id = 0; $id < ($self -> numCommands()); $id ++) {
    if (not ($self -> get($id) -> isComplete())) { return 0 };
  }
  return 1;
}

sub semaphore {
  my $self     = shift;
  my $id       = shift;

  return ($self -> {SEMAPHORE_PREFIX} . ".$id");
}

sub resetAll {
  my $self       = shift;
  
  for ($id = 0; $id < ($self -> numCommands()); $id ++) {
    $self -> get($id) -> setStatus(PENDING);
    unlink $self -> semaphore($id);
  }
  return 1;
}

my $MAX_CYCLES=1000;
my $SLEEP=10;

sub run {
  my $self       = shift;
  my $LOG        = shift;
  my $DEBUG      = shift;

  my $counter = 0;
  my $somethingChanged = 0;

  for (;;) {
    ++$counter;
    
    print $LOG "new cycle ... $counter\n" if $DEBUG;
    if ($counter > $MAX_CYCLES) {print $LOG "MAX_CYCLES exceeded - something is very wrong.\n"; return 0};

    if ($self -> allComplete()) {last;};
  
    for (my $id = 0; $id < ($self -> numCommands()); $id ++) {
      my $cmd = $self -> get($id);
      my $semaphore = $self -> semaphore($id);

      if ($cmd -> isReady()) {
        print $LOG " ... Launching $id \n" if $DEBUG;
	print $LOG "     ...   " . $cmd -> {COMMAND} . "\n"  if $DEBUG;
        $cmd -> setStatus(RUNNING);
	
        $cmd -> launch($semaphore);
	$somethingChanged = 1;
      }
      elsif ($cmd -> isRunning()) {
        if (-e $semaphore) {
          print $LOG " ... Command $id has completed.\n" if $DEBUG;
          unlink $semaphore;
          $cmd -> setStatus(COMPLETE);
	$somethingChanged = 1;
        }
      }
      elsif ($cmd -> isComplete()) {
      }
      else {
        print $LOG " ... $id is waiting on dependencies.\n " if $DEBUG;
      }
    }

    if ($somethingChanged) {$somethingChanged =0}
    else { print $LOG " ... sleeping for a while ... \n" if $DEBUG;
	   sleep $SLEEP; }

  }
# Do any final task
  if ($self -> {FINAL_COMMAND}) {
        print $LOG "Final task ...\n";
	$self -> {FINAL_COMMAND} -> launch(semaphore("/dev/nul"))
	};

  return 1;
}

1;
