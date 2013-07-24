#! /usr/bin/env perl

=pod

=head1 NAME

Schedule::SGELK

=head1 SYNOPSIS

A module for submitting jobs to an SGE queue.
  use Schedule::SGELK
  my $sge=Schedule::SGELK->new(-verbose=>1,-numnodes=>5,-numcpus=>8,-workingdir=>"SGE/",-waitForEachJobToStart=>1);
  $sge->set("jobname","thisisaname");
  # run a series of jobs and wait for them all to finish
  my $job=$sge->pleaseExecute("sleep 60");
  my $job2=$sge->pleaseExecute("sleep 60");
  $sge->wrapItUp();
  # or you can specify which jobs to wait for
  $sge->waitOnJobs([$job,$job2],1); # 1 means wait for all jobs to finish; 0 to wait for any free node
  # or in one step
  $sge->pleaseExecute_andWait("sleep 60");

=head1 DESCRIPTION

A module for submitting jobs to an SGE queue. Monitoring is 
performed using a combination of monitoring files 
written to the hard drive and qstat.
Submitting is performed internally by making a perl script.

=head1 AUTHOR

Author: Lee Katz <lkatz@cdc.gov>

=cut

package Schedule::SGELK;
use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/basename/;
use File::Spec;
use File::Slurp qw/read_file write_file/;
use File::Temp qw/tempdir/;
use String::Escape qw/escape/;

sub logmsg {local $0=basename $0;my $FH = *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}
local $SIG{'__DIE__'} = sub { my $e = $_[0]; cleanAllJobs(); $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

# to be called when the script exits
my @jobsToClean=();
my @jobsToMonitor=();
sub cleanAllJobs{
  logmsg "TODO: Figure out why this sub isn't being called when the script dies";
  logmsg "Cleaning all jobs";
  for (@jobsToClean){
    cleanAJob($_);
  }
}

=pod

=head2 METHODS

=over

=item sub new

create a new instance of a scheduler.
Arguments and their defaults:
  numnodes=>1 maximum nodes to use
  numcpus=>1 maximum cpus that will be used per node
  verbose=>1
  workingdir=>$ENV{PWD} a directory that all nodes can access to read/write jobs and log files
  waitForEachJobToStart=>0 Allow each job to start as it's run (0), or to wait until the qstat sees the job before continuing (1)
  jobname=>... This is the name given to the job when you view it with qstat. By default, it will be named after the script that calls this module.

=back

=cut

# args can be numnodes=>6,numthreads=>8,verbose=>0
sub new{
  my($class,%args)=@_;
  my $self=bless{},$class;

  # just start with verbosity = 0. This avoids an error in the >= checks
  $self->{'settings'}={};
  $self->{'error'}="";
  $self->{'exit_code'}=0;

  # load up if we know what we have
  foreach my $key (keys %args) {
    my $nodash=$key;
    $nodash =~ s/^\-//;
    if ($nodash eq "executable") {
      $self->executable($args{$key});
    }
    elsif ($nodash eq "use_cwd") {
      $self->use_cwd($args{$key});
    }
    elsif ($nodash eq "name") {$self->name($args{$key})}
    else {
      #$self->{$nodash}=$args{$key};
      $self->set($nodash,$args{$key});
    }
  }

  # set defaults if they are not set
  my %default=(numnodes=>1,numcpus=>1,verbose=>1,workingdir=>$ENV{PWD},waitForEachJobToStart=>0);
  while(my($key,$value)=each(%default)){
    $self->settings($key,$value) if(!defined($self->settings($key)));
  }

  return $self;
}

=pod

=over

=item sub error($msg,$exit_code) or error()

Get or set the error. Can set the error code too, if provided.

=back

=cut

# Sets the error and returns the previous error message.
# Or, simply returns the current error message.
sub error{
  my($self,$msg,$exit_code)=@_;
  return $self->{error} if(!defined($msg));
  my $oldmsg=$self->{error};
  $self->{error}=$msg;
  $self->{exit_code}=$exit_code if(defined($exit_code));
  return $oldmsg;
}
=pod

=over

=item sub set() get() settings()

Get or set a setting. All settings are listed under sub new().
If a setting is provided without a value, then nothing new will be set,
and only the value of the specified setting will be returned.

=back

=cut

sub set{
  my($self,$key,$value)=@_;
  if(defined($key) && defined($value)){
    $self->{'settings'}{$key}=$value;
  } elsif (defined($key)){
    return $self->{'settings'}{$key};
  }
  return %{$self->{'settings'}} if wantarray;
  return $self->{'settings'};
}
# renaming of sub set()
sub get{
  my($self,@args)=@_;
  return $self->set(@args);
}
# renaming of sub set()
sub settings{
  my($self,@args)=@_;
  return $self->set(@args);
}

=pod

=over

=item pleaseExecute()

This is the main method. It will submit a command to the cluster.
  $sge->set("jobname","a_nu_start");
  $sge->pleaseExecute("someCommand with parameters");

If you are already occupying more than numnodes, then it will pause before 
executing the command. It will also create many files under workingdir, so
be sure to specify it. The workingdir is the current directory by default.

=back

=cut

sub pleaseExecute{
  my($self,$cmd)=@_;
  local $0=basename $0;
  my $settings=$self->settings;
  my $jobid=-1; # -1 is an error state
  $$settings{jobname}||="run$0";
  $$settings{logfile}||="$0.log";
  $$settings{numcpus}||=1;
  $$settings{timeout}||=60; # how long we will wait for qsub to start
  return 0 if($cmd=~/^\s*$/); # if there's no command, then no worries

  $self->waitOnJobs(\@jobsToMonitor,0); # wait until the job can be submitted

  my $rand=int(rand(999999));
  my $tempdir=$$settings{workingdir};
  # create a bash script with the literal command in it
  my $script="$tempdir/qsub.$rand.pl";

  my $output="$$settings{workingdir}/$0.$rand.log";
  open(SCRIPT,">",$script) or die "Could not write to temporary script: $!";
  print SCRIPT "#!/usr/bin/env perl\n\n";
  #   it has SGE params in it
  print SCRIPT "#\$ -N $$settings{jobname}\n";
  print SCRIPT "#\$ -V\n";
  print SCRIPT "#\$ -wd $ENV{PWD}\n";
  print SCRIPT "#\$ -pe smp $$settings{numcpus}\n";
  print SCRIPT "#\$ -o $output\n";
  print SCRIPT "#\$ -e $output\n";
  print SCRIPT "use strict;\nuse warnings;\n";
  print SCRIPT "use File::Slurp qw/read_file write_file/;\n";

  my($submitted,$running,$finished)=("$$settings{workingdir}/$0.$rand.submitted", "$$settings{workingdir}/$0.$rand.running", "$$settings{workingdir}/$0.$rand.finished");
  # announces that it was submitted
  my $sanitized=escape('qqbackslash',$cmd);
  print SCRIPT "write_file('$submitted',$sanitized\n);\n";
  # it runs the command
  print SCRIPT "write_file('$running',$sanitized\n);\n";
  print SCRIPT "system($sanitized);\n";
  # print a parsable error if the script dies. This error will be in $output
  print SCRIPT "die \"QSUB ERROR\n\$?\n\$!\" if(\$?);\n";
  # announces when it is finished
  print SCRIPT "write_file('$finished',$sanitized\n);\n";
  close SCRIPT;
  system("touch $script"); die if $?; # make the system close the script. Why isn't Perl closing it?
  #system("cat $script");sleep 60;die;

  # now run the script and get the jobid
  my $out=`qsub $script`; chomp($out);
  #logmsg $out."\n  Waiting on a confirmation";
  if($out=~/Your job (\d+)/){
    $jobid=$1;
  } else {
    logmsg "WARNING: the last job submitted did not have an obvious jobid. It can't be tracked!";
  }

  # monitor for the script to be running before moving on
  my $started=time;
  while(!-e $submitted){
    last if(!$self->settings("waitForEachJobToStart"));
    sleep 1;
    die "Command timed out!\n  $cmd" if((time-$started)>$$settings{timeout});
    die "Command resulted in an error. qstat -j $jobid for more info\n  $cmd" if($self->jobStatus($jobid) eq 'Eqw');
  }

  # TODO create a link from the jobid to the random id
  
  my %return=(jobid=>$jobid,submitted=>$submitted,running=>$running,finished=>$finished,tempdir=>$tempdir,output=>$output,cmd=>$cmd,script=>$script);
  push(@jobsToClean,\%return) if(!$self->settings("keep"));
  push(@jobsToMonitor,\%return);
  return %return if wantarray;
  return \%return;
}

=pod

=over

=item pleaseExecute_andWait()

Exact same as pleaseExecute(), except it will wait for the command to finish 
before continuing. Internally calls pleaseExecute() and then waitOnJobs().
However one key difference between pleaseExecute() and this sub is that you can
give a list of commands.

  # this will take 100 seconds because all commands have to finish.
  $sge->pleaseExecute_andWait(["sleep 60","sleep 100","sleep 3"]);

=back

=cut

sub pleaseExecute_andWait{
  my($self,$cmd)=@_;
  my %settings=$self->settings;
  my $mustfinish=$self->settings("mustfinish"); # should be restored later
  $self->set("mustfinish",0);
  $cmd=[$cmd] if(ref($cmd) eq ""); # let cmd be a string but turn it into a list internally
  my(@jobid);
  for(@$cmd){
    my $jobid=$self->pleaseExecute($_);
    push(@jobid,$jobid);
    $self->waitOnJobs(\@jobid);
  }
  $self->waitOnJobs(\@jobid,1);
}

=pod

=over

=item checkJob($jobHash)

Checks the status of a given job. The job variable is obtained from pleaseExecute().
$self->error can be set if there is an error in the job. Return values:
1 for finished; 0 for still running or hasn't started; -1 for error.

=back

=cut 

sub checkJob{
  my($self,$job)=@_;
  # see what the job status is {jobid}
  my $status=$self->jobStatus($$job{jobid});
  return 0 if($status eq 'qw'); # queued but not running
  # see if the job has even started: {submitted}
  return 0 if(!-e $$job{submitted});
  # if the job finished, then great! {finished}
  return 1 if(-e $$job{finished});
  # see if the job is running: {running}
  if(-e $$job{running}){
    if($status eq 'Eqw'){
      $self->error("Command resulted in an error. qstat -j $$job{jobid} for more info\n  $$job{cmd}");
      return -1;
    } elsif($status =~/[rt]/){
      return 0;
    } elsif($status=~/0/){
      logmsg "Warning: job $$job{jobid} indicates it has started running but hasn't finished. However, qstat indicates it is no longer running. I will check intermittently to see what the status is." if((caller(1))[3] !~/checkJob/);
      for my $i(0..50){
        sleep (($i+1)*2); # the sleep timer's interval goes slower each time
        my $newstatus=$self->checkJob($job);
        return $newstatus if($newstatus);
      }
    } else { 
      logmsg "WARNING: I don't know how to interpret status $status for jobid $$job{jobid}. I'm not going to do anything about it."; 
    }
    my @output=read_file($$job{output});
    if($output[-3] =~/QSUB ERROR/){
      $self->error("ERROR: $output[-2], exit code: $output[-1]");
    }
    return -1;
  }
  logmsg "ERROR: Could not understand what the status is of job $$job{jobid}!\n".Dumper($job);
}

=pod

=over

=item jobStatus(jobid)

Given an SGE job id, it returns its qstat status

=back

=cut

#TODO cache the qstat results so that it is only run every 1-2 seconds. 
#  I don't want to overload the machine with system calls.
sub jobStatus{
  my($self,$jobid)=@_;
  my $state=0;
  open(QSTAT,"qstat|") or die "ERROR: could not execute qstat! $!";
  while(<QSTAT>){
    s/^\s+|\s+$//g;
    my @F=split /\s+/;
    if($F[0] eq $jobid){
      $state=$F[4];
    }
  }
  close QSTAT;
  return $state;
}

=pod

=over

=item wrapItUp()

Waits on all jobs to finish before pausing the program.
Calls waitOnJobs internally. Does not take any parameters.

=back

=cut

# Wait on all jobs to finish and clear out the queue.
sub wrapItUp{
  my($self)=@_;
  $self->waitOnJobs(\@jobsToMonitor,1);
  return 1;
}

=pod

=over

=item waitOnJobs($jobList,[$mustFinish])

Waits on all given jobs to finish. The job list are jobs as given by pleaseExecute().
If $mustFinish evaluates to true, then the program will pause until
all jobs are finished.
Calls on checkJob() internally. Will die with an error message if a job dies.

=back

=cut

# Wait on enough jobs to finish before returning.
# If a job finishes, splice it from the job array.
sub waitOnJobs{
  my($self,$job,$mustfinish)=@_;
  my %settings=$self->settings;
  $settings{mustfinish}=$mustfinish if(defined($mustfinish));
  logmsg "We have reached node capacity ($settings{numnodes})! Waiting for a job to finish." if(@$job >= $settings{numnodes} && $settings{verbose});
  while(@$job > 0){
    for(my $i=0;$i<@$job;$i++){
      my $state=$self->checkJob($$job[$i]);
      if($state==1){
        logmsg "A job finished: $$job[$i]{jobid}" if($settings{verbose});
        splice(@$job,$i,1);
        last;
      } elsif($state==-1){
        die "A job failed! Look at $$job[$i]{output} for more details.\nError message was ".$self->error()."\n".Dumper($$job[$i]);
      }
    }
    sleep 1;
    # break out if you don't have to finish yet but you can still add in another job
    last if(!$settings{mustfinish} && @$job<$settings{numnodes});
  }
  return @$job;
}

=pod

=over

=item cleanAJob

This is internally used for cleaning up files after a job is done. 
Do not use externally.

=back

=cut

sub cleanAJob{
  my($job)=@_;
  for (qw(running submitted finished output script)){
    logmsg $$job{$_};
    unlink $$job{$_};
  }
  system("qdel $$job{jobid}");
  die "Internal error" if $?;
  return 1;
}

sub mktempdir{
  my ($self,$settings) = @_;
  $settings||={};
  # SGELK.22623.XXXXX
  my $tempdir_path = File::Spec->join(File::Spec->tmpdir(), (split("::",(caller(1))[3]))[1].".$$.XXXXX");
  my $tempdir = tempdir($tempdir_path, CLEANUP => !($$settings{keep}));
  return $tempdir;
}


1;
