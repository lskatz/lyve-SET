#! /usr/bin/env perl

=pod
=head1

=over

=item

Schedule::SGELK
A module for submitting jobs to an SGE queue. Monitoring is 
performed using a combination of monitoring files 
written to the hard drive and qstat.
Submitting is performed internally by making a perl script.
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

=over

=item

sub new: create a new instance of a scheduler.
Arguments and their defaults:
  numnodes=>1 maximum nodes to use
  numcpus=>1 maximum cpus that will be used per node
  verbose=>1
  workingdir=>$ENV{PWD} a directory that all nodes can access to read/write jobs and log files
  waitForEachJobToStart=>0 Allow each job to start as it's run (0), or to wait until the qstat sees the job before continuing (1)
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

# submit to the cluster
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
  
  my %return=(jobid=>$jobid,submitted=>$submitted,running=>$running,finished=>$finished,tempdir=>$tempdir,output=>$output,cmd=>$cmd,script=>$script);
  push(@jobsToClean,\%return) if(!$self->settings("keep"));
  push(@jobsToMonitor,\%return);
  return %return if wantarray;
  return \%return;
}

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

# return 1 for finished. 0 for still running.
# if error, return -1 and set $self->error(jobid)
sub checkJob{
  my($self,$job)=@_;
  # see if the job has even started: {submitted}
  return 0 if(!-e $$job{submitted});
  # if the job finished, then great! {finished}
  return 1 if(-e $$job{finished});
  # see if the job is running: {running}
  if(-e $$job{running}){
    # see what the job status is {jobid}
    my $status=$self->jobStatus($$job{jobid});
    if($status eq 'Eqw'){
      $self->error("Command resulted in an error. qstat -j $$job{jobid} for more info\n  $$job{cmd}");
      return -1;
    } elsif($status =~/[rt]/){
      return 0;
    } else { warn "WARNING: I don't know how to interpret status $status for jobid $$job{jobid}. I'm not going to do anything about it."; }
    my @output=read_file($$job{output});
    if($output[-3] =~/QSUB ERROR/){
      $self->error("ERROR: $output[-2], exit code: $output[-1]");
    }
    return -1
  }
  die "Could not understand what the status is of job!\n".Dumper($job);
}
# return the job status of the id
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

# Wait on all jobs to finish and clear out the queue.
sub wrapItUp{
  my($self)=@_;
  $self->waitOnJobs(\@jobsToMonitor,1);
  return 1;
}
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

# remove all the files and everything associated with a job
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
