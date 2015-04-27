package LyveSET;
use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw/fileparse basename dirname/;
use Term::ANSIColor;
use Data::Dumper;

our @EXPORT_OK = qw(logmsg @fastqExt @fastaExt @bamExt);

local $0=basename $0;

######
# CONSTANTS

our @fastqExt=qw(.fastq.gz .fastq .fq .fq.gz);
our @fastaExt=qw(.fasta .fna .faa .mfa .fa);
our @bamExt=qw(.sorted.bam .bam);

###########################################
### COMMON SUBS (not object subroutines) ##
###########################################
# Redefine how the Lyve-SET scripts die
$SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
# Centralized logmsg
sub logmsg {print STDERR "$0: ".(caller(1))[3].": @_\n";}

###########################################
### CLASS DEFINITION ######################
###########################################

# Settings is a hash
#  min, max => range of values expected, inclusive
#  fp       => where LyveSET prints out to. Default: STDERR
# Invocation: $set=LyveSET->new(option1=>...)
sub new{
  my($class,%settings)=@_;
  my $self={settings=>\%settings};
  bless($self,$class);

  $self->progressCharacterSet($settings{min},$settings{max},\%settings);

  return $self;
}

sub progressCharacterSet{
  my($self,$min,$max,$settings)=@_;
  if($max < $min){
    die "ERROR: the max is larger than the min";
  }

  # Greyscale colors; ascii characters
  my @char=qw(_ . o 0 O);
  my @color=reverse(10..23); $_="grey$_" for(@color);

  # Make an array of characters to choose from
  my @characterSet;
  for(my $i=0;$i<@char;$i++){
    for(my $j=0;$j<@color;$j++){
      push(@characterSet,[$char[$i],$color[$j]]);
    }
  }

  # Save into the object
  $$self{progressbar}{character}=\@characterSet;
  $$self{progressbar}{numcategories}=scalar(@char) * scalar(@color);
  $$self{progressbar}{numPerCategory}=int(($max-$min+1)/$$self{progressbar}{numcategories});
}

# Give an integer, print out a progress character
sub printProgress{
  my($self,$int)=@_;
  my $settings=$$self{settings};
  # Where to print to
  my $fp = $$settings{fp}||*STDERR;
  # What index of the character set to use
  my $percentile=($int-$$settings{min})/($$settings{max}-$$settings{min});
  my $index=int($$self{progressbar}{numcategories}*$percentile)-1;
  $index=0 if($index < 0);

  # Print something pretty
  my $charStruct=$$self{progressbar}{character}[$index];
  print $fp color $$charStruct[1]; # change colors
  print $fp $$charStruct[0];       # progress
  print $fp color "reset";         # change back
}

1; # gotta love how we we return 1 in modules. TRUTH!!!

