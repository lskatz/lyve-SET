#!/usr/bin/env perl
# Author:Lee Katz <lkatz@cdc.gov>
# Thanks: Darlene Wagner for giving me this idea

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use List::MoreUtils qw(uniq);
use FindBin;
use Bio::Perl;

$ENV{PATH}="$FindBin::RealBin/../lib/phispy/phiSpyNov11_v2.3:".$ENV{PATH};

local $0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}

exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help species=s numcpus=i));
  $$settings{numcpus}||=1;

  my $fasta=$ARGV[0];
  die usage() if(!$fasta || $$settings{help});

  my $phispyCoordinates=phispy($fasta,$settings);
  my $phastCoordinates=phast($fasta,$settings);

  my %seq=readFasta($fasta,$settings);

  for ($phispyCoordinates,$phastCoordinates){
    maskSequences(\%seq,$_,$settings);
  }

  # print the masked sequences
  while(my ($id,$seq)=each(%seq)){
    $seq=~s/(.{60})/$1\n/g;
    chomp($seq);
    print "$id\n$seq\n";
  }

  return 0;
}

sub phispy{
  return [];
}

sub phast{
  my($fasta,$settings)=@_;

  #blastx -query 2014L-6068.fasta -db ~/bin/lyve-SET/lib/phast/phast.faa -evalue 0.05 -outfmt 6 -num_threads 8
  my $db="$FindBin::RealBin/../lib/phast/phast.faa";
  my $allResults=`blastx -query '$fasta' -db $db -evalue 0.05 -outfmt 6 -num_threads $$settings{numcpus}`;
  die "ERROR with blastx: $!" if $?;

  my(@range);
  for my $result(split(/\n/,$allResults)){
    $result=~s/^\s+|\s+$//g; # trim
    my ($contig,$hit,$identity,$length,$gaps,$mismatches,$sstart,$send,$qstart,$qend,$e,$score)=split /\t/, $result;
    next if($score < 50);
    
    push(@range,"$contig:$qstart-$qend");
  }
  @range=uniq(@range);
  @range=sort{
    my($contigA,$rangeA)=split(/:/,$a);
    my($startA,$stopA)=split(/\-/,$rangeA);
    my($contigB,$rangeB)=split(/:/,$b);
    my($startB,$stopB)=split(/\-/,$rangeB);

    return $contigA cmp $contigB if($contigA ne $contigB);
    $startA<=>$startB; 
  } @range;

  return \@range;
}

sub maskSequences{
  my($seq,$ranges,$settings)=@_;
  
  for my $r(@$ranges){
    my($contig,$range)=split(/:/,$r);
    my($start,$stop)=split(/\-/,$range);
    die "ERROR: could not find contig $contig in the input fasta sequence!" if(!$$seq{$contig});
    #($contig,$start,$stop)=("NODE38length2280cov11.704", 5, 10); # debug
    my $length=abs($stop-$start)+1;

    # convert to zero-base
    $start--;
    $stop--;

    # use substr to switch out phages with N
    $$seq{$contig}=substr($$seq{$contig},0,$start) . ("N" x $length) . substr($$seq{$contig},$stop+1);
  }
}

sub readFasta{
  my($fasta,$settings)=@_;
  my $in=Bio::SeqIO->new(-file=>$fasta);
  my %seq;
  while(my $seq=$in->next_seq){
    $seq{$seq->id}=$seq->seq;
  }
  return %seq;
}

sub usage{
  "Finds phages in a fasta file using phast and PhiSpy
  Usage: $0 file.fasta
  --numcpus 1
  "
}
