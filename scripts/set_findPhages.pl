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

local $0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}

exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help species=s numcpus=i));
  $$settings{numcpus}||=1;

  my $fasta=$ARGV[0];
  die usage() if(!$fasta || $$settings{help});

  my $phastCoordinates=phast($fasta,$settings);

  my %seq=readFasta($fasta,$settings);

  my $ranges=findRangeUnion(\%seq,$phastCoordinates,$settings);

  # Print the ranges to stdout
  for my $r(@$ranges){
    print join("\t",@$r)."\n";
  }
  return 0;

  # print the masked sequences
  while(my ($id,$seq)=each(%seq)){
    $seq=~s/(.{60})/$1\n/g;
    chomp($seq);
    print ">$id\n$seq\n";
  }

  return 0;
}

sub phast{
  my($fasta,$settings)=@_;

  my $db="$FindBin::RealBin/../lib/phast/phast.faa";
  logmsg "Running blastx against $db";
  my $allResults=`blastx -query '$fasta' -db $db -evalue 0.05 -outfmt 6 -num_threads $$settings{numcpus}`;
  #die $allResults;
  die "ERROR with blastx: $!" if $?;

  logmsg "Parsing results";
  my(@range);
  for my $result(split(/\n/,$allResults)){
    $result=~s/^\s+|\s+$//g; # trim
    my ($contig,$hit,$identity,$length,$gaps,$mismatches,$sstart,$send,$qstart,$qend,$e,$score)=split /\t/, $result;
    next if($score < 50 || $length < 20);
    
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

sub findRangeUnion{
  my($seq,$ranges,$settings)=@_;

  my $points=simplifyRangesToPoints($ranges,$settings);
  return pointsToRanges($points,$settings);
  
}

sub simplifyRangesToPoints{
  my($ranges,$settings)=@_;
  
  # 1) Convert to a better format
  # 2) Find the union as we go
  my (%unionIndex);
  my %range;
  for my $r(@$ranges){
    my($contig,$range)=split(/:/,$r);
    my($start,$stop)=split(/\-/,$range);
    $start--; $stop--;
    
    # Convert the range to a true perl range, with a contig prefix
    my @range=($start..$stop);
    $_=$contig.":".$_ for(@range);

    # Index the range
    @unionIndex{@range} = (1) x @range;
  }
  
  # sort the points by contig/pos
  my @union=sort {
    my($contigA,$posA)=split(/:/,$a);
    my($contigB,$posB)=split(/:/,$b);
    return $contigA cmp $contigB if($contigA ne $contigB);
    return $posA<=>$posB;
  } keys(%unionIndex);
  return \@union;
}

# Points should already be sorted
sub pointsToRanges{
  my($points,$settings)=@_;
  my @range;
  my ($currentContig,$start)=split(/:/,$$points[0]);
  my $end=$start;
  my $numPoints=@$points;
  for(my $i=1;$i<$numPoints-1;$i++){
    my($contig,$pos)=split(/:/,$$points[$i]);
    my($lastContig,$lastPos)=split(/:/,$$points[$i-1]);
    my($nextContig,$nextPos)=split(/:/,$$points[$i+1]);

    # The range ends if the next integer is not subsequent.
    if($contig ne $lastContig || $pos-$lastPos!=1){
      $end=$lastPos; # cut off the range

      push(@range,[$currentContig,$start,$end]);

      $currentContig=$contig;
      $start=$pos;
    }

    # Record one final range at the end
    if($i==$numPoints-2 && $contig eq $nextContig){
      push(@range,[$nextContig,$start,$nextPos]);
    }
  }

  return \@range;
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
