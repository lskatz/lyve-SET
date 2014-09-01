#!/usr/bin/env perl
# Find the distribution of distances between SNPs in a set of VCF files

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Statistics::Descriptive;
use Math::Round qw/nearest/;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help distribution=s binning=i outputType=s)) or die $!;
  $$settings{binning}||=100;
  $$settings{outputType}||="default";
  $$settings{distribution}||="nonparametric";
  $$settings{distribution}=lc($$settings{distribution});
  die usage() if(!@ARGV || $$settings{help});

  # find distances for each VCF
  my @distance;
  for my $vcf(@ARGV){
    my $vcfHash=readVcf($vcf);
    # get the distance distribution for the genome
    while(my($seqname,$posHash)=each(%$vcfHash)){
      # distances for the contig
      my @pos=sort {$a <=> $b} keys(%$posHash);
      for(my $i=0;$i<@pos-1;$i++){
        push(@distance,$pos[$i+1]-$pos[$i]);
      }
    }
  }
  # sort the distances - not necessary but good for pre-emptive debugging
  @distance=sort {$a <=> $b} @distance;
  # Bin the distances according to the binning size.
  my @binnedDistance=binData(\@distance,$settings);

  # alternative output types
  if($$settings{outputType} eq 'distances'){
    print $_."\n" for(@distance);
    return 0;
  }elsif($$settings{outputType} eq 'binned-distances'){
    print $_."\n" for (@binnedDistance);
    return 0;
  }

  die "I do not understand the output type you requested: $$settings{outputType}" if($$settings{outputType} ne 'default');

  my $cutoff;
  if($$settings{distribution} eq 'bimodal'){
    $cutoff=bimodalCutoff(\@binnedDistance,\@distance,$settings);
  } else {
    $cutoff=nonparametricCutoff(\@binnedDistance,\@distance,$settings);
  }
  $cutoff=nearest(1,$cutoff);
  print "$cutoff\n";
  return 0;
}

sub nonparametricCutoff{
  my($binnedDistance,$distance,$settings)=@_;
  my $stat=Statistics::Descriptive::Full->new;
  $stat->add_data(@$distance);
  my $fivePercent=$stat->percentile(5);
  return $fivePercent;
}

sub bimodalCutoff{
  my($binnedDistance,$distance,$settings)=@_;

  # Find the two modes of the binned distances
  my($larger,$smaller)=findModes($binnedDistance,$settings);
  # Split the data into two distributions
  my $midpoint=($smaller+$larger)/2;
  my(@smaller,@larger);
  for(@$distance){
    if($_ >= $midpoint){
      push(@larger,$_);
    } else {
      push(@smaller,$_);
    }
  }

  # Report the smaller distribution's average +2 stdev
  my($avg,$stdev)=avgStdev(@smaller);
  my $cutoff=$avg+2*$stdev;
  return $cutoff;
}

sub binData{
  my($data,$settings)=@_;
  my @binned;
  for(my $i=0;$i<@$data;$i++){
    push(@binned,nearest($$settings{binning},$$data[$i]));
  }
  return @binned;
}

sub findModes{
  my($data,$settings)=@_;
  # Count how many times each data point occurs.
  my %data;
  for(@$data){
    $data{$_}++;
  }
  
  # Find the modes with the largest and then second largest count.
  my @sortedCount=sort {$a<=>$b} values(%data);
  my($largest,$secondLargest)=(0,0);
  while(my($dataPoint,$count)=each(%data)){
    $largest=$dataPoint if($count==$sortedCount[-1]);
    $secondLargest=$dataPoint if($count==$sortedCount[-2]);
  }
  return ($largest,$secondLargest);
}

sub avgStdev{
  my @data=@_;
  my $stat=Statistics::Descriptive::Full->new;
  $stat->add_data(@data);
  my($mean,$stdev)=(
    sprintf("%0.2f",$stat->mean),
    sprintf("%0.2f",$stat->standard_deviation)
  );
  return($mean,$stdev);
}

sub readVcf{
  my($vcf,$settings)=@_;
  my $vcfHash={};
  #$diskIoStick->down; # mark that one process is using the disk
  open(VCF,"<",$vcf) or die "ERROR: could not open vcf file $vcf:$!";
  while(<VCF>){
    next if(/^#/);
    chomp;
    my($contig,$pos,undef,$ref,$alt)=split /\t/;
    $$vcfHash{$contig}{$pos}=$alt;

    # Indels will just be an N because it is too difficult to deal with those.
    if($ref eq '*' || $alt eq '*' || length($ref)>1 || length($alt)>1){
      # lowercase N to mask it later, when trying to recover other bases for
      # low coverage. This N is not due to low coverage.
      $$vcfHash{$contig}{$pos}='n';
    }
  }
  close VCF;
  #$diskIoStick->up; # mark that one process is no longer using the disk
  return $vcfHash;
}

sub usage{
  "Given a set of VCFs and their coordinates, determines how close a SNP must be in order for it to be 'clustered.'  The 'clustered' distance is at the 95% cutoff.
  Usage: $0 *.vcf
  -b binning (default: 10)
  -o default
    choices: default, distances, binned-distances
  -d distribution
    choices: nonparametric (default), bimodal
  "
}
