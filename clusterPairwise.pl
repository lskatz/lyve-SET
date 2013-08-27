#!/usr/bin/env perl

# clusters genomes based on distances.  Uses hierarchical bottom-up clustering.

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

sub logmsg{$|++;print STDERR "@_\n";$|--}
exit(main());
sub main{
  my $settings={};
  GetOptions($settings,qw(help maxdist=s perCluster=i histogram));
  die usage() if($$settings{help});
  $$settings{perCluster}||=2;

  return histogram($settings) if($$settings{histogram});

  $$settings{maxdist}||=die "ERROR: need a maximum distance\n".usage();
  my $neighbor=identifyCloseNeighbors($settings);
  my $cluster=findClusters($neighbor,$settings);
  logmsg "Found ".scalar(@$cluster)." clusters";
  return 0;
}

# parse the input file and retain only close neighbors
sub identifyCloseNeighbors{
  my($settings)=@_;
  my %neighbor;
  my %dist;
  my $maxdist=$$settings{maxdist};
  while(<>){
    chomp;
    my @F=split /\t/;
    next if($F[2] > $maxdist);
    push(@{$neighbor{$F[0]}},$F[1]);
    push(@{$neighbor{$F[1]}},$F[0]);
    $dist{$F[0]}{$F[1]}=$F[2];
    $dist{$F[1]}{$F[0]}=$F[2];
  }
  return (\%neighbor,\%dist) if wantarray;
  return \%neighbor;
}

sub findClusters{
  my($neighbor,$settings)=@_;
  my %neighbor=%$neighbor; # make a copy so that the internal pointer isn't messed up

  my @cluster;
  while(my($centerNode,undef)=each(%neighbor)){
    my $seen={};
    my $cluster=findClusterAroundNode($centerNode,\%neighbor,$seen,$settings);
    next if(!@$cluster);
    unshift(@$cluster,$centerNode);
    next if(@$cluster<$$settings{perCluster});
    print join("\t",scalar(@$cluster),@$cluster)."\n";
    push(@cluster,$cluster);
  }
  return \@cluster;
}

# given a "center" node, find the cluster around it
sub findClusterAroundNode{
  my($centerNode,$neighbor,$seen,$settings)=@_;
  my @cluster; # the output
  
  # Failsafe: don't keep processing if you are recursing too deep
  return \@cluster if($$settings{cluster_level}++ > 9999);

  my $neighbors=$$neighbor{$centerNode};
  for my $n(@$neighbors){
    next if(!$n || $$seen{$n});
    $$seen{$n}=1;
    deleteNode($centerNode,$neighbor,$settings);
    my $secondaries=findClusterAroundNode($n,$neighbor,$seen,$settings);
    push(@cluster,$n,@$secondaries);
  }

  return \@cluster;
}

# delete a node once it has been used, so that we don't backtrack
sub deleteNode{
  my($node,$neighbor,$settings)=@_;
  delete($$neighbor{$node});
  for my $otherNode(keys(%$neighbor)){
    for(my $i=0;$i<@{$$neighbor{$otherNode}};$i++){
      my $otherOtherNode=@{$$neighbor{$otherNode}}[$i] || '';
      delete $$neighbor{$otherNode}[$i] if($otherOtherNode eq $node);
    }
  }
  return $neighbor;
}

sub histogram{
  my($settings)=@_;
  my(%hist);
  while(<>){
    chomp;
    my $dist=(split /\t/)[2];
    my $rounded=100*int($dist/100);
    $hist{$rounded}++;
  }
  my @hist=sort {$a<=>$b} keys(%hist);
  for my $num(@hist){
    print "$num $hist{$num}\n";
  }
  return 0;
}


## histogram
# cut -f 3 pairwiseDist-MT.txt | perl -lane '$rounded=100*int($F[0]/100);$num{$rounded}++;END{while(($num,$count)=each(%num)){print "$count\t$num";}}' | sort -k2,2n 
sub usage{
  "Produces clusters given a three column file from pairwiseDistances.pl
  genome1  genome2  distance
  Usage: $0 -m maxdist < distances.txt
  -m the maximum distance between two nodes in a cluster.
  -p the minimum number per cluster, before reporting it (default: 2)
  --help this help menu
  --histogram to display a histogram and then exit
  "
}
