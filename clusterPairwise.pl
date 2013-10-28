#!/usr/bin/env perl

# clusters genomes based on distances.  Uses hierarchical bottom-up clustering.

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
no warnings 'recursion';
use POSIX;

sub logmsg{$|++;print STDERR "@_\n";$|--}
exit(main());
sub main{
  my $settings={};
  GetOptions($settings,qw(help maxdist=s perCluster=i oerCluster=i histogram));
  die usage() if($$settings{help});
  $$settings{perCluster}||=2;
  $$settings{oerCluster}||=LONG_MAX; # max size per cluster

  return histogram($settings) if($$settings{histogram});

  $$settings{maxdist}||=die "ERROR: need a maximum distance\n".usage();
  my ($neighbor,$dist)=identifyCloseNeighbors($settings);
  my $cluster=findClusters($neighbor,$dist,$settings);
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
  #$$settings{oerCluster}||=scalar(keys(%dist))**2; # a large number that can't be reached
  #die Dumper $settings;
  return (\%neighbor,\%dist) if wantarray;
  return \%neighbor;
}

sub findClusters{
  my($neighbor,$pDist,$settings)=@_;
  my %neighbor=%$neighbor; # make a copy so that the internal pointer isn't messed up
  my $defaultMaxDist=$$settings{maxdist}; # back up this value

  my @cluster;
  my %seen; # keep track of which nodes have been parsed
  while(my($centerNode,undef)=each(%neighbor)){
    next if($seen{$centerNode}); # %seen is set within findClusterAroundNode()
    $$settings{maxdist}=$defaultMaxDist; # set this back to normal
    my $cluster=findClusterAroundNode($centerNode,\%neighbor,\%seen,$settings);
    next if(!@$cluster);
    unshift(@$cluster,$centerNode);
    
    # if there are too many members in the cluster, make the maxdist more strict or release the cluster back to the wild.
    while(scalar(@$cluster) > $$settings{oerCluster} && $$settings{maxdist} > 50){
      $seen{$_}-- for(@$cluster); # unmark the cluster as seen
      $$settings{maxdist}-=100; 
      $$settings{maxdist}=5 if($$settings{maxdist}<5);
      $cluster=findClusterAroundNode($centerNode,\%neighbor,\%seen, $settings);
      $seen{$_}++ for(@$cluster); # mark it back as seen
    }
    # if there are too few members, just don't report it.
    next if(@$cluster<$$settings{perCluster});

    # the cluster passed all checks: output it.
    my $avgPDist=distances($cluster,$pDist,$settings);
    print join("\t",scalar(@$cluster),$avgPDist,@$cluster)."\n";
    push(@cluster,$cluster);
  }
  return \@cluster;
}

# given a "center" node, find the cluster around it
sub findClusterAroundNode{
  my($centerNode,$neighbor,$seen,$settings)=@_;
  my @cluster=(); # the output

  # Don't process this cluster if you've already processed it.
  # Lock $seen just to make sure!
  { 
    lock $seen;
    return \@cluster if($$seen{$centerNode});
    $$seen{$centerNode}++; # and now we've seen it
    # Failsafe: don't keep processing if you are recursing too deep
    return \@cluster if($$settings{cluster_level}++ > 9999);
  }

  my $neighbors=$$neighbor{$centerNode};
  for my $n(@$neighbors){
    next if($$seen{$n});
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

sub distances{
  my($cluster,$pdist,$settings)=@_;
  my $max=0; my $min=1e999;
  my @data;
  my $num=@$cluster;
  for(my $i=0;$i<$num;$i++){
    my $first=$$cluster[$i];
    for(my $j=$i+1;$j<$num;$j++){
      my $second=$$cluster[$j];
      my $pd=$$pdist{$first}{$second} || $$settings{maxdist};
      push(@data,$pd);

      $max=$pd if($max < $pd);
      $min=$pd if($min > $pd);
    }
  }
  my $avg=int(average(\@data));
  my $stdev=int(stdev(\@data));
  my $plusminus="+/-";
  return "$avg $plusminus $stdev [$min,$max]";
}

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}



## histogram
# cut -f 3 pairwiseDist-MT.txt | perl -lane '$rounded=100*int($F[0]/100);$num{$rounded}++;END{while(($num,$count)=each(%num)){print "$count\t$num";}}' | sort -k2,2n 
sub usage{
  "Produces clusters given a three column file from pairwiseDistances.pl
  genome1  genome2  distance
  Usage: $0 -m maxdist < distances.txt
  -m the maximum distance between two nodes in a cluster.
  -p the minimum number per cluster, before reporting it (default: 2)
  -o the maximum per cluster (default: no maximum)
  --help this help menu
  --histogram to display a histogram and then exit
  "
}
