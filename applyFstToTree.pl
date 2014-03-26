#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/sum/;
use File::Basename;
use Time::HiRes qw/time/;
use Bio::TreeIO;
use Statistics::Descriptive;

use threads;
use Thread::Queue;

$0=basename $0;
sub logmsg{$|++;print STDERR "TID".threads->tid." $0: @_\n";$|--;}
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(tree=s help probability format=s)) or die;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{format}||="newick";
  
  my $fst=readFst(\@ARGV,$settings);
  printTree($fst,$$settings{tree},$settings);
  return 0;
}

sub readFst{
  my($tsv,$settings)=@_;
  my %Fst;
  for(@$tsv){
    open(TSV,$_) or die "ERROR: Could not read TSV $_: $!";
    while(my $line=<TSV>){
      chomp $line;
      my($fst,$taxa)=split /\t/,$line;
      my @taxon=sort {$a cmp $b} split(/,/,$taxa);
      my $key=join(",",@taxon);
      $Fst{$key}=$fst;
    }
  }

  transformToPValues(\%Fst,$settings) if($$settings{probability});

  $_=int($_*100) for(values %Fst);

  return \%Fst;
}

# TODO transform in a method that pays attention to the numbers in each group
sub transformToPValues{
  my($Fst,$settings)=@_;
  my $zscoreFactory = Statistics::Zed->new(
    ccor=>1,
    tails=>1,
    precision_s=>2,
    precision_p=>4,
  );
  my $stat=Statistics::Descriptive::Full->new;
  $stat->add_data(values %$Fst);
  my $mean=$stat->mean;
  my $stdev=$stat->standard_deviation;
  
  while(my($key,$value)=each(%$Fst)){
    my $Z=($value-$mean)/$stdev;
    my $p=$zscoreFactory->z2p($Z);
    $$Fst{$key}=$p;
    #logmsg "$value => $Z => $p";
  }
}

sub printTree{
  my($fst,$tree,$settings)=@_;
  my $in=new Bio::TreeIO(-file=>$tree);
  my $out=Bio::TreeIO->new(-format=>$$settings{format});
  while(my $treeIn = $in->next_tree){
    for my $node($treeIn->get_nodes){
      # get all descendents of the node and alpha-sort them
      my @d=sort{$a->id cmp $b->id} $node->each_Descendent;
      next if($node->is_Leaf);
      # make a comma-separated key
      my $key; $key.=$_->id."," for(@d);
      $key=~s/,$//;
      logmsg $key;
      # next if that key isn't in $fst: not all $fst combos were calculated.
      next if(!$$fst{$key});
      # apply the Fst with that key
      logmsg $node->bootstrap;
      $node->bootstrap($$fst{$key});
      $node->id($$fst{$key});
      #logmsg $node->bootstrap;
    }
    $out->write_tree($treeIn);
  }
}
      

sub usage{
  "Applies Fst to a tree and replaces any bootstrap values it might have
  Usage: $0 fst.tsv -t tree.dnd > new.dnd
  fst.tsv is a tab-delimited file with two columns: fst and a comma-separated list of taxa
  -t tree.dnd is a tree to apply values to
  -p To apply normal-distribution p-values instead of Fst
  -f format (default: newick)
  "
}
