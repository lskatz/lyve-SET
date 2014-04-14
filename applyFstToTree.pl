#!/usr/bin/env perl
# Applies the fixation index to a tree

use strict;
use warnings;
use Data::Dumper;
use Bio::TreeIO;
use Getopt::Long;
use File::Basename;
use List::Util qw/min max sum/;
use FindBin;
use Statistics::Descriptive;
use Math::Round qw/nearest/;

use lib "$FindBin::Bin/lib";
use Statistics::Data;
use Statistics::Zed;

my $zed=Statistics::Zed->new(tails=>1,precision_p=>4,precision_s=>4);
$0=fileparse $0;
sub logmsg{print STDERR (caller(1))[3].": @_\n";}
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(tree=s pairwise=s outprefix=s fst=s help numcpus=i reps=i outputType=s)) or die $!;
  die usage($settings) if($$settings{help});
  $$settings{numcpus}||=1;
  $$settings{reps}||=100;
  $$settings{outputType}||="samples";
  
  # required parameters
  for(qw(tree pairwise outprefix)){
    die "ERROR: need parameter for $_\n".usage() if(!$$settings{$_});
  }
  # check for file existence
  for(qw(tree pairwise)){
    die "ERROR: file $_ does not exist!" if(!-e $$settings{$_});
  }

  my ($fstAvg,$fstStdev);
  if($$settings{fst}){
    logmsg "Reading Fsts to get a background distribution.";
    ($fstAvg,$fstStdev)=readFst($$settings{fst},$settings);
  }
  logmsg "Reading the pairwise file";
  my $pairwise=readPairwise($$settings{pairwise},$settings);
  logmsg "Reading the tree and using the pairwise distances to calculate Fst on the tree";
  print join("\t",qw(Fst num1 num2 group1 group2))."\n"; # headers for stdout
  my ($fstTreeArr,$fstSample)=applyFst($$settings{tree},$pairwise,$fstAvg,$fstStdev,$settings);
  #my $pvalueTreeArr=applyFst($$settings{tree},$pairwise,$fstAvg,$fstStdev,{%$settings,transform=>1});

  outputTrees([$fstTreeArr],$$settings{outprefix},$settings);
  outputStatisticsThresholds($fstSample,$$settings{outprefix},$settings);

  return 0;
}

sub readFst{
  my($file,$settings)=@_;
  my @fst;
  open(FST,$file) or die "ERROR: could not read fst file:$!";
  while(<FST>){
    chomp;
    push(@fst,(split(/\t/,$_))[0]);
  }
  close FST;
  my $stat=Statistics::Descriptive::Full->new;
  $stat->add_data(@fst);
  return($stat->mean,$stat->standard_deviation);
}

# Read the pairwise file.
sub readPairwise{
  my($pairwise,$settings)=@_;
  my %dist;
  open(PAIRWISE,$pairwise) or die "ERROR: could not read $pairwise:$!";
  while(<PAIRWISE>){
    chomp;
    my($g1,$g2,$dist)=split /\t/;
    $dist{$g1}{$g2}=$dist;
    $dist{$g2}{$g1}=$dist;
  }
  close PAIRWISE;
  return \%dist;
}

# Count how many ancestor nodes are in a tree
sub countAncNodes{
  my($infile,$settings)=@_;
  my $in=new Bio::TreeIO(-file=>$infile);
  my $i=0;
  my $numTrees=0;
  while(my $treeIn=$in->next_tree){
    $numTrees++;
    for my $node($treeIn->get_nodes){
      next if($node->is_Leaf);
      next if(defined($node->bootstrap) && $node->bootstrap < 0.7);
      $i++;
    }
  }
  $in->close; # make sure to close it since it will be opened again soon

  logmsg "WARNING: the number of trees in the input file is $numTrees, whereas I expected only one." if($numTrees > 1);
  return $i;
}

sub applyFst{
  my($infile,$pairwise,$fstAvg,$fstStdev,$settings)=@_;
  $$settings{transform}||=0;
  my (@fstTree);
  my $i=0;
  my @ID=keys(%$pairwise); # all taxa

  if($$settings{transform}){
    logmsg "Making a tree by calculating new Fst values and transforming them to p-values";
  } else {
    logmsg "Making a tree with Fst";
  }

  logmsg "Counting nodes in the tree before analyzing...";
  my $numAncestorNodes=countAncNodes($infile,$settings);
  logmsg "I estimate about $numAncestorNodes ancestor nodes (not leaves) that will be analyzed";

  logmsg "Defining all possible outgroups";
  my %group2=defineGroup2s($infile,$pairwise,$settings);

  logmsg "Calculating Fst for each high-confidence node";
  my $in=new Bio::TreeIO(-file=>$infile);
  my %fstSample;
  while(my $treeIn=$in->next_tree){
    for my $node($treeIn->get_nodes){
      next if($node->is_Leaf);

      # get all descendents of the node and alpha-sort them
      my @d;
      for my $taxon($node->get_all_Descendents){
        my $taxonId=$taxon->id;
        next if(!$taxonId);
        push(@d,$taxon) if(defined $$pairwise{$taxonId});
      }
      @d=sort{$a->id cmp $b->id} @d;
      my $numDesc=@d;

      # The bootstrap is the boostrap if it exists, or it is the node identifier because of Newick weirdness.
      my $bootstrap=$node->bootstrap || $node->id;
         $bootstrap=0 if(!$bootstrap);
      # If the bootstrap is low then don't bother doing fst on this node -- it's not supported!
      #if($bootstrap < 70){
      #  logmsg "Skipping node '$bootstrap' due to the low bootstrap values";
      #  $node->bootstrap(1); # give it a really low fst
      #  $node->id(1);        # give it a really low fst
      #  next;
      #}

      # find Fst
      my @group1=map($_->id,@d);
      # min/max allowed in a group that group1 is randomly compared against. 5% down/5% up
      my($minNum,$maxNum)=(nearest(1,$numDesc * 0.95),nearest(1,$numDesc / 0.95));
      # create a sample set of all these groups
      my @clades;
      for($i=$minNum;$i<$maxNum;$i++){
        push(@clades,@{ $group2{$i} }) if(defined($group2{$i}));
      }
      my($avg,$stdev,$rawData)=fstDistribution(\@group1,\@clades,$pairwise,$settings);
      push(@{ $fstSample{$numDesc} },@$rawData);
      my $identifier=sprintf("%0.2f",$avg);

      # transform the Fst to p-value if background distribution is set
      #if(defined($fstAvg) && defined($fstStdev)){
      #  my $z=(($avg) - ($fstAvg))/$fstStdev;
      #  my $p=$zed->z2p($z,tails=>1);
      #  $identifier=$p;
      #}

      # apply the Fst with that key
      $node->bootstrap($identifier);
      $node->id($identifier);

      if(++$i % 100 ==0){
        logmsg "Finished with $i ancestor nodes";
      }
      print join("\t",sprintf("%0.4f",$avg),scalar(@group1),scalar(@group1),join(",",@group1),'-')."\n" if($$settings{outputType} eq 'averages');
    }
    push(@fstTree,$treeIn);
  }
  return (\@fstTree,\%fstSample);
}

sub outputTrees{
  my($treeArr,$prefix,$settings)=@_;
  my $numTrees=@{$$treeArr[0]};
  my $fstOut=new Bio::TreeIO(-file=>">$prefix.fst.dnd");
  #my $pOut=  new Bio::TreeIO(-file=>">$prefix.pvalue.dnd");
  for(my $i=0;$i<$numTrees;$i++){
    my $fstTree=$$treeArr[0][$i];
    $fstOut->write_tree($fstTree);

    my $pvalueTree=$$treeArr[1][$i] || next;
    #$pOut->write_tree($pvalueTree);
  }
  $fstOut->close;
  #$pOut->close;
  return 1;
}

sub outputStatisticsThresholds{
  my($fstSample,$prefix,$settings)=@_;
  my $stat=Statistics::Descriptive::Full->new;

  my $out="$prefix.stats.tsv";
  my ($min,$max)=(min(keys(%$fstSample)),max(keys(%$fstSample)));
  open(STATS,">",$out) or die "ERROR: could not write to file $out:$!";
  print STATS join("\t",qw(numInGroup threshold average standardDeviation))."\n";
  for(my $i=$min;$i<=$max;$i++){
    next if(!$$fstSample{$i} || @{$$fstSample{$i}} < 2);
    $stat->add_data(@{$$fstSample{$i}});
    my ($avg,$stdev)=(sprintf("%0.4f",$stat->mean),sprintf("%0.4f",$stat->standard_deviation));
    my $threshold=sprintf("%0.4f",$stat->percentile(95));
    my $threshold2=sprintf("%0.4f",$stat->percentile(75));
    print STATS join("\t",$i,"$threshold/$threshold2",$avg,$stdev)."\n";
  }
  close STATS;
  return 1;
}

# Return a group with the same number of elements as group1.
# Will be randomized but will not have repeat elements nor will it intersect with group1.
sub group2{
  my($inGroup,$totalGroup,$settings)=@_;
  my %inGroup;
  $inGroup{$_}=1 for(@$inGroup); # index
  my @outGroup;
  my $ingroupNum=@$inGroup;
  my $indexFromTotal=@$totalGroup-1;
  my $j=0;
  for(my $i=0;$i<=$indexFromTotal;$i++){
    my $randEl=$$totalGroup[rand($indexFromTotal)];
    if(!$inGroup{$randEl}){
      push(@outGroup,$randEl);
      last if(++$j >= $ingroupNum);
      $inGroup{$randEl}=1; # this will make it be random without replacement
    } else {
      #logmsg "WARNING ($i): $randEl is in inGroup:".Dumper $inGroup;
    }
  }
  return \@outGroup;
}

# Find all phylogroups; bin them by the number of members.
# Returns @arr=>(3=>[[g1,g2,g3],[g4,g5,g6]],...4=>[...
#  Array=>[numInGroup=>[\@group1,\@group2...],numInGroup2=>[...
sub defineGroup2s{
  my($tree,$pairwise,$settings)=@_;
  my %group;
  my $in=Bio::TreeIO->new(-file=>$tree);
  while(my $treeIn=$in->next_tree){
    for my $node($treeIn->get_nodes){
      next if($node->is_Leaf);

      # get all descendents of the node and alpha-sort them
      my @d;
      for my $taxon($node->get_all_Descendents){
        my $taxonId=$taxon->id;
        next if(!$taxonId);
        push(@d,$taxon) if(defined $$pairwise{$taxonId});
      }
      @d=sort{$a->id cmp $b->id} @d;
      @d=map{$_->id} @d;
      my $numDesc=@d;
      next if(!$numDesc || $numDesc < 2 || !@d);

      # Define the within distance in the last element
      my $dist=averageGroupDistance(\@d,\@d,$pairwise,1,$settings);
      push(@d,$dist);
      push(@{$group{$numDesc}},\@d);
    }
  }
  #print Dumper \%group;die;
  $in->close;
  return %group;
}

# find an average and standard deviation of an Fst of a group1
sub fstDistribution{
  my($group1,$allTaxa,$distance,$settings)=@_;
  die "ERROR: reps can't be less than 1" if($$settings{reps}<1);
  my $stat=Statistics::Descriptive::Full->new;
  my $group1Str=join("",@$group1);

  my @fst;
  return(0.01,0,\@fst) if(@$allTaxa < 5);# no point in making fst if there are only a few other groups
  my $within1=averageGroupDistance($group1,$group1,$distance,1,$settings);
  for(my $i=0;$i<$$settings{reps};$i++){
    my @group2=@{ $$allTaxa[rand(@$allTaxa - 1)] };
    if($group1Str eq join("",@group2)){
      $i--;
      next;
    }
    my $within2=pop(@group2);
    my $between=averageGroupDistance($group1,\@group2,$distance,0,$settings);
    my $within=($within1+$within2)/2;
    my $fst=($between - $within)/$between;
    push(@fst,$fst);
    print join("\t",sprintf("%0.4f",$fst),scalar(@$group1),scalar(@group2),join(",",@$group1),join(",",@group2))."\n" if($$settings{outputType} eq 'samples');
  }
  $stat->add_data(@fst);
  return($stat->mean,$stat->standard_deviation,\@fst);
}

sub fst{
  my($group1,$group2,$distance,$settings)=@_;
  #warn Dumper [$group1,$group2];
  my $within1=averageGroupDistance($group1,$group1,$distance,1,$settings);
  my $within2=averageGroupDistance($group2,$group2,$distance,1,$settings);
  # for some reason, this takes a very long time to run when multithreaded
  my $between=averageGroupDistance($group1,$group2,$distance,0,$settings);

  my $numGroup1=@$group1;
  my $numGroup2=@$group2;
  my $total=$numGroup1+$numGroup2;
  my $within=($within1*$numGroup2/$total + $within2*$numGroup1/$total)/2;
  return 0 if(!$between);
  my $fst=($between-$within)/$between;
  return $fst;
}

# specify the $is_same parameter to note that the two groups are really the same group.
sub averageGroupDistance{
  my($id1,$id2,$distance,$is_same,$settings)=@_;
  my @distance;
  my %seen;
  for(my $i=0;$i<@$id1;$i++){
    my $j=0;
    $j=$i+1 if($is_same);
    for(my $j=$j;$j<@$id2;$j++){
      next if($$id1[$i] eq $$id2[$j]);
      die "ERROR: $$id1[$i] - $$id2[$j] distance was not found" if(!defined($$distance{$$id1[$i]}{$$id2[$j]}));
      push(@distance,$$distance{$$id1[$i]}{$$id2[$j]});
    }
  }
  return 0 if(!@distance);
  return sum(@distance)/@distance;
}


sub usage{
  my($settings)=@_;
  my $help="Applies the fixation index (Fst) to a tree using pairwise distances.
  Usage: $0 -t in.dnd -p pairwise.tsv --outprefix out --outputType averages > fstaverages.tsv
  Outputs groups with fst > 0.7
  -t in.dnd The tree file, which will be parsed by BioPerl. Format determined by extension.
  -p pairwise.tsv The pairwise distances file. NOTE: you can get pairwise distances from pairwiseDistances.pl
  --outprefix prefix The output prefix. Prefix can have a directory name encoded in it also, e.g. 'tmp/prefix'
    Output files are: prefix.fst.dnd, prefix.pvalue.dnd
  -h for additional help";
  return $help if(!$$settings{help}); # return shorthand help
  $help.="
  OPTIONAL
  --outputType samples The type of stdout you want.  samples==every random sample; averages==average per node; >/dev/null for no output
  --numcpus 1 The number of cpus to use (not currently multithreaded anyway)
  --reps 100 The number of repetitions for each node to calculate an average Fst
  Example: $0 applyFstToTree.pl -p pairwise.tsv -t 4b.dnd --outprefix out --outputType samples > fstsampling.tsv
  ";
  return $help;
}
