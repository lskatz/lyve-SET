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

use lib "$FindBin::Bin/lib";
#use Statistics::Zed;

$0=fileparse $0;
sub logmsg{print STDERR (caller(1))[3].": @_\n";}
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(tree=s pairwise=s outprefix=s fst=s help numcpus=i)) or die $!;
  die usage() if($$settings{help});
  $$settings{numcpus}||=1;
  
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
  my ($fstTreeArr,$pvalueTreeArr)=applyFst($$settings{tree},$pairwise,$fstAvg,$fstStdev,$settings);

  outputTrees([$fstTreeArr,$pvalueTreeArr],$$settings{outprefix},$settings);

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
  my (@fstTree,@pvalueTree);
  my $i=0;
  my @ID=keys(%$pairwise); # all taxa

  logmsg "Counting nodes in the tree before analyzing...";
  my $numAncestorNodes=countAncNodes($infile,$settings);
  logmsg "I estimate about $numAncestorNodes ancestor nodes (not leaves) that will be analyzed";

  logmsg "Defining all possible outgroups";
  my %group2=defineGroup2s($infile,$pairwise,$settings);

  logmsg "Calculating Fst for each high-confidence node";
  my $in=new Bio::TreeIO(-file=>$infile);
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

      # If the bootstrap is low then don't bother doing fst on this node -- it's not supported!
      if(defined($node->bootstrap) && $node->bootstrap < 0.7){
        logmsg "Skipping ".$node->bootstrap." due to the low bootstrap values";
        $node->bootstrap(0.01); # give it a really low fst
        next;
      }

      # find Fst
      my @group1=map($_->id,@d);
      # min/max allowed in a group that group1 is randomly compared against. 5% down/5% up
      my($minNum,$maxNum)=(int($numDesc/10*9.5),int($numDesc/9.5*10));
      # create a sample set of all these groups
      my @clades;
      for($i=$minNum;$i<$maxNum;$i++){
        push(@clades,@{ $group2{$i} }) if(defined($group2{$i}));
      }
      my($avg,$stdev)=fstDistribution(\@group1,\@clades,$pairwise,{reps=>100});
      my $identifier=sprintf("%0.2f",$avg);

      # transform the Fst to p-value if background distribution is set
      if(defined($fstAvg) && defined($fstStdev)){
        my $z=(($avg) - ($fstAvg))/sqrt(($fstStdev*$fstStdev)/100);
        #my $p=$zed->p_value($z);
        #$identifier=$p;
      }

      # apply the Fst with that key
      my $bootstrap=$node->bootstrap || $node->id;
         $bootstrap=0 if(!$bootstrap);
      $node->bootstrap($identifier);
      $node->id($identifier);

      # The Fst is potentially significant if greater than 50%.
      # Otherwise leave it up to the user to sort or filter these results by using stdout.
      print join("\t",sprintf("%0.4f",$avg),$bootstrap,join(",",@group1))."\n" if($avg > 0.5);

      if(++$i % 100 ==0){
        logmsg "Finished with $i ancestor nodes";
      }
    }
    push(@fstTree,$treeIn);
  }
  return (\@fstTree,\@pvalueTree);
}

sub outputTrees{
  my($treeArr,$prefix,$settings)=@_;
  my $numTrees=@{$$treeArr[0]};
  my $fstOut=new Bio::TreeIO(-file=>">$prefix.fst.dnd");
  my $pOut=  new Bio::TreeIO(-file=>">$prefix.pvalue.dnd");
  for(my $i=0;$i<$numTrees;$i++){
    my $fstTree=$$treeArr[0][$i];
    my $pvalueTree=$$treeArr[1][$i] || $$treeArr[0][$i];
    $fstOut->write_tree($fstTree);
    $pOut->write_tree($pvalueTree);
  }
  $fstOut->close;
  $pOut->close;
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

sub fstDistribution{
  my($group1,$allTaxa,$distance,$settings)=@_;
  $$settings{reps}||=100;
  die "ERROR: reps can't be less than 1" if($$settings{reps}<1);
  my $stat=Statistics::Descriptive::Full->new;

  my @fst;
  my $within1=averageGroupDistance($group1,$group1,$distance,1,$settings);
  for(my $i=0;$i<$$settings{reps};$i++){
    last if(@$allTaxa < 2); # no point in making fst if there is no other group
    my @group2=@{ $$allTaxa[rand(@$allTaxa - 1)] };
    if(join("",@$group1) eq join("",@group2)){
      $i--;
      next;
    }
    my $within2=pop(@group2);
    my $between=averageGroupDistance($group1,\@group2,$distance,0,$settings);
    my $within=($within1+$within2)/2;
    my $fst=($between - $within)/$between;
    #die Dumper ["($within1 $within2)=>$within $between => $fst",$group1,\@group2] if(@group2 < 30);
    #my $fst=fst($group1,\@group2,$distance);
    push(@fst,$fst);
  }
  return(0.01,0) if(@fst < 2);
  $stat->add_data(@fst);
  return($stat->mean,$stat->standard_deviation);
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
  "Applies the fixation index (Fst) to a tree using pairwise distances.
  Usage: $0 -t in.dnd -p pairwise.tsv -o prefix > fstgroups.tsv
  Outputs groups with fst > 0.7
  -t in.dnd The tree file, which will be parsed by BioPerl. Format determined by extension.
  -p pairwise.tsv The pairwise distances file. NOTE: you can get pairwise distances from pairwiseDistances.pl
  -o prefix The output prefix. Prefix can have a directory name encoded in it also, e.g. 'tmp/prefix'
    Output files are: prefix.fst.dnd, prefix.pvalue.dnd
  OPTIONAL
  -f fst.tsv A background fixation index file, generated from fstFromPairwise.pl
  --numcpus 1 The number of cpus to use
  "
}
