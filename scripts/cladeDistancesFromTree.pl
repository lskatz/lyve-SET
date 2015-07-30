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
use Statistics::Basic qw(:all);

use lib "$FindBin::Bin/../lib";

$0=fileparse $0;
sub logmsg{print STDERR (caller(1))[3].": @_\n";}
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(tree=s pairwise=s outprefix=s help numcpus=i reps=i outputType=s)) or die $!;
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
    die "ERROR: file $$settings{$_} does not exist!" if(!-e $$settings{$_});
  }

  logmsg "Reading the pairwise file";
  my $pairwise=readPairwise($$settings{pairwise},$settings);
  
  logmsg "Reading the tree and assigning names to high bootstrap nodes";
  my ($treeArr, $labelHash, $allLeafArr) = readSigClades($$settings{tree}, $settings);

  logmsg "Reading high bootstrap nodes and outputing statistics";
  my ($pairwiseCladeStats, $singletonCladeStats, $singletonTaxaStats) = generateSigCladeStats($labelHash, $pairwise, $allLeafArr, $settings);
  
  outputTrees([$treeArr],$$settings{outprefix},$settings);
  #outputStatisticsThresholds($fstSample,$$settings{outprefix},$settings);

  return 0;
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

# 	Generates unique node labels for significant nodes (>70%)
#	For each significant unique node label, generates a list of descendent leaves 
# Input:		Newick tree with bootstrap values
# Returns : 	Newick tree relabeled with unique node labels
#			Hash keyed to unique node labels pointing to lists of descendent leaves for each label
sub readSigClades{
  my($infile,$settings)=@_;
  my @bootTree;
  my $i=0;
  my %labels = ();	# Descendent taxa leaf name lists keyed to each significant node's label
  my @allLeaves;	# List of all leaf ids in tree

  logmsg "Counting nodes in the tree before analyzing...";
  my $numAncestorNodes=countAncNodes($infile,$settings);
  logmsg "I estimate about $numAncestorNodes ancestor nodes (not leaves) that will be analyzed";

  logmsg "Compiling leaf list for each high-confidence node";
  my $in=new Bio::TreeIO(-file=>$infile);
  while(my $treeIn=$in->next_tree){
    my $labelLetter = 'A';
    for my $node($treeIn->get_nodes){
      next if($node->is_Leaf);

      # get all leaves of the node and alpha-sort them
      my @d;
      for my $taxon($node->get_all_Descendents){
        my $taxonId=$taxon->id;
        next if(!$taxonId);
        next if(!$taxon->is_Leaf);
        push(@d,$taxonId); #if(defined $$pairwise{$taxonId});
      }
      @d=sort{$a cmp $b} @d;
      my $numDesc=@d;
      
      # The bootstrap is the boostrap if it exists, 
      # 	or it is the node identifier because of Newick weirdness.
      my $bootstrap=$node->bootstrap || $node->id;
         $bootstrap=0 if(!$bootstrap);
      # If the bootstrap is low then skip
      if($bootstrap < 70){
        logmsg "Skipping node '$bootstrap' due to the low bootstrap values";
        $node->bootstrap(1); # give it a really low fst
        $node->id(1);        # give it a really low fst
        next;
      }

	 # Relabel significant tree nodes with node name and bootstrap value
	 #my $labelName = sprintf "%s %s", $labelLetter, $bootstrap;
	 my $labelName = $labelLetter;
      $node->bootstrap($labelName);
      $node->id($labelName);

      # Insert descendent node list into label hash
      $labels{ $labelLetter } = \@d;
	 
	 ++$labelLetter;
      
      if(++$i % 100 ==0){
        logmsg "Finished with $i ancestor nodes";
      }
    }
    push(@bootTree,$treeIn);
    
    logmsg "Compiling list of all leaves in tree.";
    foreach my $leaf ($treeIn->get_leaf_nodes){
    		push(@allLeaves, $leaf->id);
    }
    my $numLeaves = @allLeaves;
    logmsg "There are $numLeaves leafs in this tree.";
  }
  
  return (\@bootTree,\%labels,\@allLeaves);
}

# Iterates across nodes and generates distance metrics pairwise 
#	with other nodes, and between the node and everything else
sub generateSigCladeStats{
	my ($labelHashRef,$pairwiseRef,$allLeafArrRef, $settings)=@_;
	my @pairwiseDistances;
	my @nodeSingletonDistances;
	my @taxaSingletonDistances;
	
	logmsg "Returning significant clade stats";
	foreach my $node1 (keys %$labelHashRef) {
		# 1. Between and within node distances
		foreach my $node2 (keys %$labelHashRef) {
			my $distances = distancesBetweenNodes( $$labelHashRef{$node1}, 
											$$labelHashRef{$node2}, 
							   				$pairwiseRef, $settings );
			push(@pairwiseDistances, sprintf("%s\t%s\t%d\t%d\t%d\t%d",
										$node1, $node2, median(@$distances),
										min(@$distances), max(@$distances),
										medianAbsoluteDeviation($distances)));
		}
		
		# 2. Between node leaves and all other leaves
		my @nodeLeafCopy = @{$$labelHashRef{$node1}};
		my @allLeafCopy = @$allLeafArrRef;
		my %diff;
		
		# We set the keys of the diff hash to each label in allLeafCopy
		# Then we delete out all the keys that match labels in nodeLeafCopy
		# And what's left in diff is the difference between the two arrays
		@diff{ @allLeafCopy } = undef;
		delete @diff{ @nodeLeafCopy };
		my @diffLeaf = keys %diff;
		my $distances = distancesBetweenNodes( $$labelHashRef{$node1}, \@diffLeaf, 
						   				$pairwiseRef, $settings );
		push(@nodeSingletonDistances, sprintf("%s\t%s\t%d\t%d\t%d\t%d",
										$node1, "Rest", median(@$distances),
										min(@$distances), max(@$distances),
										medianAbsoluteDeviation($distances)));				   				
	}
	
	# 3. For each leaf, we also calculate distances to all other leafs
	foreach my $taxa (@$allLeafArrRef){
		my @taxaArr = ($taxa);
		my @allLeafCopy = @$allLeafArrRef;
		my %diff;
		
		# We set the keys of the diff hash to each label in allLeafCopy
		# Then we delete out all the keys that match labels in nodeLeafCopy
		# And what's left in diff is the difference between the two arrays
		@diff{ @allLeafCopy } = undef;
		delete @diff{ @taxaArr };
		my @diffLeaf = keys %diff;
		my $distances = distancesBetweenNodes( \@taxaArr, \@diffLeaf, 
						   				$pairwiseRef, $settings );
		push(@taxaSingletonDistances, sprintf("%s\t%s\t%d\t%d\t%d\t%d",
										$taxa, "Rest", median(@$distances),
										min(@$distances), max(@$distances),
										medianAbsoluteDeviation($distances)));
	}
	return (\@pairwiseDistances, \@nodeSingletonDistances, \@taxaSingletonDistances);
}

# Distances from the pairwise distance list between any two leaf lists
sub distancesBetweenNodes{
	my ($leafListRef1, $leafListRef2, $pairwiseRef, $settings)=@_;
	my @distanceArr = ();
	
	# Gather distances between members of list1 and list2
	foreach my $name1 (@$leafListRef1) {
		foreach my $name2 (@$leafListRef2){
			push (@distanceArr, ${$pairwiseRef}{$name1}{$name2}) if $name1 ne $name2;
		}
	}
	return \@distanceArr;
}

# Median absolute deviation (MAD):
#   1. All the absolute distances from the median
#   2. Take the median of those absolute distances
# An LSKatz original
sub medianAbsoluteDeviation{
  my($data,$settings)=@_;
  my $stat=Statistics::Descriptive::Full->new();
  $stat->add_data(@$data);
  my $median=$stat->median();
  my @deviation;
  for(my $i=0;$i<@$data;$i++){
    push(@deviation,abs($$data[$i] - $median));
  }
  my $stat2=Statistics::Descriptive::Full->new();
  $stat2->add_data(@deviation);
  return $stat2->median;
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
