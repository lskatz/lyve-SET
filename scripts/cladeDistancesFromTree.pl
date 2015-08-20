#!/usr/bin/env perl
# Labels nodes with significant bootstrap values and returns useful pairwise distance statistics

use strict;
use warnings;
use Data::Dumper;
use Bio::TreeIO;
use Getopt::Long;
use File::Basename;
use List::Util qw/min max sum/;
use FindBin;
use Math::Round qw/nearest/;

use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/lib/perl5";
use Statistics::Descriptive;
use Statistics::Basic qw(median);
use LyveSET qw(logmsg);

$0=fileparse $0;
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(tree=s pairwise=s outprefix=s help outputType=s bootstrap=i nameTrim=i)) or die $!;
  die usage($settings) if($$settings{help});
  $$settings{outputType}||="merged";
  $$settings{bootstrap}||=70;
  $$settings{nameTrim}||=30;
  
  # required parameters
  for(qw(tree pairwise outprefix)){
    die "ERROR: need parameter for $_\n".usage() if(!$$settings{$_});
  }
  # check for file existence
  for(qw(tree pairwise)){
    die "ERROR: file $$settings{$_} does not exist!" if(!-e $$settings{$_});
  }

  logmsg "Reading the pairwise file";
  my $pairwise=readPairwise($$settings{pairwise},$$settings{nameTrim},$settings);
  
  logmsg "Reading the tree and assigning names to high bootstrap nodes";
  my ($treeArr, $labelHash, $allLeafArr) = readSigClades($$settings{tree},$$settings{bootstrap},$$settings{nameTrim},$settings);

  logmsg "Reading high bootstrap nodes and outputing statistics";
  my ($pairwiseCladeStats, $singletonCladeStats, $singletonTaxaStats, $treeStats) = generateSigCladeStats($labelHash, $pairwise, $allLeafArr, $settings);
  
  outputTrees([$treeArr],$$settings{outprefix},$settings);
  outputStatistics($pairwiseCladeStats,$singletonCladeStats,$singletonTaxaStats,$treeStats,$$settings{outprefix},$settings);
  outputLabels($labelHash,$$settings{outprefix},$settings);

  return 0;
}

# Read the pairwise file.
sub readPairwise{
  my($pairwise,$nameTrim,$settings)=@_;
  my %dist;
  open(PAIRWISE,$pairwise) or die "ERROR: could not read $pairwise:$!";
  while(<PAIRWISE>){
    chomp;
    my($g1,$g2,$dist)=split /\t/;
    $g1=substr($g1,0,$nameTrim);
    $g2=substr($g2,0,$nameTrim);
    next if($g1 eq $g2);
    
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
  my($infile,$bootstrapthreshold,$nameTrim,$settings)=@_;
  my @bootTree;	# Tree with significant nodes labeled with letters
  my $i=0;
  my %labels = ();	# Descendent taxa leaf name lists keyed to each significant node's label
  my @allLeaves;	# List of all leaf ids in tree
  my @labelList;	# Each leaf and its leaf names

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
        $taxonId=substr($taxonId,0,$nameTrim);
        push(@d,$taxonId); #if(defined $$pairwise{$taxonId});
      }
      @d=sort{$a cmp $b} @d;
      my $numDesc=@d;
      
      # The bootstrap is the boostrap if it exists, 
      # 	or it is the node identifier because of Newick weirdness.
      my $bootstrap=$node->bootstrap || $node->id;
         $bootstrap=0 if(!$bootstrap);
      # If the bootstrap is low then skip
      if($bootstrap < $bootstrapthreshold){
        logmsg "Skipping node '$bootstrap' due to the low bootstrap values";
        $node->bootstrap(''); # give it a really low fst
        $node->id('');        # give it a really low fst
        next;
      }

	 # Relabel significant tree nodes with node name and bootstrap value
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
	my @treeSummary;
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
	
	# 4.	Calculate overall distances across the entire tree
	my $treeDistances = distancesBetweenNodes( $allLeafArrRef, $allLeafArrRef, 
						   				$pairwiseRef, $settings );
	push(@treeSummary, sprintf("%s\t%s\t%d\t%d\t%d\t%d",
						"Tree", "Tree", median(@$treeDistances),
						min(@$treeDistances), max(@$treeDistances),
						medianAbsoluteDeviation($treeDistances)));
	return (\@pairwiseDistances, \@nodeSingletonDistances, \@taxaSingletonDistances,\@treeSummary);
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
  my $fstOut=new Bio::TreeIO(-file=>">$prefix.labelednodes.dnd");
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

sub outputStatistics{
  my($pairwiseCladeStats,$singletonCladeStats,$singletonTaxaStats,$treeStats,$prefix,$settings)=@_;
  
  if($$settings{outputType} eq "merged"){
  	my $out="$prefix.merged.cladestats.tsv";
  
  	open(STATS,">",$out) or die "ERROR: could not write to file $out:$!";
  	print STATS join("\t",qw(Node1 Node2 median min max MAD))."\n";
	print STATS join("\n",@$pairwiseCladeStats)."\n";
	print STATS join("\n",@$singletonCladeStats)."\n";
	print STATS join("\n",@$singletonTaxaStats)."\n";
	print STATS join("\n",@$treeStats)."\n";
  	close STATS;
  }else{
  	my $outpairwise="$prefix.pairwise.cladestats.tsv";
  	my $outsingClade="$prefix.singletonClades.cladestats.tsv";
  	my $outsingTaxa="$prefix.singletonTaxa.cladestats.tsv";
   	my $outTree="$prefix.tree.cladestats.tsv";
  
  	open(STATS1,">",$outpairwise) or die "ERROR: could not write to file $outpairwise:$!";
  	print STATS1 join("\t",qw(Node1 Node2 median min max MAD))."\n";
	print STATS1 join("\n",@$pairwiseCladeStats)."\n";
	close STATS1;
	
	open(STATS2,">",$outsingClade) or die "ERROR: could not write to file $outsingClade:$!";
  	print STATS2 join("\t",qw(Node1 Node2 median min max MAD))."\n";
	print STATS2 join("\n",@$singletonCladeStats)."\n";
	close STATS2;
	
	open(STATS3,">",$outsingTaxa) or die "ERROR: could not write to file $outsingTaxa:$!";
  	print STATS3 join("\t",qw(Node1 Node2 median min max MAD))."\n";
	print STATS3 join("\n",@$singletonTaxaStats);
  	close STATS3;
  	
  	open(STATS4,">",$outTree) or die "ERROR: could not write to file $outTree:$!";
  	print STATS4 join("\t",qw(Node1 Node2 median min max MAD))."\n";
	print STATS4 join("\n",@$treeStats) ."\n";
  	close STATS4;
  }
  return 1;
}

sub outputLabels{
  	my($labelHashRef,$prefix,$settings)=@_;
  
  	my $out="$prefix.labels.tsv";
  	open(LABELS,">",$out) or die "ERROR: could not write to file $out:$!";
  	print LABELS join("\t",qw(Label Leaves))."\n";
  	while ( my ( $key, $value ) = each %$labelHashRef ){
  		printf LABELS "%s: %s\n", $key, join(", ", @$value);
  	}
  	close LABELS;
  	return 1;
}

sub usage{
  my($settings)=@_;
  my $help="Labels nodes with significant bootstrap values and returns useful
   pairwise distance statistics for those nodes and the tree.
  
  Usage: $0 -t in.dnd -p pairwise.tsv --outprefix out
  
  Outputs Newick tree, list of labels and leaves, and tab delimited files
	with pairwise distance statistics
  
  -t in.dnd The tree file, which will be parsed by BioPerl. 
  	Format determined by extension.
  -p pairwise.tsv The pairwise distances file. NOTE: you can get pairwise 
  	distances from pairwiseDistances.pl
  --outprefix prefix The output prefix. Prefix can have a directory name 
  	encoded in it also, e.g. 'tmp/prefix'.
    	Output files are: prefix.labelednodes.dnd, prefix.merged.stats.tsv,
    	and prefix.labels.tsv 
    	(using the split option will output: prefix.pairwise.stats.tsv, 
    	prefix.singletonClades.stats.tsv, prefix.singletonTaxa.stats.tsv, 
    	and prefix.tree.stats.tsv )
  -h for additional help";
  return $help if(!$$settings{help}); # return shorthand help
  $help.="
  OPTIONAL
  --outputType merged (default) The type of distance stats output you want. 
  	'merged' for one big file (prefix.merged.stats.tsv)
  	'split' for separate files for pairwise by node (prefix.pairwise.stats.tsv), 
  	individual nodes vs other leaves (prefix.singletonClades.stats.tsv), 
  	individual taxa vs other leaves (singletonTaxa.stats.tsv), and overall 
  	stats for the entire tree (prefix.tree.stats.tsv)
  --bootstrap 70 (default) The percent threshold for significant bootstrap 
  	values. Only nodes with bootstrap values higher than this number will 
  	be included in pairwise by node and individual node statistic calculations.
  --nameTrim 30 (default) The character limit for taxa sample ID names both in
     the pairwise distance file and tree file. If a name in either file is 
     longer than this value, it will be trimmed.
  
  Example: 
  $0 applyFstToTree.pl -p pairwise.tsv -t 4b.dnd --outprefix out --outputType split --bootstrap 80 --nameTrim 20
  ";
  return $help;
}
