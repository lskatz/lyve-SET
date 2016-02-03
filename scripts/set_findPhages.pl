#!/usr/bin/env perl
# Author:Lee Katz <lkatz@cdc.gov>
# Thanks: Darlene Wagner for giving me this idea

require 5.12.0;

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse/;
use File::Temp qw/tempdir/;
use List::Util qw/min max/;
use List::MoreUtils qw(uniq);
use Bio::Perl;
use threads;
use Thread::Queue;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg/;
use lib "$FindBin::RealBin/../lib/lib/perl5";
use Array::IntSpan;

local $0=fileparse $0;
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i tempdir=s flanking=i db|database=s));
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("phastXXXXXX",CLEANUP=>1,TMPDIR=>1);
  $$settings{flanking}||=3000;
  $$settings{db}||="$FindBin::RealBin/../lib/phast/phast.faa";
  logmsg "Running blastx against $$settings{db}";

  my $fasta=$ARGV[0];
  die usage() if(!$fasta || $$settings{help});

  logmsg "Tempdir is $$settings{tempdir}";
  my $ranges=phast($fasta,$settings);

  # Print the ranges to stdout
  for my $r(@$ranges){
    print join("\t",@$r)."\n";
  }

  return 0;
}

sub phast{
  my($fasta,$settings)=@_;

  # longest gene in phast is 8573bp, and all regions produced should have 
  # at least that length just in case.
  my $regions=`makeRegions.pl $fasta --numcpus $$settings{numcpus} --numchunks $$settings{numcpus} --overlapby 9999`;
  die "ERROR: problem with makeRegions.pl" if $?;
  my @regions=split(/\n/,$regions);
  logmsg "Regions are: ".join(", ",@regions);

  # Better parallelization: one fasta entry per cpu.
  # Split the query into multiple files and then figure out
  # how many cpus per blast job we need.
  my %seq;
  my $seqin=Bio::SeqIO->new(-file=>$fasta);
  while(my $seq=$seqin->next_seq){
    $seq{$seq->id}=$seq;
  }
  $seqin->close;

  # Spawn threads to take care of each fasta file
  my $regionQ=Thread::Queue->new(@regions);
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&blastWorker,\%seq,$regionQ,$settings);
    $regionQ->enqueue(undef); # one terminator per thread
  }

  my %range;
  # Join together the threads.
  for(@thr){
    my $tRange=$_->join;
    if($_->error()){
      die "ERROR in TID".$_->tid;
    }
    
    # Figure out ranges from the threads
    # Add these coordinates to ranges
    my %seen; # ranges that have already been added.
    for my $r(@$tRange){
      my($contig,$lo,$hi)=@$r;
      $range{$contig}||=Array::IntSpan->new();

      # Set up the "soft" flanking.
      my $softLo=$lo;
      my $softHi=$hi;
      for(my $i=1;$i<$$settings{flanking};$i+=50){
        my $loTmp=$lo-$i;
        my $hiTmp=$hi+$i;

        if($range{$contig}->lookup($loTmp)){
          $softLo=$loTmp;
        }
        if($range{$contig}->lookup($hiTmp)){
          $softHi=$hiTmp;
        }
      }

      ## Don't add this range if was also added.
      ## If we can do anything to avoid ->addrange(),
      ## then we should. It is slow.
      ##next if($seen{$lo}{$hi}++);

      $range{$contig}->set_range($softLo,$softHi,1);
    }
  }

  # Range objects to arrays [contig,start,stop]
  my @range;
  while(my($contig,$obj)=each(%range)){
    ##for my $r($obj->rangeList){
    $obj->consolidate();
    for my $r($obj->get_range_list){
      # Phages are generally longer than just a gene, 
      # so skip any puny ranges.
      my $length=$$r[1] - $$r[0] + 1;
      #next if($length < 2000);
      push(@range,[$contig,@$r]);
    }
  }

  return \@range;
}

sub blastWorker{
  my($seqHash,$regionQ,$settings)=@_;
  my $tempdir=tempdir("$$settings{tempdir}/phastXXXXXX");

  my @range;
  my $i=2; # uniq counter for modified blast files
  while(defined(my $region=$regionQ->dequeue)){
    # Write the sequence to a file
    my($contig,$startStop)=split(/:/,$region);
    my($start,$stop)=split(/-/,$startStop);
    open(BLASTIN,'>',"$tempdir/in.fna") or die "ERROR: could not open $tempdir/in.fna: $!";
    print BLASTIN '>'.$$seqHash{$contig}->id."\n".$$seqHash{$contig}->subseq($start,$stop)."\n";
    close BLASTIN;

    # Blast the sequence
    system("blastx -query $tempdir/in.fna -db $$settings{db} -outfmt 6 -num_threads 1 -max_target_seqs 200000 > $tempdir/bls.out");
    die if $?;
    
    # Read the blast results
    open(BLASTOUT,'<',"$tempdir/bls.out") or die "ERROR: could not open $tempdir/bls.out: $!";
    open(BLASTOUTMOD,'>',"$tempdir/bls.$i.out") or die "ERROR: could not open $tempdir/bls.2.out: $!";
    while(<BLASTOUT>){
      my ($contig,$hit,$identity,$length,$gaps,$mismatches,$qstart,$qend,$sstart,$send,$e,$score)=split /\t/;
        # Think of a good minimum identity. Phages are very diverse
        # but you still want to weed out the awful hits. Remember:
        # this is in protein space and so lower numbers are still
        # okay.
        next if($identity < 50 || $e > 0.05);

      # Reevaluate where the coordinates start based on the subseq
      $qstart+=$start;
      $qend+=$start;

      # Record what happened
      print BLASTOUTMOD join("\t",$contig,$hit,$identity,$length,$gaps,$mismatches,$qstart,$qend,$sstart,$send,$e,$score);

      # Get hi and lo coordinates
      my $lo=min($qstart,$qend);
      my $hi=max($qstart,$qend);
      push(@range,[$contig,$lo,$hi]);
    }
    close BLASTOUT;
    close BLASTOUTMOD;

    $i++;
  }

  return \@range;
}

sub usage{
  "Finds phages in a fasta file using phast
  Usage: $0 file.fasta
  --numcpus  1
  --tempdir  /tmp
  --flanking 3000 Give 'soft' edges to ranges. If blast hits are this many
                  nt away from another blast hit, then join the ranges and
                  include any intermediate positions. If ranges cannot be
                  joined, then do not extend them by this flanking length.
  "
}

