#!/usr/bin/env perl
# Author:Lee Katz <lkatz@cdc.gov>
# Thanks: Darlene Wagner for giving me this idea

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse/;
use File::Temp qw/tempdir/;
use List::Util qw/min max/;
use List::MoreUtils qw(uniq);
use Bio::Perl;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg/;
use lib "$FindBin::RealBin/../lib/lib/perl5";
use Number::Range;

local $0=fileparse $0;
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i tempdir=s flanking=i));
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("phastXXXXXX",CLEANUP=>1,TMPDIR=>1);
  $$settings{flanking}||=0;

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

  my $tempdir=tempdir("$$settings{tempdir}/phastXXXXXX",CLEANUP=>1,PREFIX=>$$settings{tempdir});
  my $db="$FindBin::RealBin/../lib/phast/phast.faa";
  logmsg "Running blastx against $db";

  # longest gene in phast is 8573bp, and all regions produced should have 
  # at least that length just in case.
  my $regions=`makeRegions.pl $fasta --numcpus $$settings{numcpus} --numchunks $$settings{numcpus} --overlapby 8573`;
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

  my $i=0;
  for my $region(@regions){
    my($contig,$coordinates)=split(/:/,$region);
    my($start,$stop)=split(/\-/,$coordinates);
      $stop||=$start;

    # The file name will have the start coordinate as its name
    my $file="$tempdir/$start.fna";
    open(SEQOUT,">",$file) or die "ERROR: could not write seq to temp file $file: $!";
    print SEQOUT ">".$contig."\n".$seq{$contig}->subseq($start,$stop)."\n";
    close SEQOUT;
    $i++;
  }
  my $threadsPerBlast=int($$settings{numcpus}/$i);
  $threadsPerBlast=1 if($threadsPerBlast<1);

  # Better parallelization: one fasta entry per cpu.
  # Split the query into multiple files and then figure out
  # how many cpus per blast job we need.

  # Perform blast on these split files.
  logmsg "Created blast input query files under $tempdir/*.fna";
  my $xargsCommand=qq(ls $tempdir/*.fna | xargs -P $$settings{numcpus} -n 1 sh -c '
    offset=\$(basename \$0 .fna); # Find the coordinate offset from the filename
    blastx -query \$0 -db $db -evalue 0.05 -outfmt 6 -num_threads $threadsPerBlast |\\
    perl -lane \"
      \\\$F[6]+=\$offset; # query coordinate offset
      \\\$F[7]+=\$offset; # query coordinate offset
      print join(\\\"\\\t\\\",\@F); # lots of backslashes because of language-inception
    \" > \$0.bls
  ');
  #die $xargsCommand;
  system($xargsCommand);
  die "ERROR with blastx: $!" if $?;
  #my $allResults=`blastx -query '$fasta' -db $db -evalue 0.05 -outfmt 6 -num_threads $$settings{numcpus}`;
  my @allResults=`cat $tempdir/*.fna.bls`;
  die "ERROR with cat on $tempdir/*.fna.bls" if($?);
  warn "No results were returned by blastx" if(!@allResults);

  my $flanking=$$settings{flanking}; #bp
  logmsg "Parsing results with a soft flanking distance of $flanking";
  my(%range, %seenRange);
  for my $result(@allResults){
    $result=~s/^\s+|\s+$//g; # trim
    my ($contig,$hit,$identity,$length,$gaps,$mismatches,$qstart,$qend,$sstart,$send,$e,$score)=split /\t/, $result;
    next if($score < 50 || $length < 20);
    # Don't bother the Range object if this is a range we've already seen.
    # This will speed things up since Number::Range is slow, and 
    # it will happen because each genome region can hit multiple phages.
    next if($seenRange{$contig}{$qstart}{$qend}++);

    # Make sure there is a range object for this contig.
    $range{$contig}||=Number::Range->new;
    my $lo=min($qstart,$qend);
    my $hi=max($qstart,$qend);

    # Add some coordinates between close hits based on
    # the flanking distance. Start from high to low
    # flanking numbers so that the longest range possible
    # is caught.

    # Flanking for lo
    my $loSoftFlank=$lo;
    for(my $i=$flanking;$i>0;$i--){
      if($range{$contig}->inrange($lo-$i)){
        $loSoftFlank=$lo-$i;
        last;
      }
    }
    # Flanking for hi
    my $hiSoftFlank=$hi;
    for(my $i=$flanking;$i>0;$i--){
      if($range{$contig}->inrange($hi+$i)){
        $hiSoftFlank=$hi+$i;
        last;
      }
    }

    # Add these coordinates to ranges
    no warnings;
    $range{$contig}->addrange($loSoftFlank..$hiSoftFlank);
  }

  logmsg "Finished adding coordinate ranges. Now creating a bed format";

  # Translate the ranges found in the Range objects into 
  # an array of [contig,start,stop]
  my @range;
  while(my($contig,$rangeObj)=each(%range)){
    my @rangeList=$rangeObj->rangeList;
    for(@rangeList){
      push(@range,[$contig,@$_]);
    }
  }

  return \@range;
}

sub usage{
  "Finds phages in a fasta file using phast
  Usage: $0 file.fasta
  --numcpus  1
  --tempdir  /tmp
  --flanking 0    Give 'soft' edges to ranges. If blast hits are this many
                  nt away from another blast hit, then join the ranges and
                  include any intermediate positions. If ranges cannot be
                  joined, then do not extend them by this flanking length.
  "
}
