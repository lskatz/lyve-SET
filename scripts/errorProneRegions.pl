#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw/basename/;
use Data::Dumper;
use Bio::SeqIO;
use Bio::FeatureIO;
use Bio::SeqFeature::Annotated;
use Bio::Annotation::SimpleValue;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg/;

local $0=basename($0);

# List all the error regexes as specifically as possible.
# The value of each regex is an array:
#   1. A score 1 to 100 representing how sure we are this is an error prone region
#   2. type - snp or indel
#   3. read - R1 or R2.  Use S for single-end
#
# References:
#   https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0976-y

# TODO: would it be faster if I put the /g into these regex declarations?
my %regionRegex=(
  GAII        => {
                   qr/(CGG)/i      => [qw(90 snp R1)],
                   qr/(GGG)/i      => [qw(90 snp R1R2)],
                   qr/(GCG)/i      => [qw(90 snp R1)],
                   qr/(CGG)/i      => [qw(90 snp R2)],
                   qr/(CCG)/i      => [qw(90 snp R2)],
                 },
  MiSeq       => {
                   qr/(GGG)/i      => [qw(90 snp R1R2)],
                   qr/(CGG)/i      => [qw(90 snp R1R2)],
                   qr/(TGG)/i      => [qw(90 snp R1)],
                   qr/(GCG)/i      => [qw(90 snp R2)],
                 },
  HiSeq       => {
                   qr/(GGG)/i      => [qw(90 snp R1R2)],
                   qr/(CGG)/i      => [qw(90 snp R1R2)],
                   qr/(AGG)/i      => [qw(90 snp R1R2)],
                 },
  iontorrent  => {
                   qr/(A{5,})/i    => [qw(90 indel S)],
                   qr/(C{5,})/i    => [qw(90 indel S)],
                   qr/(G{5,})/i    => [qw(90 indel S)],
                   qr/(T{5,})/i    => [qw(90 indel S)],
                   # Give homopolymers of 4 moderate importance
                   qr/(A{4})/i    => [qw(50 indel S)],
                   qr/(C{4})/i    => [qw(50 indel S)],
                   qr/(G{4})/i    => [qw(50 indel S)],
                   qr/(T{4})/i    => [qw(50 indel S)],
                   # homopolymers of three might have errors but are less important
                   qr/(A{3})/i     => [qw(10 indel S)],
                   qr/(C{3})/i     => [qw(10 indel S)],
                   qr/(G{3})/i     => [qw(10 indel S)],
                   qr/(T{3})/i     => [qw(10 indel S)],
                 },
);

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(platform=s format=s min-score=i help)) or die $!;
  die usage($settings) if($$settings{help} || !@ARGV);

  $$settings{'min-score'}||=1;
  $$settings{format}||="bed";
  # The platform is actually an array
  $$settings{platform}||=join(",",keys(%regionRegex));
  $$settings{platform}=[split(/,/,$$settings{platform})];

  my @infile=@ARGV;

  my $seq    =readSeq(\@infile,$settings);
  my $region =errorProneRegions($seq,$settings);

  # Write to stdout
  my $featout=Bio::FeatureIO->new(-format=>$$settings{format}, -version=>3);
  for my $feat(@$region){
    $featout->write_feature($feat);
  }

  return 0;
}

sub readSeq{
  my($infile,$settings)=@_;
  my %seq;
  
  for(@$infile){
    my $in=Bio::SeqIO->new(-file=>$_);
    while(my $seq=$in->next_seq){
      if($seq{$seq->id}){
        die "ERROR: found the same contig twice: ".$seq->id;
      }
      $seq{$seq->id}=$seq->seq;
    }
  }
  return \%seq;
}

# Search for error prone regions in each sequence
sub errorProneRegions{
  my($seqHash,$settings)=@_;
  my @region;
  while(my($seqname,$seq)=each(%$seqHash)){
    my $pos=searchForPositions($seqname,$seq,$settings);
    push(@region, @$pos);
  }
  return \@region;
}

# Search for error prone regions within a sequence
sub searchForPositions{
  my($seqname,$seq,$settings)=@_;
  logmsg $seqname;

  my $source_tag=basename($0);
  
  my @feat;
  for my $platform(@{ $$settings{platform} }){ 
    logmsg $platform;
    for my $regex(keys(%{ $regionRegex{$platform} })){
      my $regexValue=$regionRegex{$platform}{$regex};
      next if($$regexValue[0] < $$settings{'min-score'});
      while($seq=~/$regex/g){
        my($start,$end)=($-[0]+1,$+[0]); # 1-based
        my $match = $1;
        my $feat=Bio::SeqFeature::Annotated->new(
            -seq_id  => $seqname,
            -start   => $start,
            -end     => $end,
            -primary => "errorProne_$platform",
            -score   => $$regexValue[0],
            -source  => $source_tag,
            # http://song.cvs.sourceforge.net/viewvc/song/ontology/sofa.obo?revision=1.217
            -type    =>"possible_base_call_error",
        );
        $feat->add_Annotation(Bio::Annotation::SimpleValue->new(
          -tagname=>"Name",
          -value  => sprintf("Platform:%s;match:%s,type:%s,read:%s",$platform,$match,$$regexValue[1],$$regexValue[2]),
        ));
        push(@feat,$feat);
      }
    }
  }

  # Sort by contig name, or by
  # start position, or by
  # the score
  return [sort { 
    $a->seq_id cmp $b->seq_id || 
    $a->start <=> $b->start || 
    $b->score <=> $a->score
  } @feat];
}


sub usage{
  my($settings)=@_;
  my @platform=keys(%regionRegex);
  my $platforms=join(", ",@platform);
  my $platformEx=join(",",@platform[0,1]);
  "$0: finds error-prone regions in a sequence file
  Usage: $0 assembly.fasta > out.bed
  --platform   ''    Which platform to define errors.
                     Values: $platforms
                     Can define multiple platforms by using
                     commas, e.g., --platform $platformEx
                     Use no value to include all regions
  --min-score   1    Minimum score to accept (1-100)
  --format    bed    Output format.  E.g., bed, gff
  "
}

