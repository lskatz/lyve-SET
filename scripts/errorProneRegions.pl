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
# Capturing parentheses:
#   $1 - the motif
#   $2 - the error-prone base
#
# References:
#   https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0976-y

# TODO: would it be faster if I put the /g into these regex declarations?
my %regionRegex=(
  GAII        => {
                   qr/(CG(G))/i      => [qw(90 snp R1)],
                   qr/(GG(G))/i      => [qw(90 snp R1R2)],
                   qr/(GC(G))/i      => [qw(90 snp R1)],
                   qr/(CG(G))/i      => [qw(90 snp R2)],
                   qr/(CC(G))/i      => [qw(90 snp R2)],
                 },
  MiSeq       => {
                   qr/(GG(G))/i      => [qw(90 snp R1R2)],
                   qr/(CG(G))/i      => [qw(90 snp R1R2)],
                   qr/(TG(G))/i      => [qw(90 snp R1)],
                   qr/(GC(G))/i      => [qw(90 snp R2)],
                 },
  HiSeq       => {
                   qr/(GG(G))/i      => [qw(90 snp R1R2)],
                   qr/(CG(G))/i      => [qw(90 snp R1R2)],
                   qr/(AG(G))/i      => [qw(90 snp R1R2)],
                 },
  iontorrent  => {
                   # Homopolymers of 5+ are important
                   qr/(A{4,}(A))/i    => [qw(90 indel S)],
                   qr/(C{4,}(C))/i    => [qw(90 indel S)],
                   qr/(G{4,}(G))/i    => [qw(90 indel S)],
                   qr/(T{4,}(T))/i    => [qw(90 indel S)],
                   # Give homopolymers of 4 moderate importance
#                  qr/(A{3}(A))/i    => [qw(50 indel S)],
#                  qr/(C{3}(C))/i    => [qw(50 indel S)],
#                  qr/(G{3}(G))/i    => [qw(50 indel S)],
#                  qr/(T{3}(T))/i    => [qw(50 indel S)],
                   # homopolymers of three might have errors but are less important
#                  qr/(A{2}(A))/i     => [qw(10 indel S)],
#                  qr/(C{2}(C))/i     => [qw(10 indel S)],
#                  qr/(G{2}(G))/i     => [qw(10 indel S)],
#                  qr/(T{2}(T))/i     => [qw(10 indel S)],
                 },
);

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(platform=s format=s min-score=i help motifs)) or die $!;
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
        my($motifStart,$motifEnd)=($-[1]+1,$+[1]); # 1-based
        my($errorStart,$errorEnd)=($-[2]+1,$+[2]); # 1-based
        my $motif = $1;
        my $match = $2;

        # Add the match
        my $errorFeat=Bio::SeqFeature::Annotated->new(
            -seq_id  => $seqname,
            -start   => $errorStart,
            -end     => $errorEnd,
            -primary => "errorProne_$platform",
            -score   => $$regexValue[0],
            -source  => $source_tag,
            # http://song.cvs.sourceforge.net/viewvc/song/ontology/sofa.obo?revision=1.217
            -type    =>"possible_base_call_error",
        );
        $errorFeat->add_Annotation(Bio::Annotation::SimpleValue->new(
          -tagname=>"Name",
          -value  => sprintf("error-prone:%s",$match),
        ));
        push(@feat,$errorFeat);

        # Add the motif if requested
        if($$settings{motifs}){
          my $motifFeat=Bio::SeqFeature::Annotated->new(
              -seq_id  => $seqname,
              -start   => $motifStart,
              -end     => $motifEnd,
              -primary => "errorProne_$platform",
              -score   => $$regexValue[0],
              -source  => $source_tag,
              # http://song.cvs.sourceforge.net/viewvc/song/ontology/sofa.obo?revision=1.217
              -type    =>"nucleotide_motif",
          );
          $motifFeat->add_Annotation(Bio::Annotation::SimpleValue->new(
            -tagname=>"Name",
            -value  => sprintf("Platform:%s;motif:%s,type:%s,read:%s",$platform,$motif,$$regexValue[1],$$regexValue[2]),
          ));
          push(@feat,$motifFeat);
        }
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
  --motifs           Output motif regions in addition to
                     error-prone sites
  "
}

