#!/usr/bin/env perl

#Author: Lee Katz, Eishita Tyagi
##########################################

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/fileparse basename dirname/;
use File::Spec;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg @vcfExt/;
#use lib "$FindBin::RealBin/../lib/lib/perl5";
#use File::Slurp qw/write_file/;

# Figure out the script name and directory
my ($name,$scriptsdir,$suffix)=fileparse($0);
$scriptsdir=File::Spec->rel2abs($scriptsdir);
my $programdir = ""; 
if ($scriptsdir =~ /(.*)\/scripts/) {$programdir = $1;}

# Check for the conf file
my $conf="$programdir/config/snpEff.conf";
if(!-e $conf){
  logmsg "WARNING: Could not find $conf. Attempting to install it.";
  system("cd $programdir && make install-snpEff");
  die "ERROR: could not find and could not install $conf" if $?;
}

# Check for the jar file
my $JAR = "$scriptsdir/../lib/snpEff.jar";
if(!-e $JAR){
  logmsg "WARNING: Could not find $JAR.  Attempting to install it.";
  system("cd $scriptsdir/.. && make install-snpEff");
  die "ERROR: could not find and could not install SnpEff" if $?;
}

exit(main());

sub main{
  #Get run paramaters
  my $settings={numcpus=>1,force=>0,sumspecies=>0,deleteTmp=>0};
  GetOptions($settings,qw(outdir=s reference=s genbank=s organism=s codon=s force! deleteTmp! numcpus=i help));
   die usage()."\n" if($$settings{help} || !@ARGV);  

  # Some defaults
  $$settings{outdir} ||= basename($0).".out";
    $$settings{outdir}=File::Spec->rel2abs($$settings{outdir});
  $$settings{datadir} ||= "$$settings{outdir}/reference"; 
  $$settings{organism} ||=  "reference";
  $$settings{codon} ||=  "Standard";  

  my @VCF=@ARGV;


  #check database build params
  if ((not defined $$settings{reference}) && (not defined $$settings{genbank})){ 
    die "\nMust specify reference name for snpEff to search for OR provide genbank file of desired reference\n".usage()."\n";
  }
  #checking if name provided was correct
  if (defined $$settings{reference}) {
    if ($$settings{reference} =~ /(.+).genome/){
      $$settings{reference} = $1;
    } 
    my $refcheck = `grep -w "$$settings{reference}.genome" $conf`; 
    die "ERROR: was not able find $$settings{reference} in $conf" if $?;
    if ($refcheck !~ /.genome/){
      die "\nReference name does not exist in $conf. Please check or provide --genbank\n".usage()."\n";
    }
  }

  # Or checking if gbk provided
  if ((defined $$settings{organism}) && (not defined $$settings{genbank})){
    die "\nProvide genbank file of desired reference\n".usage()."\n";
  }  
  if (defined $$settings{organism}){
    $$settings{organism} =~ s/ +/_/;
  } 
  if (!$$settings{codon}) {
    #checking if codon table type provided was correct
    if ($$settings{codon} =~ /codon.(.+)/){
      $$settings{reference} = $1;
    } 
    my $refcheck = `grep -w "codon.$$settings{codon}" $conf`; 
    die "ERROR: could not look for codon.$$settings{codon} in $conf" if $?;
    if ($refcheck !~ /codon./){
      die "\nCodon table type does not exist in snpEff config. Please check or modify snpEff.config\n".usage()."\n"
    }
  }
      
  #creating required run directories
  system("mkdir -pv $$settings{outdir} >&2");
  die "ERROR: could not make directory $$settings{outdir}" if $?;
  system("mkdir -pv $$settings{datadir} >&2");
  die "ERROR: could not make directory $$settings{datadir}" if $?;
  
  #starting program
  #Step 1: build snpeff database 
  BuildSnpEffDB($settings);

  #Step 2: Run snpeff against database 
  RunSnpEff(\@VCF,$settings);
  `mv $$settings{outdir}/snpEff.config $$settings{datadir}/.`;

  if ($$settings{deleteTmp} == 1){`rm -rf $$settings{datadir}`;}
  return 0;
}

##########################################################################################################
sub BuildSnpEffDB{
  my ($settings) = @_;
  my $configFile="$$settings{outdir}/snpEff.config";
  system("cp -v $conf $configFile");
  die "ERROR: could not copy $conf to $configFile" if $?;

  system("sed -i 's/data.dir.*/data.dir = reference/' $configFile");
  die "ERROR: could not run sed on $configFile" if $?;

  # Download the reference database if it exists and if it was defined
  if (defined $$settings{reference}) {
    system("java -jar $JAR download -c $configFile -v $$settings{reference}");
    die "ERROR: could not run java -jar $JAR" if $?;
  }

  if (defined $$settings{genbank}) {
    # Set up the SnpEff database directory
    system("mkdir -pv $$settings{datadir}/$$settings{organism} >&2");
    die "ERROR: could not make folder $$settings{datadir}/$$settings{organism}" if $?;
    system("cp -v $$settings{genbank} $$settings{datadir}/$$settings{organism}/genes.gbk >&2");
    die "ERROR could not copy $$settings{genbank} to $$settings{datadir}/$$settings{organism}/genes.gbk" if $?;

    # Write the new config line
    open(SNPEFFCONF,">>","$configFile") or die "ERROR: could not write to file $configFile: $!";
    print SNPEFFCONF "\n# $$settings{organism}, version 1\n$$settings{organism}.genome : $$settings{organism}\n";

    # Find what seqnames we have in the genbank file
    my $seqin=Bio::SeqIO->new(-file=>"$$settings{datadir}/$$settings{organism}/genes.gbk",-verbose=>-1);
    my @chr=();
    while(my $seq=$seqin->next_seq){
      push(@chr,$seq->id);
    }
    my $chrStr=join(" ",@chr);

    # List the seqnames in the index line in the config file
    print SNPEFFCONF "\t$$settings{organism}.chromosomes : $chrStr\n";
    # List each seqname entry with the correct codon table
    for (my $i=0; $i < @chr; $i++) {
      print SNPEFFCONF "\t$$settings{organism}.$chr[$i].codonTable : $$settings{codon}\n";
    }
    close SNPEFFCONF;

    # Create the database using the new config file.
    system("java -jar $JAR build -c $configFile -v -genbank $$settings{organism} >&2");
    die "ERROR building the SnpEff DB" if $?;
  }
}
##########################################################################################################

##########################################################################################################
sub RunSnpEff{
  my ($VCF,$settings) = @_; 
  my $name = "";
  my @list = @$VCF;
  for my $vcf(@$VCF){
    my $b=basename($vcf,@vcfExt);
    my $outfile="$$settings{outdir}/$b.vcf";
    my $command="java -jar $JAR -stats $$settings{outdir}/snpEff_summary.html -v -c $$settings{outdir}/snpEff.config $$settings{organism} $vcf > $outfile";
    system($command);
    die "ERROR: running SnpEff\n  $command" if $?;

    # tabix, bgzip
    system("bgzip -f $outfile >&2 && tabix -f $outfile.gz >&2");
    die "ERROR with bgzip and/or tabix" if $?;
  }
}
##########################################################################################################

sub usage{ 
  local $0=basename($0);
  "$0: Runs SNP annotation on a VCF file
  Usage: $0 --genbank in.gbk --outdir dir *.vcf[.gz]
  
      --genbank      Full reference genbank file with sequence. Cannot be used with --reference
      --outdir       output directory
      --organism     name of organism you'd like entered in your snpEff config. Else will call it 'reference'
      --codon        codon table type as per snpEff.config settings. Default 'Standard'
      --numcpus      1 (default)
      --deleteTmp    (default do not delete temporary files)  
  "
      #--reference    name of reference for snpEff to build from SnpEff.config
      #--force (default: don't overwrite output directory if exists. Currently this option will only overwrite files and does not create a fresh environment)
}

