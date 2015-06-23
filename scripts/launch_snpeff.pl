#!/usr/bin/perl 
#Author: Lee Katz, Eishita Tyagi #Revised: 06/15/2015
#########################################################################################################################

use FindBin;
use lib "$FindBin::RealBin/../lib";
use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use File::Spec;

sub logmsg {local $0=basename $0;my $FH = *STDOUT; print $FH "\n$0: ".(caller(1))[3].": @_\n";}
my ($name,$scriptsdir,$suffix)=fileparse($0);
$scriptsdir=File::Spec->rel2abs($scriptsdir);my $programdir = ""; if ($scriptsdir =~ /(.*)\/scripts/) {$programdir = $1;}
my $JAR = "$programdir/lib/snpEff.jar";  if (! -e $JAR) { die "ERROR: jar file does not exist at $JAR\nPlease edit $0 to reflect where the jar file is.\n";}

exit(main());

sub main{
	#Get run paramaters
	my $settings={numcpus=>8,threads=>5,force=>0,sumspecies=>0,deleteTmp=>0};
	GetOptions($settings,qw(indir=s outdir=s reference=s genbank=s organism=s codon=s force! deleteTmp! threads=s numcpus=s help));
 	die usage()."\n" if($$settings{help});	

	#check params
	$$settings{indir} || die "\nMissing parameter: -indir\n".usage()."\n";
  	$$settings{outdir} ||= $$settings{indir}."_out";

	#check database build params
	if ((not defined $$settings{reference}) && (not defined $$settings{genbank})){ 
		die "\nMust specify reference name for snpEff to search for OR provide genbank file of desired reference\n".usage()."\n";}
	if (defined $$settings{reference}) {my $refcheck = ""; if ($$settings{reference} =~ /(.+).genome/){$$settings{reference} = $1;} #checking if name provided was correct
		$refcheck = `grep -w "$$settings{reference}.genome" $programdir/config/original/snpEff.conf`; 
		if ($refcheck !~ /.genome/){die "\nReference name does not exist in $programdir/config/original/snpEff.conf. Please check or provide --genbank\n".usage()."\n"}}
	if ((defined $$settings{organism}) && (not defined $$settings{genbank})) #or checking if gbk provided 
		{die "\nProvide genbank file of desired reference\n".usage()."\n";}	
	if (defined $$settings{organism}){$$settings{organism} =~ s/ +/_/;} $$settings{organism} ||=  "reference";
	$$settings{codon} ||=  "Standard";	
	if (defined $$settings{codon}) {my $codoncheck = ""; if ($$settings{codon} =~ /codon.(.+)/){$$settings{reference} = $1;} #checking if codon table type provided was correct
		my $refcheck = `grep -w "codon.$$settings{codon}" $programdir/config/original/snpEff.conf`; 
		if ($refcheck !~ /codon./){die "\nCodon table type does not exist in snpEff config. Please check or modify snpEff.config\n".usage()."\n"}}
			
	#creating required run directories
	if ((-d $$settings{outdir}) && ($$settings{force} == 1)){`rm -rf $$settings{outdir}`; `mkdir -p $$settings{outdir}`;}
	if ((-d $$settings{outdir}) && ($$settings{force} == 0)){die "\noutput directory exists. Specify new name or --force\n".usage()."\n";}
	if (! -d $$settings{outdir}){`mkdir -p $$settings{outdir}`;} $$settings{datadir} ||= "$$settings{outdir}/reference"; `mkdir -p $$settings{datadir}`;
	
	for my $param (qw(indir outdir)){ my $b=$param; $b=~s/dir$//; $$settings{$param}||=$b;
  		die "ERROR: Could not find $param under $$settings{$param}/ \n".usage() if(!-d $$settings{$param});
		$$settings{$param}=File::Spec->rel2abs($$settings{$param});}

	#starting program
	#Step 1: build snpeff database 
	BuildSnpEffDB($settings);

	#Step 1: Run snpeff against database 
	RunSnpEff($settings);
	`mv $$settings{outdir}/snpEff.config $$settings{datadir}/.`;

	if ($$settings{deleteTmp} == 1){`rm -rf $$settings{datadir}`;}
	return 0;}

##########################################################################################################
sub usage{ #subroutine to print usage of program
	"
 	Usage: launch_snpEff.pl --indir=<Input directory with .vcf files (one directory at a time)> --reference=<name of reference for snpEff to build from SnpEff.config.> 
                [ --genbank=<Full reference genbank file (genes.gbk) WITH SEQUENCE, if building custom snpEff reference database. Cannot be used with --reference> 
		  --organism=<name of organism you'd like entered in your snpEff config. Else will call it 'reference'>	
		  --codon=<codon table type as per snpEff.config settings. Default 'Standard'>
		  --outdir=<indir>_out
                  --numcpus=8 (default)
		  --threads=3 (default)
		  --force (default: don't overwrite output directory if exists. Currently this option will only overwrite files and does not create a fresh environment)
		  --deleteTmp (default do not delete temporary files)  ]	
	"}
##########################################################################################################

##########################################################################################################
sub BuildSnpEffDB{
	my ($settings) = @_;
	`cp $programdir/config/original/snpEff.conf $$settings{outdir}/snpEff.config`;
	`sed -i 's/data.dir.*/data.dir = reference/' $$settings{outdir}/snpEff.config`;
	if (defined $$settings{reference}) {`java -jar $JAR download -c $$settings{outdir}/snpEff.config -v $$settings{reference}`;}
	if (defined $$settings{genbank}) {
		`mkdir -p $$settings{datadir}/$$settings{organism}`;`cp $$settings{genbank} $$settings{datadir}/$$settings{organism}/genes.gbk`;
		`echo "\n# $$settings{organism}, version 1\n$$settings{organism}.genome : $$settings{organism}" >> $$settings{outdir}/snpEff.config`;
		my $chr = `grep "LOCUS" $$settings{datadir}/$$settings{organism}/genes.gbk`;
		my @chr = split ("\n",$chr); my @temp = ();
		for (my $i=0; $i < @chr; $i++) {if ($chr[$i] =~ /^LOCUS\s+(.+)\s+\d+\s+bp.*/){$chr[$i] = $1; $chr[$i] =~ s/\s+//g; push (@temp,"$chr[$i],");}}
		$temp[$#temp] =~ s/,//g; `echo "\t$$settings{organism}.chromosomes : @temp" >> $$settings{outdir}/snpEff.config`;
		for (my $i=0; $i < @chr; $i++) {`echo "\t$$settings{organism}.$chr[$i].codonTable : $$settings{codon}" >> $$settings{outdir}/snpEff.config`;}
		`java -jar $JAR build -c $$settings{outdir}/snpEff.config -v -genbank $$settings{organism}`;}
	}
##########################################################################################################

##########################################################################################################
sub RunSnpEff{
	my ($settings) = @_; my $name = "";
	my @list = <$$settings{indir}/*.vcf>;
	for (my $i = 0; $i <= @list-1; $i = $i + 1){
		if ($list[$i] =~ /$$settings{indir}\/(.*)\.vcf/){$name = $1} 
		`java -jar $JAR -stats $$settings{outdir}/snpEff_summary.html -v -c $$settings{outdir}/snpEff.config $$settings{organism} $$settings{indir}/$name.vcf > $$settings{outdir}/$name.snpEFF.vcf`;}
	}
##########################################################################################################





















