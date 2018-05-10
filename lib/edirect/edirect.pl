#!/usr/bin/env perl

# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#            National Center for Biotechnology Information (NCBI)
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government do not place any restriction on its use or reproduction.
#  We would, however, appreciate having the NCBI and the author cited in
#  any work or product based on this material.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
# ===========================================================================
#
# File Name:  edirect.pl
#
# Author:  Jonathan Kans
#
# Version Creation Date:   8/20/12
#
# ==========================================================================

# Entrez Direct - EDirect

# use strict;
use warnings;

my ($LibDir, $ScriptName);

use File::Spec;

# EDirect version number

$version = "8.60";

BEGIN
{
  my $Volume;
  ($Volume, $LibDir, $ScriptName) = File::Spec->splitpath($0);
  $LibDir = File::Spec->catpath($Volume, $LibDir, '');
  if (my $RealPathname = eval {readlink $0}) {
    do {
      $RealPathname = File::Spec->rel2abs($RealPathname, $LibDir);
      ($Volume, $LibDir, undef) = File::Spec->splitpath($RealPathname);
      $LibDir = File::Spec->catpath($Volume, $LibDir, '')
    } while ($RealPathname = eval {readlink $RealPathname});
  } else {
    $LibDir = File::Spec->rel2abs($LibDir)
  }
  $LibDir .= '/aux/lib/perl5';
}
use lib $LibDir;

# usage - edirect.pl -function arguments

use Data::Dumper;
use Encode;
use Getopt::Long;
use HTML::Entities;
use LWP::Simple;
use LWP::UserAgent;
use Net::hostent;
use POSIX;
use Time::HiRes;
use URI::Escape;

# required first argument is name of function to run

$fnc = shift or die "Must supply function name on command line\n";

# get starting time

$begin_time = Time::HiRes::time();

# definitions

use constant false => 0;
use constant true  => 1;

# URL address components

$base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/";

$ecitmat  = "ecitmatch.cgi";
$efetch   = "efetch.fcgi";
$einfo    = "einfo.fcgi";
$elink    = "elink.fcgi";
$epost    = "epost.fcgi";
$esearch  = "esearch.fcgi";
$espell   = "espell.fcgi";
$esummary = "esummary.fcgi";

# utility subroutines

sub clearflags {
  @rest = ();
  %labels = ();
  %macros = ();
  $alias = "";
  $author = "";
  $basx = "";
  $batch = false;
  $chr_start = -1;
  $chr_stop = -1;
  $clean = false;
  $cmd = "";
  $complexity = 0;
  $db = "";
  $dbase = "";
  $dbs = "";
  $dbto = "";
  $debug = false;
  $drop = false;
  $dttype = "";
  $emaddr = "";
  $email = "";
  $err = "";
  $extend = -1;
  $extrafeat = -1;
  $field = "";
  $fields = false;
  $feature = "";
  $gtype = "";
  $help = false;
  $holding = "";
  $http = "";
  $id = "";
  $input = "";
  $journal = "";
  $just_num = false;
  $key = "";
  $lbl = "";
  $links = false;
  $location = "";
  $log = false;
  $max = 0;
  $meadow = "";
  $min = 0;
  $mndate = "";
  $mode = "";
  $molecule = "";
  $mxdate = "";
  $name = "";
  $neighbor = false;
  $num = "";
  $organism = "";
  $output = "";
  $page = "";
  $pair = "";
  $pipe = false;
  $pub = "";
  $query = "";
  $raw = false;
  $related = false;
  $rldate = 0;
  $seq_start = 0;
  $seq_stop = 0;
  $showgi = false;
  $silent = false;
  $sort = "";
  $source = "";
  $spell = false;
  $split = "";
  $status = "";
  $stp = "";
  $strand = "";
  $style = "";
  $tool = "";
  $trim = false;
  $trunc = false;
  $tuul = "";
  $type = "";
  $verbose = false;
  $volume = "";
  $web = "";
  $word = false;
  $year = "";

  $stop_words="#a#about#again#all#almost#also#although#always#among#an#and#" .
  "another#any#are#as#at#be#because#been#before#being#between#both#but#by#can#" .
  "could#did#do#does#done#due#during#each#either#enough#especially#etc#for#" .
  "found#from#further#had#has#have#having#here#how#however#i#if#in#into#is#it#" .
  "its#itself#just#kg#km#made#mainly#make#may#mg#might#ml#mm#most#mostly#" .
  "must#nearly#neither#no#nor#obtained#of#often#on#our#overall#perhaps#pmid#" .
  "quite#rather#really#regarding#seem#seen#several#should#show#showed#shown#" .
  "shows#significantly#since#so#some#such#than#that#the#their#theirs#them#" .
  "then#there#therefore#these#they#this#those#through#thus#to#upon#use#used#" .
  "using#various#very#was#we#were#what#when#which#while#with#within#without#would#";

  $os = "$^O";

  $api_key = "";
  $api_key = $ENV{NCBI_API_KEY} if defined $ENV{NCBI_API_KEY};
}

sub do_sleep {

  if ( $api_key ne "" ) {
    if ( $log ) {
      print STDERR "sleeping 1/10 second\n";
    }
    Time::HiRes::usleep(110);
    return;
  }

  if ( $log ) {
    print STDERR "sleeping 1/3 second\n";
  }
  Time::HiRes::usleep(350);
}

# gets a live UID for any database

sub get_zero_uid {

  my $db = shift (@_);

  my $val = "";

  %zeroUidHash = (
    'assembly'         =>  '443538',
    'bioproject'       =>  '146229',
    'biosample'        =>  '3737421',
    'biosystems'       =>  '1223165',
    'blastdbinfo'      =>  '1023214',
    'books'            =>  '1371014',
    'cdd'              =>  '277499',
    'clinvar'          =>  '10510',
    'clone'            =>  '18646800',
    'dbvar'            =>  '6173073',
    'gap'              =>  '872875',
    'gds'              =>  '200022309',
    'gencoll'          =>  '398148',
    'gene'             =>  '3667',
    'genome'           =>  '52',
    'genomeprj'        =>  '72363',
    'geoprofiles'      =>  '16029743',
    'gtr'              =>  '558757',
    'homologene'       =>  '510',
    'medgen'           =>  '162753',
    'mesh'             =>  '68007328',
    'nlmcatalog'       =>  '0404511',
    'nuccore'          =>  '1322283',
    'nucest'           =>  '338968657',
    'nucgss'           =>  '168803471',
    'nucleotide'       =>  '1322283',
    'omim'             =>  '176730',
    'pcassay'          =>  '1901',
    'pccompound'       =>  '16132302',
    'pcsubstance'      =>  '126522451',
    'pmc'              =>  '209839',
    'popset'           =>  '27228303',
    'probe'            =>  '9997691',
    'protein'          =>  '4557671',
    'proteinclusters'  =>  '2945638',
    'pubmed'           =>  '2539356',
    'pubmedhealth'     =>  '27364',
    'seqannot'         =>  '9561',
    'snp'              =>  '137853337',
    'sra'              =>  '190091',
    'structure'        =>  '61024',
    'taxonomy'         =>  '562',
    'unigene'          =>  '1132160',
  );

  if ( defined $zeroUidHash{$db} ) {
    $val = $zeroUidHash{$db};
  }

# return $val below to restore debugging test

  return 0;
}

# support for substitution of (#keyword) to full query phrase or URL component

sub map_labels {

  my $qury = shift (@_);

  if ( $query !~ /\(#/ ) {
    return $qury;
  }

  if ( scalar (keys %labels) > 0 ) {
    for ( keys %labels ) {
      $ky = $_;
      $vl = $labels{$_};
      $qury =~ s/\((#$ky)\)/\($vl#\)/g;
    }
    $qury =~ s/\((\w+)#\)/\(#$1\)/g;
  }

  return $qury;
}

sub map_macros {

  my $qury = shift (@_);

  if ( $qury !~ /\(#/ ) {
    return $qury;
  }

  if ( scalar (keys %macros) > 0 ) {
    for ( keys %macros ) {
      $ky = $_;
      $vl = $macros{$_};
      $qury =~ s/\((\#$ky)\)/$vl/g;
    }
  }

  return $qury;
}

sub read_aliases {

  if ( $alias ne "" ) {
    if (open (my $PROXY_IN, $alias)) {
      while ( $thisline = <$PROXY_IN> ) {
        $thisline =~ s/\r//;
        $thisline =~ s/\n//;
        $thisline =~ s/ +/ /g;
        $thisline =~ s/> </></g;

        if ( $thisline =~ /(.+)\t(.+)/ ) {
          $ky = $1;
          $vl = $2;
          $vl =~ s/\"//g;
          $macros{"$ky"} = "$vl";
        }
      }
      close ($PROXY_IN);
    } else {
      print STDERR "Unable to open alias file '$alias'\n";
    }
  }
}

# base EUtils URL can be overridden for access to test versions of server

sub adjust_base {

  if ( $basx eq "" ) {

    # if base not overridden, check URL of previous query, stick with main or colo site,
    # since history server data for EUtils does not copy between locations, by design

    if ( $web ne "" and $web =~ /NCID_\d+_\d+_(\d+)\.\d+\.\d+\.\d+_\d+_\d+_\d+/ ) {
      if ( $1 == 130 ) {
        $base = "https://eutils.be-md.ncbi.nlm.nih.gov/entrez/eutils/";
      } elsif ( $1 == 165 ) {
        $base = "https://eutils.st-va.ncbi.nlm.nih.gov/entrez/eutils/";
      }
    }
    return;
  }

  # shortcut for eutilstest base
  if ( $basx eq "test" ) {
    $basx = "https://eutilstest.ncbi.nlm.nih.gov/entrez/eutils";
  }

  if ( $basx =~ /\(#/ ) {
    $basx = map_macros ($basx);
  }

  if ( $basx !~ /^https:\/\// ) {
    $basx = "https://" . $basx;
  }

  if ( $basx !~ /\/$/ ) {
    $basx .= "/";
  }

  if ( $basx !~ /\/entrez\/eutils\/$/ ) {
    $basx .= "entrez/eutils/";
  }

  $base = $basx;
}

# ensure that ENTREZ_DIRECT data structure contains required fields

sub test_edirect {

  my $dbsx = shift (@_);
  my $webx = shift (@_);
  my $keyx = shift (@_);
  my $numx = shift (@_);
  my $label = shift (@_);

  if ( $dbsx eq "" ) {
    close (STDOUT);
    die "Db value not found in $label input\n";
  }
  if ( $webx eq "" ) {
    close (STDOUT);
    die "WebEnv value not found in $label input\n";
  }
  if ( $keyx eq "" ) {
    close (STDOUT);
    die "QueryKey value not found in $label input\n";
  }
  if ( $numx eq "" ) {
    close (STDOUT);
    die "Count value not found in $label input\n";
  }
}

sub get_email {

  # adapted from code provided by Aaron Ucko

  my $addr = "";
  if (defined $ENV{EMAIL}) {
    $addr = $ENV{EMAIL};
  } else {
    # Failing that, try to combine the username from USER or whoami
    # with the contents of /etc/mailname if available or the system's
    # qualified host name.  (Its containing domain may be a better
    # choice in many cases, but anyone contacting abusers can extract
    # it if necessary.)
    my $lhs = $ENV{USER} || `whoami`;
    my $rhs = "";
    if (-r '/etc/mailname') {
      $rhs = `cat /etc/mailname`;
    } else {
      my @uname = POSIX::uname();
      $rhs = $uname[1];
      if ($rhs !~ /\./) {
        # clearly unqualified, try to resolve back and forth
        my $h = gethostbyname($rhs);
        if (defined $h  &&  $h->name =~ /\./) {
          $rhs = $h->name;
        }
      }
    }
    $addr = $lhs . '@' . $rhs;
  }

  return $addr;
}

# elink and epost currently need a separate ESearch to get the correct result count

sub get_count {

  my $dbsx = shift (@_);
  my $webx = shift (@_);
  my $keyx = shift (@_);
  my $tulx = shift (@_);
  my $emlx = shift (@_);

  test_edirect ( $dbsx, $webx, $keyx, "1", "count" );

  $url = $base . $esearch;
  $url .= "?db=$dbsx&query_key=$keyx&WebEnv=$webx";
  $url .= "&retmax=0&usehistory=y";

  $url .= "&edirect=$version";

  if ( $os ne "" ) {
    $url .= "&edirect_os=$os";
  }

  if ( $api_key ne "" ) {
    $url .= "&api_key=$api_key";
  }

  if ( $tulx eq "" ) {
    $tulx = "entrez-direct";
  }
  if ( $tulx ne "" ) {
    $tulx =~ s/\s+/\+/g;
    $url .= "&tool=$tulx-count";
  }

  if ( $emlx eq "" ) {
    $emlx = get_email ();
  }
  if ( $emlx ne "" ) {
    $emlx =~ s/\s+/\+/g;
    $url .= "&email=$emlx";
  }

  $keyx = "";
  $numx = "";
  $errx = "";

  $output = get ($url);

  if ( ! defined $output ) {
    print STDERR "Failure of '$url'\n";
    return "", "";
  }

  if ( $output eq "" ) {
    print STDERR "No get_count output returned from '$url'\n";
    return "", ""
  }

  if ( $debug ) {
    print STDERR "$output\n";
  }

  $keyx = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
  $numx = $1 if ($output =~ /<Count>(\S+)<\/Count>/);
  $errx = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);

  if ( $errx ne "" ) {
    close (STDOUT);
    die "ERROR in count output: $errx\nURL: $url\n\n";
  }

  if ( $numx eq "" ) {
    die "Count value not found in count output - WebEnv $webx\n";
  }

  return $numx, $keyx;
}

sub get_uids {

  my @working = ();

  my $dbsx = shift (@_);
  my $webx = shift (@_);
  my $keyx = shift (@_);
  my $sttx = shift (@_);
  my $chkx = shift (@_);
  my $numx = shift (@_);
  my $tulx = shift (@_);
  my $emlx = shift (@_);

  test_edirect ( $dbsx, $webx, $keyx, $numx, "uids" );

  # adjust retmax if -stop has been overridden
  # e.g., efetch -format uid -start 2 -stop 9

  if ( $sttx + $chkx > $numx ) {
    $chkx = $numx - $sttx;
  }

  $url = $base . $esearch . "?db=$dbsx&query_key=$keyx&WebEnv=$webx";
  $url .= "&rettype=uilist&retmode=text";
  $url .= "&retstart=$sttx&retmax=$chkx";

  $url .= "&edirect=$version";

  if ( $os ne "" ) {
    $url .= "&edirect_os=$os";
  }

  if ( $api_key ne "" ) {
    $url .= "&api_key=$api_key";
  }

  if ( $tulx eq "" ) {
    $tulx = "edirect";
  }
  if ( $tulx ne "" ) {
    $tulx =~ s/\s+/\+/g;
    $url .= "&tool=$tulx";
  }

  if ( $emlx eq "" ) {
    $emlx = get_email ();
  }
  if ( $emlx ne "" ) {
    $emlx =~ s/\s+/\+/g;
    $url .= "&email=$emlx";
  }

  if ( $debug ) {
    print STDERR "$url\n";
  }

  $keep_trying = true;
  for ( $try = 0; $try < 3 && $keep_trying; $try++) {

    $data = get ($url);

    if ( defined $data ) {
      $keep_trying = false;
    } else {
      print STDERR "Failure of '$url'\n";
    }
  }
  if ( $keep_trying ) {
    return @working;
  }

  if ( $data eq "" ) {
    print STDERR "No get_uids output returned from '$url'\n";
    return @working;
  }

  if ( $debug ) {
    print STDERR "$data\n";
  }

  my @ids = ($data =~ /<Id>(\d+)<\/Id>/g);
  foreach $uid (@ids) {
    push (@working, $uid);
  }

  return @working;
}

# send actual query

sub do_post_yielding_ref {

  my $urlx = shift (@_);
  my $argx = shift (@_);
  my $tulx = shift (@_);
  my $emlx = shift (@_);
  my $intr = shift (@_);

  if ( $os ne "" ) {
    $argx .= "&edirect_os=$os";
  }

  if ( $api_key ne "" ) {
    $argx .= "&api_key=$api_key";
  }

  $argx .= "&edirect=$version";

  if ( $intr ) {
    if ( $tulx eq "" ) {
      $tulx = "edirect";
    }
    if ( $tulx ne "" ) {
      $tulx =~ s/\s+/\+/g;
      $argx .= "&tool=$tulx";
    }

    if ( $emlx eq "" ) {
      $emlx = get_email ();
    }
    if ( $emlx ne "" ) {
      $emlx =~ s/\s+/\+/g;
      $argx .= "&email=$emlx";
    }
  }

  my $empty = '';
  $rslt = \$empty;

  if ( $debug or $log ) {
    print STDERR "$urlx?$argx\n";
  }

  if ( $http eq "get" or $http eq "GET" ) {
    if ( $argx ne "" ) {
      $urlx .= "?";
      $urlx .= "$argx";
    }

    $rslt = get ($urlx);

    if ( ! defined $rslt ) {
      print STDERR "Failure of '$urlx'\n";
      return "";
    }

    if ( $rslt eq "" ) {
      print STDERR "No do_get output returned from '$urlx'\n";
      return "";
    }

    if ( $debug ) {
      print STDERR "$rslt\n";
    }

    return \$rslt;
  }

  $usragnt = new LWP::UserAgent (timeout => 300);
  $usragnt->env_proxy;

  $req = new HTTP::Request POST => "$urlx";
  $req->content_type('application/x-www-form-urlencoded');
  $req->content("$argx");

  $res = $usragnt->request ( $req );

  if ( $res->is_success) {
    $rslt = $res->content_ref;
  } else {
    print STDERR $res->status_line . "\n";
  }

  if ( $$rslt eq "" ) {
    if ( $argx ne "" ) {
      $urlx .= "?";
      $urlx .= "$argx";
    }
    print STDERR "No do_post output returned from '$urlx'\n";
    print STDERR "Result of do_post http request is\n";
    print STDERR Dumper($res);
    print STDERR "\n";
  }

  if ( $debug ) {
    print STDERR $$rslt, "\n";
  }

  return $rslt;
}

sub do_post {
  my $rslt = do_post_yielding_ref(@_);
  return $$rslt;
}

# read ENTREZ_DIRECT data structure

sub read_edirect {

  my $dbsx = "";
  my $webx = "";
  my $keyx = "";
  my $numx = "";
  my $stpx = "";
  my $errx = "";
  my $tulx = "";
  my $emlx = "";
  my $inlabel = false;
  my $inmacro = false;
  my $ky = "";
  my $vl = "";

  @other = ();
  $has_num = false;
  $all_num = true;

  while ( defined($thisline = <STDIN>) ) {
    $thisline =~ s/\r//;
    $thisline =~ s/\n//;
    $thisline =~ s/^\s+//;
    $thisline =~ s/\s+$//;
    if ( $thisline =~ /<Labels>/ ) {
      $inlabel = true;
    } elsif ( $thisline =~ /<\/Labels>/ ) {
      $inlabel = false;
    } elsif ( $thisline =~ /<Macros>/ ) {
      $inmacro = true;
    } elsif ( $thisline =~ /<\/Macros>/ ) {
      $inmacro = false;
    } elsif ( $inlabel ) {
      if ( $thisline =~ /<Label>/ ) {
        $ky = "";
        $vl = "";
      } elsif ( $thisline =~ /<\/Label>/ ) {
        if ( $vl ne "" and $ky ne "" ) {
          $labels{"$ky"} = "$vl";
        }
      } elsif ( $thisline =~ /<Key>(.+)<\/Key>/ ) {
        $ky = $1;
      } elsif ( $thisline =~ /<Val>(.+)<\/Val>/ ) {
        $vl = $1;
      }
    } elsif ( $inmacro ) {
      if ( $thisline =~ /<Macro>/ ) {
        $ky = "";
        $vl = "";
      } elsif ( $thisline =~ /<\/Macro>/ ) {
        if ( $vl ne "" and $ky ne "" ) {
          $macros{"$ky"} = "$vl";
        }
      } elsif ( $thisline =~ /<Key>(.+)<\/Key>/ ) {
        $ky = $1;
      } elsif ( $thisline =~ /<Val>(.+)<\/Val>/ ) {
        $vl = $1;
      }
    } elsif ( $thisline =~ /<Db>(\S+)<\/Db>/ ) {
      $dbsx = $1;
    } elsif ( $thisline =~ /<WebEnv>(\S+)<\/WebEnv>/ ) {
      $webx = $1;
    } elsif ( $thisline =~ /<QueryKey>(\S+)<\/QueryKey>/ ) {
      $keyx = $1;
    } elsif ( $thisline =~ /<Count>(\S+)<\/Count>/ ) {
      $numx = $1;
    } elsif ( $thisline =~ /<Step>(\S+)<\/Step>/ ) {
      $stpx = $1 + 1;
    } elsif ( $thisline =~ /<Error>(.+?)<\/Error>/i ) {
      $errx = $1;
    } elsif ( $thisline =~ /<Tool>(.+?)<\/Tool>/i ) {
      $tulx = $1;
    } elsif ( $thisline =~ /<Email>(.+?)<\/Email>/i ) {
      $emlx = $1;
    } elsif ( $thisline =~ /<Silent>Y<\/Silent>/i ) {
      $silent = true;
    } elsif ( $thisline =~ /<Silent>N<\/Silent>/i ) {
      $silent = false;
    } elsif ( $thisline =~ /<Verbose>Y<\/Verbose>/i ) {
      $verbose = true;
    } elsif ( $thisline =~ /<Verbose>N<\/Verbose>/i ) {
      $verbose = false;
    } elsif ( $thisline =~ /<Debug>Y<\/Debug>/i ) {
      $debug = true;
      $silent = false;
    } elsif ( $thisline =~ /<Debug>N<\/Debug>/i ) {
      $debug = false;
    } elsif ( $thisline =~ /<Log>Y<\/Log>/i ) {
      $log = true;
    } elsif ( $thisline =~ /<Log>N<\/Log>/i ) {
      $log = false;
    } elsif ( $thisline =~ /<ENTREZ_DIRECT>/i ) {
    } elsif ( $thisline =~ /<\/ENTREZ_DIRECT>/i ) {
    } elsif ( $thisline =~ /<Labels>/i ) {
    } elsif ( $thisline =~ /<\/Labels>/i ) {
    } elsif ( $thisline =~ /<Macros>/i ) {
    } elsif ( $thisline =~ /<\/Macros>/i ) {
    } elsif ( $thisline =~ /^(\d+)$/ ) {
      push (@other, $1);
      $has_num = true;
    } elsif ( $thisline =~ /^(\d+).\d+$/ ) {
      push (@other, $1);
      $has_num = true;
    } elsif ( $thisline =~ /^(.+)$/ ) {
      push (@other, $1);
      $all_num = false;
    }
  }

  return ( $dbsx, $webx, $keyx, $numx, $stpx, $errx,
           $tulx, $emlx, $has_num && $all_num, @other );
}

# write ENTREZ_DIRECT data structure

sub write_edirect {

  my $dbsx = shift (@_);
  my $webx = shift (@_);
  my $keyx = shift (@_);
  my $numx = shift (@_);
  my $stpx = shift (@_);
  my $errx = shift (@_);
  my $tulx = shift (@_);
  my $emlx = shift (@_);

  my $seconds = "";
  my $end_time = Time::HiRes::time();
  my $elapsed = $end_time - $begin_time;
  if ( $elapsed > 0.0005 ) {
    if ( $elapsed =~ /(\d+)\.(\d\d\d)\d+/ ) {
      $seconds = "$1.$2";
    }
  }

  if ( $stpx eq "" and $webx ne "" ) {
    $stpx = "1";
  }

  if ( $tulx ne "" ) {
    $tulx =~ s/\s+/\+/g;
  }
  if ( $emlx ne "" ) {
    $emlx =~ s/\s+/\+/g;
  }

  if ( $verbose ) {
    print STDERR "\n";
    print STDERR "edirutil";

    if ( $dbsx ne "" ) {
      print STDERR " -db $dbsx";
    }
    if ( $webx ne "" ) {
      print STDERR " -web $webx";
    }
    if ( $keyx ne "" ) {
      print STDERR " -key $keyx";
    }
    if ( $numx ne "" ) {
      print STDERR " -count $numx";
    }
    if ( $stpx ne "" ) {
      print STDERR " -step $stpx";
    }
    if ( $seconds ne "" ) {
      print STDERR " -seconds $seconds";
    }

    print STDERR "\n\n";
  }

  if ( false ) {
    print STDERR "\n";
    print STDERR "<ENTREZ_DIRECT>\n";

    if ( $dbsx ne "" ) {
      print STDERR "  <Db>$dbsx</Db>\n";
    }
    if ( $webx ne "" ) {
      print STDERR "  <WebEnv>$webx</WebEnv>\n";
    }
    if ( $keyx ne "" ) {
      print STDERR "  <QueryKey>$keyx</QueryKey>\n";
    }
    if ( $numx ne "" ) {
      print STDERR "  <Count>$numx</Count>\n";
    }
    if ( $stpx ne "" ) {
      print STDERR "  <Step>$stpx</Step>\n";
    }
    if ( $errx ne "" ) {
      print STDERR "  <Error>$errx</Error>\n";
    }
    if ( $tulx ne "" ) {
      print STDERR "  <Tool>$tulx</Tool>\n";
    }
    if ( $emlx ne "" ) {
      print STDERR "  <Email>$emlx</Email>\n";
    }

    print STDERR "</ENTREZ_DIRECT>\n";
    print STDERR "\n";
  }

  print "<ENTREZ_DIRECT>\n";

  if ( $dbsx ne "" ) {
    print "  <Db>$dbsx</Db>\n";
  }
  if ( $webx ne "" ) {
    print "  <WebEnv>$webx</WebEnv>\n";
  }
  if ( $keyx ne "" ) {
    print "  <QueryKey>$keyx</QueryKey>\n";
  }
  if ( $numx ne "" ) {
    print "  <Count>$numx</Count>\n";
  }
  if ( $stpx ne "" ) {
    print "  <Step>$stpx</Step>\n";
  }
  if ( $errx ne "" ) {
    print "  <Error>$errx</Error>\n";
  }
  if ( $tulx ne "" ) {
    print "  <Tool>$tulx</Tool>\n";
  }
  if ( $emlx ne "" ) {
    print "  <Email>$emlx</Email>\n";
  }
  if ( $silent ) {
    print "  <Silent>Y</Silent>\n";
  }
  if ( $verbose ) {
    print "  <Verbose>Y</Verbose>\n";
  }
  if ( $debug ) {
    print "  <Debug>Y</Debug>\n";
  }
  if ( $log ) {
    print "  <Log>Y</Log>\n";
  }
  if ( scalar (keys %labels) > 0 ) {
    print "  <Labels>\n";
    for ( keys %labels ) {
      print "    <Label>\n";
      print "      <Key>$_</Key>\n";
      print "      <Val>$labels{$_}</Val>\n";
      print "    </Label>\n";
    }
    print "  </Labels>\n";
  }
  if ( scalar (keys %macros) > 0 ) {
    print "  <Macros>\n";
    for ( keys %macros ) {
      print "    <Macro>\n";
      print "      <Key>$_</Key>\n";
      print "      <Val>$macros{$_}</Val>\n";
      print "    </Macro>\n";
    }
    print "  </Macros>\n";
  }

  print "</ENTREZ_DIRECT>\n";
}

# wrapper to detect command line errors

sub MyGetOptions {

  my $help_msg = shift @_;

  if ( !GetOptions(@_) ) {
    die $help_msg;
  } elsif (@ARGV) {
    die ("Entrez Direct does not support positional arguments.\n"
         . "Please remember to quote parameter values containing\n"
         . "whitespace or shell metacharacters.\n");
  }
}

# subroutines for each -function

# ecntc prepares the requested tool and email arguments for an EUtils pipe

my $cntc_help = qq{
  -email    Contact person's address
  -tool     Name of script or program

};

sub ecntc {

  # ... | edirect.pl -contact -email darwin@beagle.edu -tool edirect_test | ...

  clearflags ();

  MyGetOptions(
    $cntc_help,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "econtact $version\n";
    print $cntc_help;
    return;
  }

  if ( -t STDIN and not @ARGV ) {
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
}

# efilt performs ESearch query refinement on the EUtils history server

my $filt_help = qq{
Query Specification

  -query       Query string

Document Order

  -sort        Result presentation order

Date Constraint

  -days        Number of days in the past
  -datetype    Date field abbreviation
  -mindate     Start of date range
  -maxdate     End of date range

Limit by Field

  -field       Query words individually in field
  -pairs       Query overlapping word pairs

Spell Check

  -spell       Correct misspellings in query

Publication Filters

  -pub         abstract, clinical, english, free, historical,
               journal, last_week, last_month, last_year,
               medline, preprint, published, review, structured

Sequence Filters

  -feature     gene, mrna, cds, mat_peptide, ...
  -location    mitochondrion, chloroplast, plasmid, plastid
  -molecule    genomic, mrna, trna, rrna, ncrna
  -organism    animals, archaea, bacteria, eukaryotes, fungi,
               human, insects, mammals, plants, prokaryotes,
               protists, rodents, viruses
  -source      genbank, insd, pdb, pir, refseq, swissprot, tpa

Gene Filters

  -status      alive
  -type        coding, pseudo

Miscellaneous Arguments

  -label       Alias for query step

};

sub process_extras {

  my $frst = shift (@_);
  my $publ = shift (@_);
  my $fkey = shift (@_);
  my $locn = shift (@_);
  my $bmol = shift (@_);
  my $orgn = shift (@_);
  my $sorc = shift (@_);
  my $stat = shift (@_);
  my $gtyp = shift (@_);

  $publ = lc($publ);
  $fkey = lc($fkey);
  $bmol = lc($bmol);
  $locn = lc($locn);
  $orgn = lc($orgn);
  $sorc = lc($sorc);
  $stat = lc($stat);
  $gtyp = lc($gtyp);

  %pubHash = (
    'abstract'     =>  'has abstract [FILT]',
    'clinical'     =>  'clinical trial [FILT]',
    'english'      =>  'english [FILT]',
    'free'         =>  'freetext [FILT]',
    'historical'   =>  'historical article [FILT]',
    'journal'      =>  'journal article [FILT]',
    'last_month'   =>  'published last month [FILT]',
    'last month'   =>  'published last month [FILT]',
    'last_week'    =>  'published last week [FILT]',
    'last week'    =>  'published last week [FILT]',
    'last_year'    =>  'published last year [FILT]',
    'last year'    =>  'published last year [FILT]',
    'medline'      =>  'medline [FILT]',
    'preprint'     =>  'ahead of print [FILT]',
    'review'       =>  'review [FILT]',
    'structured'   =>  'hasstructuredabstract [WORD]',
    'trial'        =>  'clinical trial [FILT]',
  );

  @featureArray = (
    "-10_signal",
    "-35_signal",
    "3'clip",
    "3'utr",
    "5'clip",
    "5'utr",
    "allele",
    "assembly_gap",
    "attenuator",
    "c_region",
    "caat_signal",
    "cds",
    "centromere",
    "conflict",
    "d_segment",
    "d-loop",
    "enhancer",
    "exon",
    "gap",
    "gc_signal",
    "gene",
    "idna",
    "intron",
    "j_segment",
    "ltr",
    "mat_peptide",
    "misc_binding",
    "misc_difference",
    "misc_feature",
    "misc_recomb",
    "misc_rna",
    "misc_signal",
    "misc_structure",
    "mobile_element",
    "modified_base",
    "mrna",
    "mutation",
    "n_region",
    "ncrna",
    "old_sequence",
    "operon",
    "orit",
    "polya_signal",
    "polya_site",
    "precursor_rna",
    "prim_transcript",
    "primer_bind",
    "promoter",
    "propeptide",
    "protein_bind",
    "rbs",
    "regulatory",
    "rep_origin",
    "repeat_region",
    "repeat_unit",
    "rrna",
    "s_region",
    "satellite",
    "scrna",
    "sig_peptide",
    "snorna",
    "snrna",
    "source",
    "stem_loop",
    "sts",
    "tata_signal",
    "telomere",
    "terminator",
    "tmrna",
    "transit_peptide",
    "trna",
    "unsure",
    "v_region",
    "v_segment",
    "variation"
  );

  %locationHash = (
    'mitochondria'   =>  'mitochondrion [FILT]',
    'mitochondrial'  =>  'mitochondrion [FILT]',
    'mitochondrion'  =>  'mitochondrion [FILT]',
    'chloroplast'    =>  'chloroplast [FILT]',
    'plasmid'        =>  'plasmid [FILT]',
    'plastid'        =>  'plastid [FILT]',
  );

  %moleculeHash = (
    'genomic'  =>  'biomol genomic [PROP]',
    'mrna'     =>  'biomol mrna [PROP]',
    'trna'     =>  'biomol trna [PROP]',
    'rrna'     =>  'biomol rrna [PROP]',
    'ncrna'    =>  'biomol ncrna [PROP]',
  );

  %organismHash = (
    'animal'           =>  'animals [FILT]',
    'animals'          =>  'animals [FILT]',
    'archaea'          =>  'archaea [FILT]',
    'archaeal'         =>  'archaea [FILT]',
    'archaean'         =>  'archaea [FILT]',
    'archaebacteria'   =>  'archaea [FILT]',
    'archaebacterial'  =>  'archaea [FILT]',
    'bacteria'         =>  'bacteria [FILT]',
    'bacterial'        =>  'bacteria [FILT]',
    'bacterium'        =>  'bacteria [FILT]',
    'eubacteria'       =>  'bacteria [FILT]',
    'eubacterial'      =>  'bacteria [FILT]',
    'eukaryota'        =>  'eukaryota [ORGN]',
    'eukaryote'        =>  'eukaryota [ORGN]',
    'eukaryotes'       =>  'eukaryota [ORGN]',
    'fungal'           =>  'fungi [FILT]',
    'fungi'            =>  'fungi [FILT]',
    'fungus'           =>  'fungi [FILT]',
    'human'            =>  'human [ORGN]',
    'humans'           =>  'human [ORGN]',
    'insect'           =>  'insecta [ORGN]',
    'insecta'          =>  'insecta [ORGN]',
    'insects'          =>  'insecta [ORGN]',
    'mammal'           =>  'mammals [FILT]',
    'mammalia'         =>  'mammals [FILT]',
    'mammalian'        =>  'mammals [FILT]',
    'mammals'          =>  'mammals [FILT]',
    'man'              =>  'human [ORGN]',
    'metaphyta'        =>  'plants [FILT]',
    'metazoa'          =>  'animals [FILT]',
    'monera'           =>  'prokaryota [ORGN]',
    'plant'            =>  'plants [FILT]',
    'plants'           =>  'plants [FILT]',
    'prokaryota'       =>  'prokaryota [ORGN]',
    'prokaryote'       =>  'prokaryota [ORGN]',
    'prokaryotes'      =>  'prokaryota [ORGN]',
    'protist'          =>  'protists [FILT]',
    'protista'         =>  'protists [FILT]',
    'protists'         =>  'protists [FILT]',
    'rodent'           =>  'rodents [ORGN]',
    'rodentia'         =>  'rodents [ORGN]',
    'rodents'          =>  'rodents [ORGN]',
    'viral'            =>  'viruses [FILT]',
    'virus'            =>  'viruses [FILT]',
    'viruses'          =>  'viruses [FILT]',
  );

  %sourceHash = (
    'ddbj'       =>  'srcdb ddbj [PROP]',
    'embl'       =>  'srcdb embl [PROP]',
    'genbank'    =>  'srcdb genbank [PROP]',
    'insd'       =>  'srcdb ddbj/embl/genbank [PROP]',
    'pdb'        =>  'srcdb pdb [PROP]',
    'pir'        =>  'srcdb pir [PROP]',
    'refseq'     =>  'srcdb refseq [PROP]',
    'swissprot'  =>  'srcdb swiss prot [PROP]',
    'tpa'        =>  'srcdb tpa ddbj/embl/genbank [PROP]',
  );

  %statusHash = (
    'alive'   =>  'alive [PROP]',
    'live'    =>  'alive [PROP]',
    'living'  =>  'alive [PROP]',
  );

  %typeHash = (
    'coding'  =>  'genetype protein coding [PROP]',
    'pseudo'  =>  'genetype pseudo [PROP]',
  );

  my @working = ();

  my $suffix = "";

  if ( $frst ne "" ) {
    push (@working, $frst);
  }

  if ( $publ ne "" ) {
    if ( defined $pubHash{$publ} ) {
      $val = $pubHash{$publ};
      push (@working, $val);
    } elsif ( $publ eq "published" ) {
      $suffix = "published";
    } else {
      die "\nUnrecognized -pub argument '$publ', use efilter -help to see available choices\n\n";
    }
  }

  if ( $fkey ne "" ) {
    if ( grep( /^$fkey$/, @featureArray ) ) {
      $val = $fkey . " [FKEY]";
      push (@working, $val);
    } else {
      die "\nUnrecognized -feature argument '$fkey', use efilter -help to see available choices\n\n";
    }
  }

  if ( $locn ne "" ) {
    if ( defined $locationHash{$locn} ) {
      $val = $locationHash{$locn};
      push (@working, $val);
    } else {
      die "\nUnrecognized -location argument '$locn', use efilter -help to see available choices\n\n";
    }
  }

  if ( $bmol ne "" ) {
    if ( defined $moleculeHash{$bmol} ) {
      $val = $moleculeHash{$bmol};
      push (@working, $val);
    } else {
      die "\nUnrecognized -molecule argument '$bmol', use efilter -help to see available choices\n\n";
    }
  }

  if ( $orgn ne "" ) {
    if ( defined $organismHash{$orgn} ) {
      $val = $organismHash{$orgn};
      push (@working, $val);
    } else {
      die "\nUnrecognized -organism argument '$orgn', use efilter -help to see available choices\n\n";
    }
  }

  if ( $sorc ne "" ) {
    if ( defined $sourceHash{$sorc} ) {
      $val = $sourceHash{$sorc};
      push (@working, $val);
    } else {
      die "\nUnrecognized -source argument '$sorc', use efilter -help to see available choices\n\n";
    }
  }

  if ( $stat ne "" ) {
    if ( defined $statusHash{$stat} ) {
      $val = $statusHash{$stat};
      push (@working, $val);
    } else {
      die "\nUnrecognized -status argument '$stat', use efilter -help to see available choices\n\n";
    }
  }

  if ( $gtyp ne "" ) {
    if ( defined $typeHash{$gtyp} ) {
      $val = $typeHash{$gtyp};
      push (@working, $val);
    } else {
      die "\nUnrecognized -type argument '$gtyp', use efilter -help to see available choices\n\n";
    }
  }

  my $xtras = join (" AND ", @working);

  if ( $suffix eq "published" ) {
    $xtras = $xtras . " NOT ahead of print [FILT]";
  }

  return $xtras;
}

# correct misspellings in query

sub spell_check_query {

  my $db = shift (@_);
  my $qury = shift (@_);

  my $url = $base . $espell;

  my $enc = uri_escape($query);
  $arg = "db=$db&term=$enc";

  my $data = do_post ($url, $arg, $tool, $email, true);

  Encode::_utf8_on($data);

  $qury = $1 if ( $data =~ /<CorrectedQuery>(.+)<\/CorrectedQuery>/ );

  return $qury;
}

sub efilt {

  # ... | edirect.pl -filter -query "bacteria [ORGN]" -days 365 | ...

  clearflags ();

  MyGetOptions(
    $filt_help,
    "query=s" => \$query,
    "sort=s" => \$sort,
    "days=i" => \$rldate,
    "mindate=s" => \$mndate,
    "maxdate=s" => \$mxdate,
    "datetype=s" => \$dttype,
    "label=s" => \$lbl,
    "field=s" => \$field,
    "spell" => \$spell,
    "pairs=s" => \$pair,
    "pub=s" => \$pub,
    "feature=s" => \$feature,
    "location=s" => \$location,
    "molecule=s" => \$molecule,
    "organism=s" => \$organism,
    "source=s" => \$source,
    "status=s" => \$status,
    "type=s" => \$gtype,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "efilter $version\n";
    print $filt_help;
    return;
  }

  # process special filter flags, add to query string
  $query = process_extras ( $query, $pub, $feature, $location, $molecule, $organism, $source, $status, $gtype );

  if ( -t STDIN ) {
    if ( $query eq "" ) {
      die "Must supply -query or -days or -mindate and -maxdate arguments on command line\n";
    }
    print "efilter -query \"$query\"\n";
    return;
  }

  ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    write_edirect ( "", "", "", "", "", $err, "", "" );
    close (STDOUT);
    if ( ! $silent ) {
      die "ERROR in filt input: $err\n\n";
    }
    return;
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  if ( $query eq "" && $rldate < 1 and $mndate eq "" and $mxdate eq "" ) {
    die "Must supply -query or -days or -mindate and -maxdate arguments on command line\n";
  }

  binmode STDOUT, ":utf8";

  if ( $dbase ne "" and $web ne "" and $key eq "" and $num eq "0" ) {
    write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
    close (STDOUT);
    die "QueryKey value not found in filter input\n";
    return;
  }

  # warn on mismatch between filter argument and database
  if ( $dbase ne "pubmed" ) {
    if ( $pub ne "" ) {
      print STDERR "\nUnexpected use of pubmed filter argument\n\n";
    }
  }
  if ( $dbase ne "nucleotide" and
       $dbase ne "nuccore" and
       $dbase ne "est" and
       $dbase ne "gss" and
       $dbase ne "protein" ) {
    if ( $feature ne "" or
         $location ne "" or
         $molecule ne "" or
         $organism ne "" or
         $source ne "" ) {
      print STDERR "\nUnexpected use of sequence filter argument\n\n";
    }
  }

  test_edirect ( $dbase, $web, $key, $num, "filter" );

  # -field combines -drop and -split (-field TITL produces same behavior as Web PubMed)
  if ( $field ne "" ) {
    $query = remove_stop_words ($query);
    $query = field_each_word ($field, $query);
  }

  # -pairs separately fields query word pairs, breaking chain at stop words
  if ( $pair ne "" ) {
    $query = remove_punctuation ($query);
    if ( $query =~ /^ +(.+)$/ ) {
      $query = $1;
    }
    if ( $query =~ /^(.+) +$/ ) {
      $query = $1;
    }
    $query =~ s/ +/ /g;
    $query = field_each_pair ($pair, $query);
  }

  # spell check each query word
  if ( $spell ) {
    $query = spell_check_query ($dbase, $query);
  }

  $url = $base . $esearch;

  $arg = "db=$dbase&query_key=$key&WebEnv=$web";
  if ( $sort ne "" ) {
    if ( $sort eq "Relevance" ) {
      $sort = "relevance";
    }
    $arg .= "&sort=$sort";
  }
  $arg .= "&retmax=0&usehistory=y";
  if ( $query ne "" ) {
    $query = map_labels ($query);
    $query = map_macros ($query);
    $enc = uri_escape($query);
    $arg .= "&term=$enc";
  }
  if ( $rldate > 0 ) {
    $arg .= "&reldate=$rldate";
    if ( $dttype eq "" ) {
      $dttype = "PDAT";
    }
  }
  if ( $mndate ne "" and $mxdate ne "" ) {
    $arg .= "&mindate=$mndate&maxdate=$mxdate";
    if ( $dttype eq "" ) {
      $dttype = "PDAT";
    }
  }
  if ( $dttype ne "" ) {
    $arg .= "&datetype=$dttype";
  }

  $wb = $web;

  $web = "";
  $key = "";
  $num = "";
  $err = "";
  my $trn = "";

  $output = "";

  for ( $tries = 0; $tries < 3 && $web eq ""; $tries++) {
    $output = do_post ($url, $arg, $tool, $email, true);

    $err = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);
    if ( $err eq "" ) {
      $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
      $key = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
      $num = $1 if ($output =~ /<Count>(\S+)<\/Count>/);
      $trn = $1 if ($output =~ /<QueryTranslation>(.+?)<\/QueryTranslation>/i);
    } else {
      if ( ! $silent ) {
        print STDERR "Retrying efilter, step $stp: $err\n";
      }
      sleep 3;
    }
  }

  if ( $err ne "" ) {
    write_edirect ( "", "", "", "", "", $err, "", "" );
    close (STDOUT);
    if ( ! $silent) {
      die "ERROR in filt output: $err\nURL: $arg\nResult: $output\n\n";
    }
    return;
  }

  if ( $web eq "" ) {
    die "WebEnv value not found in filt output - WebEnv1 $wb\n";
  }
  if ( $key eq "" ) {
    die "QueryKey value not found in filt output - WebEnv1 $wb\n";
  }

  if ( $web ne $wb ) {
    $err = "WebEnv mismatch in filt output - WebEnv1 $wb, WebEnv2 $web";
    write_edirect ( "", $wb, "", "", "", $err, "", "" );
    close (STDOUT);
    die "WebEnv value changed in filt output - WebEnv1 $wb\nWebEnv2 $web\n";
  }

  if ( $num eq "" ) {
    die "Count value not found in filt output - WebEnv1 $wb\n";
  }

  if ( $lbl ne "" and $key ne "" ) {
    $labels{"$lbl"} = "$key";
  }

  write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );

  if ( $log ) {
    if ( $trn ne "" ) {
      print STDERR "$trn\n";
    }
  }
}

# efetch -format docsum calls esmry to retrieve document summaries

sub fix_mixed_content {

  my $x = shift (@_);

  while ( $x =~ /\&amp\;/ || $x =~ /\&lt\;/ || $x =~ /\&gt\;/ ) {
    HTML::Entities::decode_entities($x);
  }
  # removed mixed content tags
  $x =~ s|<b>||g;
  $x =~ s|<i>||g;
  $x =~ s|<u>||g;
  $x =~ s|<sup>||g;
  $x =~ s|<sub>||g;
  $x =~ s|</b>||g;
  $x =~ s|</i>||g;
  $x =~ s|</u>||g;
  $x =~ s|</sup>||g;
  $x =~ s|</sub>||g;
  $x =~ s|<b/>||g;
  $x =~ s|<i/>||g;
  $x =~ s|<u/>||g;
  $x =~ s|<sup/>||g;
  $x =~ s|<sub/>||g;
  # Reencode any resulting less-than or greater-than entities to avoid breaking the XML.
  $x =~ s/</&lt;/g;
  $x =~ s/>/&gt;/g;

  return $x;
}

my %fields_to_fix = (
  'biosample' => ['SampleData'],
  'medgen'    => ['ConceptMeta'],
  'sra'       => ['ExpXml', 'Runs']
);

sub fix_one_encoding {

  my $dbase = shift (@_);
  my $data = shift (@_);

  if ( $dbase eq "pubmed" ) {
    if ( $data =~ /<Title>(.+?)<\/Title>/ ) {
      my $x = $1;
      if ( $x =~ /\&amp\;/ || $x =~ /\&lt\;/ || $x =~ /\&gt\;/ || $x =~ /\</ || $x =~ /\>/ ) {
        $x = fix_mixed_content($x);
        $data =~ s/<Title>(.+?)<\/Title>/<Title>$x<\/Title>/;
      }
    }
  } elsif ( $dbase eq "gene" ) {
    if ( $data =~ /<Summary>(.+?)<\/Summary>/ ) {
      my $x = $1;
      if ( $x =~ /\&amp\;/ ) {
        HTML::Entities::decode_entities($x);
        # Reencode any resulting less-than or greater-than entities to avoid breaking the XML.
        $x =~ s/</&lt;/g;
        $x =~ s/>/&gt;/g;
        $data =~ s/<Summary>(.+?)<\/Summary>/<Summary>$x<\/Summary>/;
      }
    }
    # $data =~ s/(\s+?)<ChrStart>(\d+)<\/ChrStart>/$1<ChrStart>$2<\/ChrStart>$1<SeqStart>${\($2 + 1)}<\/SeqStart>/g;
    # $data =~ s/(\s+?)<ChrStop>(\d+)<\/ChrStop>/$1<ChrStop>$2<\/ChrStop>$1<SeqStop>${\($2 + 1)}<\/SeqStop>/g;
  } elsif ( $dbase eq "assembly" ) {
    if ( $data =~ /<Meta>(.+?)<\/Meta>/ ) {
      my $x = $1;
      if ( $x =~ /<!\[CDATA\[\s*(.+?)\s*\]\]>/ ) {
        $x = $1;
        if ( $x =~ /</ and $x =~ />/ ) {
            # If CDATA contains embedded XML, simply remove CDATA wrapper
          $data =~ s/<Meta>(.+?)<\/Meta>/<Meta>$x<\/Meta>/;
        }
      }
    }
  } elsif (defined $fields_to_fix{$dbase}) {
    foreach $f (@{$fields_to_fix{$dbase}}) {
      if ( $data =~ /<$f>(.+?)<\/$f>/ ) {
        my $x = $1;
        if ( $x =~ /\&lt\;/ and $x =~ /\&gt\;/ ) {
          HTML::Entities::decode_entities($x);
          $data =~ s/<$f>(.+?)<\/$f>/<$f>$x<\/$f>/;
        }
      }
    }
  }

  return $data;
}

sub fix_bad_encoding {

  my $dbase = shift (@_);
  my $data = shift (@_);

  if ( $dbase eq "pubmed" || $dbase eq "gene" || $dbase eq "assembly" || defined $fields_to_fix{$dbase} ) {

    my @accum = ();
    my @working = ();
    my $prefix = "";
    my $suffix = "";
    my $docsumset_attrs = '';

    if ( $data =~ /(.+?)<DocumentSummarySet(\s+.+?)?>(.+)<\/DocumentSummarySet>(.+)/s ) {
      $prefix = $1;
      $docsumset_attrs = $2;
      my $docset = $3;
      $suffix = $4;
      my @vals = ($docset =~ /<DocumentSummary>(.+?)<\/DocumentSummary>/sg);
      foreach $val (@vals) {
        push (@working, "<DocumentSummary>");
        push (@working, fix_one_encoding ( $dbase, $val) );
        push (@working, "</DocumentSummary>");
      }
    }

    if ( scalar @working > 0 ) {
      push (@accum, $prefix);
      push (@accum, "<DocumentSummarySet$docsumset_attrs>");
      push (@accum, @working);
      push (@accum, "</DocumentSummarySet>");
      push (@accum, $suffix);
      $data = join ("\n", @accum);
      $data =~ s/\n\n/\n/g;
    }
  }

  return $data;
}

sub esmry {

  my $dbase = shift (@_);
  my $web = shift (@_);
  my $key = shift (@_);
  my $num = shift (@_);
  my $id = shift (@_);
  my $mode = shift (@_);
  my $min = shift (@_);
  my $max = shift (@_);
  my $tool = shift (@_);
  my $email = shift (@_);
  my $silent = shift (@_);
  my $verbose = shift (@_);
  my $debug = shift (@_);
  my $log = shift (@_);
  my $http = shift (@_);
  my $alias = shift (@_);
  my $basx = shift (@_);

  $dbase = lc($dbase);

  if ( $dbase ne "" and $id ne "" ) {

    if ( $id eq "0" ) {

      # id "0" returns a live UID for any database

      $id = get_zero_uid ($dbase);

      if ( $id eq "0" ) {

        # id "0" is an unrecognized accession

        return;
      }
    }

    $url = $base . $esummary;

    if ( $dbase eq "pubmed" ) {
      $id =~ s/(\d+)\.\d+/$1/g;
    }
    $arg = "db=$dbase&id=$id";
    $arg .= "&version=2.0";
    if ( $mode ne "" and $mode ne "text" ) {
      $arg .= "&retmode=$mode";
    }

    $data = do_post ($url, $arg, $tool, $email, true);

    if ($data =~ /<eSummaryResult>/i and $data =~ /<ERROR>(.+?)<\/ERROR>/i) {
      $err = $1;
      if ( $err ne "" ) {
        if ( ! $silent ) {
          print STDERR "ERROR in esummary: $err\n";
        }
      }
      if ( $err =~ "Invalid uid" ) {
        # remove Invalid uid error message from XML
        $data =~ s/<ERROR>.+?<\/ERROR>//g;
      } else {
        return;
      }
    }

    if (! $raw) {
      if ($data !~ /<Id>\d+<\/Id>/i) {
        $data =~ s/<DocumentSummary uid=\"(\d+)\">/<DocumentSummary><Id>$1<\/Id>/g;
      }
    }

    Encode::_utf8_on($data);

    if (! $raw) {
      $data = fix_bad_encoding($dbase, $data);
    }

    # remove eSummaryResult wrapper
    $data =~ s/<!DOCTYPE eSummaryResult PUBLIC/<!DOCTYPE DocumentSummarySet PUBLIC/g;
    $data =~ s/<eSummaryResult>//g;
    $data =~ s/<\/eSummaryResult>//g;

    print "$data";

    return;
  }

  if ( $dbase ne "" and $web ne "" and $key eq "" and $num eq "0" ) {
    die "QueryKey value not found in summary input\n";
    return;
  }

  test_edirect ( $dbase, $web, $key, $num, "summary" );

  $stpminusone = $stp - 1;

  if ( $max == 0 ) {
    if ( $silent ) {
      return;
    }
  }

  # use larger chunk for document summaries
  $chunk = 1000;
  for ( $start = $min; $start < $max; $start += $chunk ) {
    $url = $base . $esummary;

    $chkx = $chunk;
    if ( $start + $chkx > $max ) {
      $chkx = $max - $start;
    }

    $arg = "db=$dbase&query_key=$key&WebEnv=$web";
    $arg .= "&retstart=$start&retmax=$chkx&version=2.0";
    if ( $mode ne "" and $mode ne "text" ) {
      $arg .= "&retmode=$mode";
    }

    $data = "";
    $retry = true;

    for ( $tries = 0; $tries < 3 && $retry; $tries++) {
      $data = do_post ($url, $arg, $tool, $email, true);

      if ($data =~ /<eSummaryResult>/i and $data =~ /<ERROR>(.+?)<\/ERROR>/i ) {
        if ( ! $silent ) {
          print STDERR "Retrying esummary, step $stp: $err\n";
        }
        sleep 3;
      } else {
        $retry = false;
      }
    }

    if ($data =~ /<eSummaryResult>/i and $data =~ /<ERROR>(.+?)<\/ERROR>/i) {
      $err = $1;
      if ( $err ne "" ) {
        if ( ! $silent ) {
          my $from = $start + 1;
          my $to = $start + $chunk;
          if ( $to > $max ) {
            $to = $max;
          }
          print STDERR "ERROR in esummary ($from-$to / $max): $err\n";
          print STDERR "Replicate for debugging with:\n";
          print STDERR "  edirutil -db $dbase -web $web -key $key -count $num";
          if ( $stpminusone > 0 ) {
            print STDERR " -step $stpminusone";
          }
          my $seconds = "";
          my $end_time = Time::HiRes::time();
          my $elapsed = $end_time - $begin_time;
          if ( $elapsed > 0.0005 ) {
            if ( $elapsed =~ /(\d+)\.(\d\d\d)\d+/ ) {
              $seconds = "$1.$2";
            }
          }
          if ( $seconds ne "" ) {
            print STDERR " -seconds $seconds";
          }
          print STDERR " | efetch -format docsum -start $from -stop $to\n";
        }
      }
    } else {
      if ( $verbose ) {
        my $from = $start + 1;
        my $to = $start + $chunk;
        if ( $to > $max ) {
          $to = $max;
        }
        print STDERR "( edirutil -db $dbase -web $web -key $key -count $num";
        if ( $stpminusone > 0 ) {
          print STDERR " -step $stpminusone";
        }
        my $seconds = "";
        my $end_time = Time::HiRes::time();
        my $elapsed = $end_time - $begin_time;
        if ( $elapsed > 0.0005 ) {
          if ( $elapsed =~ /(\d+)\.(\d\d\d)\d+/ ) {
            $seconds = "$1.$2";
          }
        }
        if ( $seconds ne "" ) {
          print STDERR " -seconds $seconds";
        }
        print STDERR " | efetch -format docsum -start $from -stop $to )\n";
        $begin_time = $end_time;
      }

      if (! $raw) {
        if ($data !~ /<Id>\d+<\/Id>/i) {
          $data =~ s/<DocumentSummary uid=\"(\d+)\">/<DocumentSummary><Id>$1<\/Id>/g;
        }
      }

      Encode::_utf8_on($data);

      if (! $raw) {
        $data = fix_bad_encoding($dbase, $data);
      }

      # remove eSummaryResult wrapper
      $data =~ s/<!DOCTYPE eSummaryResult PUBLIC/<!DOCTYPE DocumentSummarySet PUBLIC/g;
      $data =~ s/<eSummaryResult>//g;
      $data =~ s/<\/eSummaryResult>//g;

      print "$data";
    }

    do_sleep ();
  }
}

# eftch can read all arguments from the command line or participate in an EUtils pipe

my $ftch_help = qq{
Format Selection

  -format        Format of record or report
  -mode          text, xml, asn.1, json
  -style         withparts, conwithfeat

Direct Record Selection

  -db            Database name
  -id            Unique identifier or accession number

Sequence Range

  -seq_start     First sequence position to retrieve
  -seq_stop      Last sequence position to retrieve
  -strand        Strand of DNA to retrieve

Gene Range

  -chr_start     Sequence range from 0-based coordinates
  -chr_stop        in gene docsum GenomicInfoType object

Sequence Flags

  -complexity    0 = default, 1 = bioseq, 3 = nuc-prot set
  -extend        Extend sequence retrieval in both directions
  -extrafeat     Bit flag specifying extra features

Miscellaneous

  -raw           Skip database-specific XML modifications

Format Examples

  -db            -format            -mode    Report Type
  ___            _______            _____    ___________

  (all)
                 docsum                      DocumentSummarySet XML
                 docsum             json     DocumentSummarySet JSON
                 full                        Same as native except for mesh
                 uid                         Unique Identifier List
                 url                         Entrez URL
                 xml                         Same as -format full -mode xml

  bioproject
                 native                      BioProject Report
                 native             xml      RecordSet XML

  biosample
                 native                      BioSample Report
                 native             xml      BioSampleSet XML

  biosystems
                 native             xml      Sys-set XML

  gds
                 native             xml      RecordSet XML
                 summary                     Summary

  gene
                 gene_table                  Gene Table
                 native                      Gene Report
                 native             asn.1    Entrezgene ASN.1
                 native             xml      Entrezgene-Set XML
                 tabular                     Tabular Report

  homologene
                 alignmentscores             Alignment Scores
                 fasta                       FASTA
                 homologene                  Homologene Report
                 native                      Homologene List
                 native             asn.1    HG-Entry ASN.1
                 native             xml      Entrez-Homologene-Set XML

  mesh
                 full                        Full Record
                 native                      MeSH Report
                 native             xml      RecordSet XML

  nlmcatalog
                 native                      Full Record
                 native             xml      NLMCatalogRecordSet XML

  pmc
                 medline                     MEDLINE
                 native             xml      pmc-articleset XML

  pubmed
                 abstract                    Abstract
                 medline                     MEDLINE
                 native             asn.1    Pubmed-entry ASN.1
                 native             xml      PubmedArticleSet XML

  (sequences)
                 acc                         Accession Number
                 est                         EST Report
                 fasta                       FASTA
                 fasta              xml      TinySeq XML
                 fasta_cds_aa                FASTA of CDS Products
                 fasta_cds_na                FASTA of Coding Regions
                 ft                          Feature Table
                 gb                          GenBank Flatfile
                 gb                 xml      GBSet XML
                 gbc                xml      INSDSet XML
                 gene_fasta                  FASTA of Gene
                 gp                          GenPept Flatfile
                 gp                 xml      GBSet XML
                 gpc                xml      INSDSet XML
                 gss                         GSS Report
                 ipg                         Identical Protein Report
                 ipg                xml      IPGReportSet XML
                 native             text     Seq-entry ASN.1
                 native             xml      Bioseq-set XML
                 seqid                       Seq-id ASN.1

  snp
                 chr                         Chromosome Report
                 docset                      Summary
                 fasta                       FASTA
                 flt                         Flat File
                 native             asn.1    Rs ASN.1
                 native             xml      ExchangeSet XML
                 rsr                         RS Cluster Report
                 ssexemplar                  SS Exemplar List

  sra
                 native             xml      EXPERIMENT_PACKAGE_SET XML
                 runinfo            xml      SraRunInfo XML

  structure
                 mmdb                        Ncbi-mime-asn1 strucseq ASN.1
                 native                      MMDB Report
                 native             xml      RecordSet XML

  taxonomy
                 native                      Taxonomy List
                 native             xml      TaxaSet XML

};

sub fix_sra_xml_encoding {

  my $data = shift (@_);

  $data =~ s/<!--[^<]+</</g;
  $data =~ s/>\s*-->/>/g;

  return $data;
}

sub fix_pubmed_xml_encoding {

  my $x = shift (@_);

  my $markup = '(?:[biu]|su[bp])';
  my $attrs = ' ?';
  # my $markup = '(?:(?:[\w.:_-]*:)?[[:lower:]-]+|DispFormula)';
  # my $attrs = '(?:\s[^>]*)?';

  # check for possible newline artifact
  $x =~ s|</$markup>\n||g;
  $x =~ s|\n<$markup$attrs>||g;

  # removed mixed content tags
  $x =~ s|</$markup>||g;
  $x =~ s|<$markup$attrs/?>||g;

  # check for encoded tags
  if ( $x =~ /\&amp\;/ || $x =~ /\&lt\;/ || $x =~ /\&gt\;/ ) {
    # remove runs of amp
    $x =~ s|&amp;(?:amp;)+|&amp;|g;
    # fix secondary encoding
    $x =~ s|&amp;lt;|&lt;|g;
    $x =~ s|&amp;gt;|&gt;|g;
    # temporarily protect encoded scientific symbols, e.g., PMID 9698410 and 21892341
    $x =~ s|(?<= )(&lt;)(=*$markup&gt;)(?= )|$1=$2|g;
    # remove encoded markup
    $x =~ s|&lt;/$markup&gt;||g;
    $x =~ s|&lt;$markup$attrs/?&gt;||g;
    # undo temporary protection of scientific symbols adjacent to space
    $x =~ s|(?<= )(&lt;)=(=*$markup&gt;)(?= )|$1$2|g;
  }

  # compress runs of horizontal whitespace
  $x =~ s/\h+/ /g;

  # remove spaces just outside of angle brackets
  $x =~ s|> |>|g;
  $x =~ s| <|<|g;

  # remove spaces just inside of parentheses
  $x =~ s|\( |\(|g;
  $x =~ s| \)|\)|g;

  return $x;
}

# for id in 9698410 16271163 17282049 20968289 21892341 22785267 25435818 27672066 28635620 28976125 29547395
# do
#   efetch -db pubmed -format xml -id "$id" |
#   xtract -pattern PubmedArticle -plg "\n\n" -sep "\n\n" -tab "\n\n" \
#     -element MedlineCitation/PMID ArticleTitle AbstractText
# done

sub eftch {

  # ... | edirect.pl -fetch -format gp | ...

  clearflags ();

  MyGetOptions(
    $ftch_help,
    "db=s" => \$db,
    "id=s" => \$id,
    "format=s" => \$type,
    "style=s" => \$style,
    "mode=s" => \$mode,
    "seq_start=i" => \$seq_start,
    "seq_stop=i" => \$seq_stop,
    "strand=s" => \$strand,
    "complexity=i" => \$complexity,
    "chr_start=i" => \$chr_start,
    "chr_stop=i" => \$chr_stop,
    "showgi" => \$showgi,
    "extend=i" => \$extend,
    "extrafeat=i" => \$extrafeat,
    "start=i" => \$min,
    "stop=i" => \$max,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "pipe" => \$pipe,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "raw" => \$raw,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "efetch $version\n";
    print $ftch_help;
    return;
  }

  # convert spaces between UIDs to commas

  $id =~ s/ /,/g;
  $id =~ s/,,/,/g;

  # "-format xml" is a shortcut for "-format full -mode xml"

  if ( $type eq "xml" and $mode eq "" ) {
    $type = "full";
    $mode = "xml";
  }

  if ( $style eq "normal" or $style eq "none" ) {
    $style = "";
  }

  if ( $style eq "conwithfeats" or $style eq "gbconwithfeat" or $style eq "gbconwithfeats" ) {
    $style = "conwithfeat";
  } elsif ( $style eq "withpart" or $style eq "gbwithpart" or $style eq "gbwithparts" ) {
    $style = "withparts";
  }

  if ( $type eq "gbconwithfeat" or $type eq "gbconwithfeats" ) {
    $type = "gb";
    $style = "conwithfeat";
  } elsif ( $type eq "gbwithparts" or $type eq "gbwithpart" ) {
    $type = "gb";
    $style = "withparts";
  }

  if ( -t STDIN and not @ARGV ) {
  } elsif ( $db ne "" and $id ne "" ) {
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    if ( ! $silent ) {
      die "ERROR in fetch input: $err\n\n";
    }
    return;
  }

  if ( $dbase eq "" ) {
    $dbase = $db;
  }

  $dbase = lc($dbase);

  if ( $type eq "" and $db ne "" ) {
    if ( get_zero_uid ($db) eq "" ) {
      die "Must supply -format report type on command line\n";
    }
    $type = "native";
  }

  if ( $mode eq "" ) {
    $mode = "text";
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  if ( $strand eq "plus" ) {
    $strand = "1";
  }
  if ( $strand eq "minus" ) {
    $strand = "2";
  }

  binmode STDOUT, ":utf8";

  # arguments can override loop start and stop

  if ( $min > 0 ) {
    $min--;
  }
  if ( $max == 0 ) {
    $max = $num
  }

  if ( $type eq "docsum" or $fnc eq "-summary" ) {

    esmry ( $dbase, $web, $key, $num, $id, $mode, $min, $max, $tool, $email,
            $silent, $verbose, $debug, $log, $http, $alias, $basx );

    return;
  }

  if ( $dbase eq "structure" and $type eq "mmdb" ) {

    emmdb ( $dbase, $web, $key, $num, $id, $tool, $email );

    return;
  }

  if ( $dbase ne "" and ( $type eq "UID" or $type eq "uid" ) ) {

    if ( $id ne "" ) {

      my @ids = split (',', $id);
      foreach $uid (@ids) {
        $uid =~ s/(\d+)\.\d+/$1/g;
        print "$uid\n";
      }

      return;
    }

    if ( $web eq "" ) {
      die "WebEnv value not found in fetch input\n";
    }

    if ( $pipe ) {
      write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
    }

    if ( $max == 0 ) {
      if ( $silent ) {
        return;
      }
    }

    # use larger chunk for UID format
    $chunk = 5000;
    for ( $start = $min; $start < $max; $start += $chunk ) {

      my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $max, $tool, $email );
      foreach $uid (@ids) {
        print "$uid\n";
      }

      do_sleep ();
    }

    return;
  }

  if ( $dbase ne "" and ( $type eq "URL" or $type eq "url" ) ) {

    if ( $id ne "" ) {

      my @ids = split (',', $id);
      $url = "https://www.ncbi.nlm.nih.gov/";
      $url .= "$dbase/";
      my $pfx = "";
      foreach $uid (@ids) {
        $uid =~ s/(\d+)\.\d+/$1/g;
        $url .= "$pfx$uid";
        $pfx = ",";
      }
      print "$url\n";

      return;
    }

    if ( $web eq "" ) {
      die "WebEnv value not found in fetch input\n";
    }

    if ( $pipe ) {
      write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
    }

    if ( $max == 0 ) {
      if ( $silent ) {
        return;
      }
    }

    # use larger chunk for URL format
    $chunk = 2000;
    for ( $start = $min; $start < $max; $start += $chunk ) {

      my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $max, $tool, $email );
      $url = "https://www.ncbi.nlm.nih.gov/";
      $url .= "$dbase/";
      my $pfx = "";
      foreach $uid (@ids) {
        $url .= "$pfx$uid";
        $pfx = ",";
      }
      print "$url\n";

      do_sleep ();
    }

    return;
  }

  if ( $dbase ne "" and ( $type eq "URLS" or $type eq "urls" ) ) {

    if ( $id ne "" ) {

      my @ids = split (',', $id);
      foreach $uid (@ids) {
        $uid =~ s/(\d+)\.\d+/$1/g;
        print "https://www.ncbi.nlm.nih.gov/$dbase/$uid\n";
      }

      return;
    }

    if ( $web eq "" ) {
      die "WebEnv value not found in fetch input\n";
    }

    if ( $pipe ) {
      write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
    }

    if ( $max == 0 ) {
      if ( $silent ) {
        return;
      }
    }

    # use larger chunk for URL format
    $chunk = 2000;
    for ( $start = $min; $start < $max; $start += $chunk ) {

      my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $max, $tool, $email );
      foreach $uid (@ids) {
        print "https://www.ncbi.nlm.nih.gov/$dbase/$uid\n";
      }

      do_sleep ();
    }

    return;
  }

  if ( $dbase ne "" and $id ne "" ) {

    if ( $id eq "0" ) {

      # id "0" returns a live UID for any database

      $id = get_zero_uid ($dbase);

      if ( $id eq "0" ) {

        # id "0" is an unrecognized accession

        return;
      }
    }

    $url = $base . $efetch;

    if ( $dbase eq "pubmed" ) {
      $id =~ s/(\d+)\.\d+/$1/g;
    }

    $arg = "db=$dbase&id=$id";

    if ( $type eq "gb" ) {
      if ( $style eq "withparts" or $style eq "master" ) {
        $arg .= "&rettype=gbwithparts";
        $arg .= "&retmode=$mode";
      } elsif ( $style eq "conwithfeat" or $style eq "withfeat" or $style eq "contigwithfeat" ) {
        $arg .= "&rettype=$type";
        $arg .= "&retmode=$mode";
        $arg .= "&gbconwithfeat=1";
      } else {
        $arg .= "&rettype=$type";
        $arg .= "&retmode=$mode";
      }
    } else {
      $arg .= "&rettype=$type";
      $arg .= "&retmode=$mode";
    }

    # -chr_start and -chr_stop are for 0-based sequence coordinates from EntrezGene
    if ( $chr_start > -1 && $chr_stop > -1 ) {
      $seq_start = $chr_start + 1;
      $seq_stop = $chr_stop + 1;
    }

    # if -seq_start > -seq_stop, swap values to normalize, indicate minus strand with -strand 2
    if ( $seq_start > 0 && $seq_stop > 0 ) {
      if ( $seq_start > $seq_stop ) {
        my $tmp = $seq_start;
        $seq_start = $seq_stop;
        $seq_stop = $tmp;
        $strand = "2";
      }
    }

    # option to show GI number (undocumented)
    if ( $showgi ) {
      $arg .= "&showgi=1";
    }

    # optionally extend retrieved sequence range in both directions
    if ( $extend > 0 ) {
      $seq_start -= $extend;
      $seq_stop += $extend;
    }

    if ( $strand ne "" ) {
      $arg .= "&strand=$strand";
    }
    if ( $seq_start > 0 ) {
      $arg .= "&seq_start=$seq_start";
    }
    if ( $seq_stop > 0 ) {
      $arg .= "&seq_stop=$seq_stop";
    }
    if ( $complexity > 0 ) {
      $arg .= "&complexity=$complexity";
    }
    if ( $extrafeat > -1 ) {
      $arg .= "&extrafeat=$extrafeat";
    }

    $data = do_post_yielding_ref ($url, $arg, $tool, $email, true);

    if ($$data =~ /<eFetchResult>/i and $$data =~ /<ERROR>(.+?)<\/ERROR>/i) {
      $err = $1;
      if ( $err ne "" ) {
        if ( ! $silent ) {
          print STDERR "ERROR in efetch: $err\n";
        }
      }
    }

    Encode::_utf8_on($$data);

    if (! $raw) {

      if ( $dbase eq "sra" and $type eq "full" and $mode eq "xml" ) {
        $$data = fix_sra_xml_encoding($$data);
      }

      if ( $dbase eq "pubmed" and $type eq "full" and $mode eq "xml" ) {
        $$data = fix_pubmed_xml_encoding($$data);
      }

      if ( $type eq "fasta" or $type eq "fasta_cds_aa" or $type eq "fasta_cds_na" or $type eq "gene_fasta" ) {
        # remove blank lines in FASTA format
        $$data =~ s/\n+/\n/g;
      }
    }

    print $$data;

    return;
  }

  if ( $dbase ne "" and $web ne "" and $key eq "" and $num eq "0" ) {
    die "QueryKey value not found in fetch input\n";
    return;
  }

  test_edirect ( $dbase, $web, $key, $num, "fetch" );

  $stpminusone = $stp - 1;

  if ( $max == 0 ) {
    if ( $silent ) {
      return;
    }
  }

  # use small chunk because fetched records could be quite large
  $chunk = 100;

  # use larger chunk for accessions
  if ( $dbase eq "nucleotide" or
       $dbase eq "nuccore" or
       $dbase eq "est" or
       $dbase eq "gss" or
       $dbase eq "protein" ) {
    if ( $type eq "ACCN" or $type eq "accn" or $type eq "ACC" or $type eq "acc" ) {
      $chunk = 4000;
    }
  }

  for ( $start = $min; $start < $max; $start += $chunk ) {
    $url = $base . $efetch;

    $chkx = $chunk;
    if ( $start + $chkx > $max ) {
      $chkx = $max - $start;
    }

    $arg = "db=$dbase&query_key=$key&WebEnv=$web";

    if ( $type eq "gb" ) {
      if ( $style eq "withparts" or $style eq "master" ) {
        $arg .= "&rettype=gbwithparts";
        $arg .= "&retmode=$mode";
      } elsif ( $style eq "conwithfeat" or $style eq "withfeat" or $style eq "contigwithfeat" ) {
        $arg .= "&rettype=$type";
        $arg .= "&retmode=$mode";
        $arg .= "&gbconwithfeat=1";
      } else {
        $arg .= "&rettype=$type";
        $arg .= "&retmode=$mode";
      }
    } else {
      $arg .= "&rettype=$type";
      $arg .= "&retmode=$mode";
    }

    $arg .= "&retstart=$start&retmax=$chkx";
    if ( $strand ne "" ) {
      $arg .= "&strand=$strand";
    }
    if ( $seq_start > 0 ) {
      $arg .= "&seq_start=$seq_start";
    }
    if ( $seq_stop > 0 ) {
      $arg .= "&seq_stop=$seq_stop";
    }
    if ( $complexity > 0 ) {
      $arg .= "&complexity=$complexity";
    }
    if ( $extrafeat > -1 ) {
      $arg .= "&extrafeat=$extrafeat";
    }

    $data = "";
    $retry = true;

    for ( $tries = 0; $tries < 3 && $retry; $tries++) {
      $data = do_post_yielding_ref ($url, $arg, $tool, $email, true);

      if ($$data =~ /<eFetchResult>/i and $$data =~ /<ERROR>(.+?)<\/ERROR>/i) {
        if ( ! $silent ) {
          print STDERR "Retrying efetch, step $stp: $err\n";
        }
        sleep 3;
      } else {
        $retry = false;
      }
    }

    if ($$data =~ /<eFetchResult>/i and $$data =~ /<ERROR>(.+?)<\/ERROR>/i) {
      $err = $1;
      if ( $err ne "" ) {
        if ( ! $silent ) {
          my $from = $start + 1;
          my $to = $start + $chunk;
          if ( $to > $num ) {
            $to = $num;
          }
          print STDERR "ERROR in efetch ($from-$to / $num): $err\n";
          print STDERR "Replicate for debugging with:\n";
          print STDERR "  edirutil -db $dbase -web $web -key $key -count $num";
          if ( $stpminusone > 0 ) {
            print STDERR " -step $stpminusone";
          }
          my $seconds = "";
          my $end_time = Time::HiRes::time();
          my $elapsed = $end_time - $begin_time;
          if ( $elapsed > 0.0005 ) {
            if ( $elapsed =~ /(\d+)\.(\d\d\d)\d+/ ) {
              $seconds = "$1.$2";
            }
          }
          if ( $seconds ne "" ) {
            print STDERR " -seconds $seconds";
          }
          print STDERR " | efetch -format $type";
          if ( $mode ne "text" ) {
            print STDERR " -mode $mode";
          }
          print STDERR " -start $from -stop $to\n";
        }
      }
    } else {
      if ( $verbose ) {
        my $from = $start + 1;
        my $to = $start + $chunk;
        if ( $to > $num ) {
          $to = $num;
        }
        print STDERR "( edirutil -db $dbase -web $web -key $key -count $num";
        if ( $stpminusone > 0 ) {
          print STDERR " -step $stpminusone";
        }
        my $seconds = "";
        my $end_time = Time::HiRes::time();
        my $elapsed = $end_time - $begin_time;
        if ( $elapsed > 0.0005 ) {
          if ( $elapsed =~ /(\d+)\.(\d\d\d)\d+/ ) {
            $seconds = "$1.$2";
          }
        }
        if ( $seconds ne "" ) {
          print STDERR " -seconds $seconds";
        }
        print STDERR " | efetch -format $type";
        if ( $mode ne "text" ) {
          print STDERR " -mode $mode";
        }
        print STDERR " -start $from -stop $to )\n";
        $begin_time = $end_time;
      }

      Encode::_utf8_on($$data);

      if (! $raw) {
        if ( $dbase eq "sra" and $type eq "full" and $mode eq "xml" ) {
          $$data = fix_sra_xml_encoding($$data);
        }

        if ( $dbase eq "pubmed" and $type eq "full" and $mode eq "xml" ) {
          $$data = fix_pubmed_xml_encoding($$data);
        }

        if ( $type eq "fasta" or $type eq "fasta_cds_aa" or $type eq "fasta_cds_na" or $type eq "gene_fasta" ) {
          # remove blank lines in FASTA format
          $$data =~ s/\n+/\n/g;
        }
      }

      print $$data;
    }

    do_sleep ();
  }
}

# einfo obtains names of databases or names of fields and links per database

my $info_help = qq{
Database Selection

  -db        Database name
  -dbs       Get all database names

Data Summaries

  -fields    Print field names
  -links     Print link names

Field Example

  <Field>
    <Name>ALL</Name>
    <FullName>All Fields</FullName>
    <Description>All terms from all searchable fields</Description>
    <TermCount>138982028</TermCount>
    <IsDate>N</IsDate>
    <IsNumerical>N</IsNumerical>
    <SingleToken>N</SingleToken>
    <Hierarchy>N</Hierarchy>
    <IsHidden>N</IsHidden>
    <IsTruncatable>Y</IsTruncatable>
    <IsRangable>N</IsRangable>
  </Field>

Link Example

  <Link>
    <Name>pubmed_protein</Name>
    <Menu>Protein Links</Menu>
    <Description>Published protein sequences</Description>
    <DbTo>protein</DbTo>
  </Link>
  <Link>
    <Name>pubmed_protein_refseq</Name>
    <Menu>Protein (RefSeq) Links</Menu>
    <Description>Link to Protein RefSeqs</Description>
    <DbTo>protein</DbTo>
  </Link>

};

sub einfo {

  # ... | edirect.pl -info -db pubmed | ...

  clearflags ();

  MyGetOptions(
    $info_help,
    "db=s" => \$db,
    "dbs" => \$dbs,
    "fields" => \$fields,
    "links" => \$links,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "einfo $version\n";
    print $info_help;
    return;
  }

  read_aliases ();
  adjust_base ();

  if ( @ARGV  ||  ($^O =~ /^MSWin/ && ! -t STDIN) ) {
    while ( defined($thisline = <STDIN>) ) {
      $tool = $1 if ( $thisline =~ /<Tool>(.+?)<\/Tool>/i );
      $email = $1 if ( $thisline =~ /<Email>(.+?)<\/Email>/i );
    }
  }

  if ( $dbase eq "" ) {
    $dbase = $db;
  }

  $dbase = lc($dbase);

  if ( $dbase eq "" and (! $dbs) ) {
    die "Must supply -db or -dbs on command line\n";
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  binmode STDOUT, ":utf8";

  $url = $base . $einfo;

  $prefix = "?";

  if ( $dbase ne "" ) {
    $url .= "$prefix" . "db=$dbase&version=2.0";
    $prefix = "&";
  }

  if ( $os ne "" ) {
    $url .= "$prefix" . "edirect_os=$os";
    $prefix = "&";
  }

  if ( $api_key ne "" ) {
    $url .= "$prefix" . "api_key=$api_key";
    $prefix = "&";
  }

  $url .= "$prefix" . "edirect=$version";
  $prefix = "&";

  if ( $tool eq "" ) {
    $tool = "edirect";
  }
  if ( $tool ne "" ) {
    $url .= "$prefix" . "tool=$tool";
    $prefix = "&";
  }

  if ( $email eq "" ) {
    $email = get_email ();
  }
  if ( $email ne "" ) {
    $url .= "$prefix" . "email=$email";
    $prefix = "&";
  }

  if ( $debug or $log ) {
    print STDERR "$url\n";
  }

  $output = get ($url);

  if ( ! defined $output ) {
    print STDERR "Failure of '$url'\n";
    return;
  }

  if ( $output eq "" ) {
    print STDERR "No einfo output returned from '$url'\n";
    return;
  }

  if ($output =~ /IsTrunc/ or $output =~ /IsRang/) {
    # print STDERR "New server is out with IsTrunc/IsRang - update executable\n";
  }

  if ( $debug ) {
    print STDERR "$output\n";
  }

  if ( $dbs and $output =~ /<DbName>/ ) {

    # -dbs now extracts database names from XML

    $output =~ s/\r//g;
    $output =~ s/\n//g;
    $output =~ s/\t//g;
    $output =~ s/ +/ /g;
    $output =~ s/> +</></g;

    my @databases = ($output =~ /<DbName>(.+?)<\/DbName>/g);
    foreach $dtbs (@databases) {
      print "$dtbs\n";
    }

    return;
  }

  if ( ( $fields or $links ) and $output =~ /<DbInfo>/ ) {

    # -db can print information directly without need to process XML with xtract

    $output =~ s/\r//g;
    $output =~ s/\n//g;
    $output =~ s/\t//g;
    $output =~ s/ +/ /g;
    $output =~ s/> +</></g;

    my $name = "";
    my $full = "";
    my $menu = "";

    if ( $fields ) {
      my @flds = ($output =~ /<Field>(.+?)<\/Field>/g);
      foreach $fld (@flds) {
        $name = "";
        $full = "_";
        if ( $fld =~ /<Name>(.+?)<\/Name>/ ) {
          $name = $1;
        }
        if ( $fld =~ /<FullName>(.+?)<\/FullName>/ ) {
          $full = $1;
        }
        if ( $name ne "" and $full ne "" ) {
          print "$name\t$full\n";
        }
      }
    }

    if ( $links ) {
      my @lnks = ($output =~ /<Link>(.+?)<\/Link>/g);
      foreach $lnk (@lnks) {
        $name = "";
        $menu = "_";
        if ( $lnk =~ /<Name>(.+?)<\/Name>/ ) {
          $name = $1;
        }
        if ( $lnk =~ /<Menu>(.+?)<\/Menu>/ ) {
          $menu = $1;
        }
        if ( $name ne "" and $menu ne "" ) {
          print "$name\t$menu\n";
        }
      }
    }

    return;
  }

  print "$output";
}

# common link history processing

sub acheck_test {

  my $dbase = shift (@_);
  my $dbto = shift (@_);
  my $name = shift (@_);
  my $web = shift (@_);
  my $key = shift (@_);
  my $num = shift (@_);
  my $stp = shift (@_);
  my $err = shift (@_);
  my $tool = shift (@_);
  my $email = shift (@_);

  if ( $num > 0 ) {
    return;
  }

  # if failure, use acheck to confirm that there are no links

  $url = $base . $elink;
  $arg = "dbfrom=$dbase&query_key=$key&WebEnv=$web&cmd=acheck";
  $data = do_post ($url, $arg, $tool, $email, true);

  # remove newlines, tabs, space between tokens, compress runs of spaces,

  $data =~ s/\r//g;
  $data =~ s/\n//g;
  $data =~ s/\t//g;
  $data =~ s/ +/ /g;
  $data =~ s/> </></g;

  if ($data =~ /<Error>(.+?)<\/Error>/i) {
    $err = $1;
    if ( $err ne "" ) {
      write_edirect ( $dbto, $web, $key, $num, $stp, $err, $tool, $email );
      close (STDOUT);
      die "ERROR in acheck test: $err\n";
    }
  }

  if ( $data !~ /<LinkInfo>/ ) {
    write_edirect ( $dbto, $web, $key, $num, $stp, $err, $tool, $email );
    close (STDOUT);
    die "Elink acheck confirmation test failed\n";
  }

  if ( $data =~ /<LinkName>$name<\/LinkName>/ ) {
    write_edirect ( $dbto, $web, $key, $num, $stp, $err, $tool, $email );
    close (STDOUT);
    die "Elink acheck test indicates non-zero count expected\n";
  }
}

sub process_history_link {

  my $arg = shift (@_);
  my $output = shift (@_);
  my $dbase = shift (@_);
  my $dbto = shift (@_);
  my $name = shift (@_);
  my $wb = shift (@_);
  my $web = shift (@_);
  my $key = shift (@_);
  my $stp = shift (@_);
  my $err = shift (@_);
  my $lbl = shift (@_);
  my $tool = shift (@_);
  my $email = shift (@_);

  my $num = "";

  if ( $wb eq "" ) {
    $wb = $web;
  }

  if ( $err ne "" ) {
    write_edirect ( "", $wb, "", "", "", $err, "", "" );
    close (STDOUT);
    if ( ! $silent ) {
      die "ERROR in link output: $err\nWebEnv: $wb\nURL: $arg\nResult: $output\n\n";
    }
    return;
  }

  if ( $web eq "" ) {
    die "WebEnv value not found in link output - WebEnv1 $wb\n";
  }
  if ( $key eq "" ) {
    write_edirect ( $dbto, $wb, $key, "0", $stp, $err, $tool, $email );
    close (STDOUT);
    # no neighbors or links can be a normal response,
    # e.g., elink -db gene -id 496376 -target medgen
    # so suppress this message
    # die "QueryKey value not found in link output - WebEnv1 $wb\n";
    return;
  }

  if ( $web ne $wb ) {
    $err = "WebEnv mismatch in link-search output - WebEnv1 $wb, WebEnv2 $web";
    write_edirect ( "", $web, "", "", "", $err, "", "" );
    close (STDOUT);
    die "WebEnv value changed in link-search output - WebEnv1 $wb\nWebEnv2 $web\n";
  }

  ( $num, $key ) = get_count ( $dbto, $web, $key, $tool, $email );

  if ( $num eq "" ) {
    $err = "Missing count in link-search output - WebEnv $web";
    write_edirect ( "", $web, "", "", "", $err, "", "" );
    close (STDOUT);
    die "Count value not found in link-search output - WebEnv $web\n";
  }

  if ( $num == 0 ) {
    acheck_test ( $dbase, $dbto, $name, $web, $key, $num, $stp, $err, $tool, $email );
  }

  if ( $lbl ne "" and $key ne "" ) {
    $labels{"$lbl"} = "$key";
  }

  write_edirect ( $dbto, $web, $key, $num, $stp, $err, $tool, $email );
}

# for large lists, break into chunks, merge result on client, post to server

sub batch_elink {

  my $dbase = shift (@_);
  my $dbto = shift (@_);
  my $name = shift (@_);
  my $web = shift (@_);
  my $key = shift (@_);
  my $num = shift (@_);
  my $stp = shift (@_);
  my $err = shift (@_);
  my $lbl = shift (@_);
  my $tool = shift (@_);
  my $email = shift (@_);
  my $auto = shift (@_);

  my %seen = ();
  my @uniq = ();

  if ( $num == 0 ) {
    if ( $silent ) {
      write_edirect ( "", "", "", "", "", $err, "", "" );
      close (STDOUT);
      return;
    }
  }

  $chunk = 200;
  for ( $start = 0; $start < $num; $start += $chunk ) {

    my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $num, $tool, $email );

    $url = $base . $elink;

    $arg = "dbfrom=$dbase&db=$dbto&cmd=neighbor&linkname=$name";
    if ( $mode ne "" ) {
      $arg .= "&retmode=$mode";
    }
    $arg .= "&id=";
    $arg .= join (',', @ids);

    $data = "";
    $retry = true;

    for ( $tries = 0; $tries < 3 && $retry; $tries++) {
      $data = do_post ($url, $arg, $tool, $email, true);

      if ($data =~ /<Error>(.+?)<\/Error>/i) {
        if ( ! $silent ) {
          print STDERR "Retrying batch elink, step $stp: $err\n";
        }
        sleep 3;
      } else {
        $retry = false;
      }
    }

    # remove newlines, tabs, space between tokens, compress runs of spaces,

    $data =~ s/\r//g;
    $data =~ s/\n//g;
    $data =~ s/\t//g;
    $data =~ s/ +/ /g;
    $data =~ s/> </></g;

    if ($data =~ /<Error>(.+?)<\/Error>/i) {
      $err = $1;
      if ( $err ne "" ) {
        if ( ! $silent ) {
          my $from = $start + 1;
          my $to = $start + $chunk;
          if ( $to > $num ) {
            $to = $num;
          }
          print STDERR "ERROR in batch elink ($from-$to / $num): $err\n";
        }
      }
    }

    if ($data !~ /<LinkSetDb>/i) {
      if ( ! $silent ) {
        print STDERR "LinkSetDb missing in batch elink result\n";
      }
    }

    while ( $data =~ /<LinkSetDb>(.*?)<\/LinkSetDb>/g ) {
      $linkset = $1;
      while ( $linkset =~ /<Id>(\d+)<\/Id>/g ) {
        $uid = $1;
        if ( ! $seen{$uid}++ ) {
          push (@uniq, $uid);
        }
      }
    }

    do_sleep ();
  }

  $dbase = $dbto;

  $dbase = lc($dbase);

  $num = scalar @uniq;

  if ( $num == 0 ) {
    acheck_test ( $dbase, $dbto, $name, $web, $key, $num, $stp, $err, $tool, $email );
    write_edirect ( $dbase, $web, $key, "0", $stp, $err, $tool, $email );
    if ( $auto ) {
      close (STDOUT);
      die "Automatic batch link failed\n";
    }
    return;
  }

  # not certain if sort is necessary - need to experiment before final release

  @sorted = sort { $a <=> $b } @uniq;

  $url = $base . $epost;

  $arg = "db=$dbase";
  if ( $key ne "" ) {
    $arg .= "&query_key=$key";
  }
  if ( $web ne "" ) {
    $arg .= "&WebEnv=$web";
  }
  $ids = join (',', @sorted);

  $arg .= "&id=$ids";

  $output = do_post ($url, $arg, $tool, $email, true);

  if ( $debug ) {
    print STDERR "$output\n";
  }

  $wb = $web;

  $web = "";
  $key = "";
  $err = "";

  $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
  $key = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
  $err = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);

  process_history_link ( $arg, $output, $dbase, $dbto, $name, $wb, $web, $key, $stp, $err, $lbl, $tool, $email );
}

# elink without a target uses the source database for neighboring

my $link_help = qq{
Destination Database

  -related    Neighbors in same database
  -target     Links in different database
  -name       Link name (e.g., pubmed_protein_refseq)

Direct Record Selection

  -db         Database name
  -id         Unique identifier(s)

Advanced Control

  -cmd        Command type (returns eLinkResult XML)
  -mode       "ref" uses LinkOut provider's web site
  -holding    Name of LinkOut provider

Batch Processing

  -batch      Bypass Entrez history mechanism

Miscellaneous Arguments

  -label      Alias for query step

Command Option Examples

  -cmd              Result
  ____              ______

  neighbor          Neighbors or links

  neighbor_score    Neighbors with computed similarity scores

  acheck            All links available

  ncheck            Existence of neighbors

  lcheck            Existence of external links (LinkOuts)

  llinks            Non-library LinkOut providers

  llinkslib         All LinkOut providers

  prlinks           Primary LinkOut provider,
                    or URL for single UID with -mode ref

};

sub elink {

  # ... | edirect.pl -link -target structure | ...

  clearflags ();

  MyGetOptions(
    $link_help,
    "db=s" => \$db,
    "id=s" => \$id,
    "format=s" => \$type,
    "target=s" => \$dbto,
    "name=s" => \$name,
    "related" => \$related,
    "neighbor" => \$neighbor,
    "cmd=s" => \$cmd,
    "mode=s" => \$mode,
    "batch" => \$batch,
    "holding=s" => \$holding,
    "label=s" => \$lbl,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "elink $version\n";
    print $link_help;
    return;
  }

  # convert spaces between UIDs to commas

  $id =~ s/ /,/g;
  $id =~ s/,,/,/g;

  if ( -t STDIN and not @ARGV ) {
  } elsif ( $db ne "" and $id ne "" ) {
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    write_edirect ( "", "", "", "", "", $err, "", "" );
    close (STDOUT);
    if ( ! $silent ) {
      die "ERROR in link input: $err\n\n";
    }
    return;
  }

  if ( $dbase eq "" ) {
    $dbase = $db;
  }

  $dbase = lc($dbase);

  my $adddbto = true;
  if ( $dbto eq "" ) {
    if ( $cmd eq "acheck" or
         $cmd eq "ncheck" or
         $cmd eq "lcheck" or
         $cmd eq "llinks" or
         $cmd eq "llinkslib" or
         $cmd eq "prlinks" ) {
      $dbto = $dbase;
      $adddbto = false;
    }
  }

  if ( $dbto eq "" and (! $related) and (! $neighbor) ) {
    die "Must supply -target or -related on command line\n";
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  binmode STDOUT, ":utf8";

  if ( $dbto eq "" ) {
    $dbto = $dbase;
  }
  if ( $name eq "" ) {
    $name = $dbase . "_" . $dbto;
  }

  if ( $cmd eq "" ) {
    $cmd = "neighbor_history";
  }

  if ( $dbase eq "nlmcatalog" ) {
    die "Entrez Direct does not support links for the nlmcatalog database\n";
  }

  if ( $dbase ne "" and $id ne "" ) {

    # process db and id command-line arguments instead of getting from history

    $url = $base . $elink;

    $arg = "dbfrom=$dbase";
    if ( $adddbto ) {
      $arg .= "&db=$dbto";
    }
    $arg .= "&cmd=$cmd&linkname=$name";
    if ( $mode ne "" ) {
      $arg .= "&retmode=$mode";
    }
    if ( $dbase eq "pubmed" ) {
      $id =~ s/(\d+)\.\d+/$1/g;
    }
    $arg .= "&id=$id";
    if ( $type eq "acc" ) {
      $arg .= "&idtype=acc";
    }

    $data = do_post ($url, $arg, $tool, $email, true);

    if ( $cmd ne "neighbor_history" ) {

      # if not neighbor_history, write eLinkResult XML instead of ENTREZ_DIRECT

      print "$data";

      return;
    }

    $wb = $web;

    $web = "";
    $key = "";
    $err = "";

    $err = $1 if ($data =~ /<Error>(.+?)<\/Error>/i);
    if ( $err eq "" ) {
      $web = $1 if ($data =~ /<WebEnv>(\S+)<\/WebEnv>/);
      $key = $1 if ($data =~ /<QueryKey>(\S+)<\/QueryKey>/);
    }

    process_history_link ( $arg, $data, $dbase, $dbto, $name, $wb, $web, $key, $stp, $err, $lbl, $tool, $email );

    return;
  }

  test_edirect ( $dbase, $web, $key, $num, "link" );

  if ( $cmd ne "neighbor_history" ) {

    # if not neighbor_history, write eLinkResult XML instead of ENTREZ_DIRECT

    if ( $num == 0 ) {
      if ( $silent ) {
        return;
      }
    }

    $chunk = 200;
    for ( $start = 0; $start < $num; $start += $chunk ) {

      my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $num, $tool, $email );

      $url = $base . $elink;

      $arg = "dbfrom=$dbase";
      if ( $adddbto ) {
        $arg .= "&db=$dbto";
      }
      $arg .= "&cmd=$cmd&linkname=$name";
      if ( $mode ne "" ) {
        $arg .= "&retmode=$mode";
      }
      $arg .= "&id=";
      $arg .= join ('&id=', @ids);
      if ( $type eq "acc" ) {
        $arg .= "&idtype=acc";
      }

      $data = do_post ($url, $arg, $tool, $email, true);

      print "$data";

      do_sleep ();
    }

    return;
  }

  if ( $batch ) {

    # large list bypass

    batch_elink ( $dbase, $dbto, $name, $web, $key, $num, $stp, "", $lbl, $tool, $email, false );
    return;
  }

  # if not breaking large lists into chunks, use neighbor_history on server

  $url = $base . $elink;

  $arg = "dbfrom=$dbase";
  if ( $adddbto ) {
    $arg .= "&db=$dbto";
  }
  $arg .= "&query_key=$key&WebEnv=$web";
  $arg .= "&cmd=$cmd&linkname=$name";
  if ( $mode ne "" ) {
    $arg .= "&retmode=$mode";
  }
  if ( $type eq "acc" ) {
    $arg .= "&idtype=acc";
  }

  $wb = $web;
  $ky = $key;
  $nm = $num;

  $web = "";
  $key = "";
  $err = "";

  $output = "";

  for ( $tries = 0; $tries < 3 && $web eq ""; $tries++) {
    $output = do_post ($url, $arg, $tool, $email, true);

    $err = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);
    if ( $err eq "" ) {
      $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
      $key = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
    } else {
      if ( ! $silent ) {
        print STDERR "Retrying elink, step $stp: $err\n";
      }
      sleep 3;
    }
  }

  # automatically fail over to batch mode under certain failure conditions

  $stpminusone = $stp - 1;

  if ( $err =~ "^Query failed" or
       $err =~ "^Timeout waiting" or
       $err =~ "^Unable to obtain" or
       $err =~ "^The read request has timed out" ) {
    if ( ! $silent ) {
      print STDERR "$err\nReplicate for debugging with:\n";
      print STDERR "  edirutil -db $dbase -web $wb -key $ky -count $nm";
      if ( $stpminusone > 0 ) {
        print STDERR " -step $stpminusone";
      }
      my $seconds = "";
      my $end_time = Time::HiRes::time();
      my $elapsed = $end_time - $begin_time;
      if ( $elapsed > 0.0005 ) {
        if ( $elapsed =~ /(\d+)\.(\d\d\d)\d+/ ) {
          $seconds = "$1.$2";
        }
      }
      if ( $seconds ne "" ) {
        print STDERR " -seconds $seconds";
      }
      if ( $dbase eq $dbto ) {
        print STDERR " | elink -related\n";
      } else {
        print STDERR " | elink -target $dbto\n";
      }
      print STDERR "Automatically switching to -batch mode\n";
    }
    batch_elink ( $dbase, $dbto, $name, $wb, $ky, $nm, $stp, "", $lbl, $tool, $email, true );
    return;
  }

  # finish reality checks, get count and key, and write ENTREZ_DIRECT structure

  process_history_link ( $arg, $output, $dbase, $dbto, $name, $wb, $web, $key, $stp, $err, $lbl, $tool, $email );
}

# emmdb downloads "Ncbi-mime-asn1 ::= strucseq" ASN.1 from MMDB unique identifiers

sub emmdb {

  # ... | edirect.pl -fetch -format mmdb | ...

  $dbase = shift (@_);
  $web = shift (@_);
  $key = shift (@_);
  $num = shift (@_);
  $id = shift (@_);
  $tool = shift (@_);
  $email = shift (@_);

  if ( $err ne "" ) {
    die "ERROR in mmdb input: $err\n\n";
  }

  $mbase = "https://www.ncbi.nlm.nih.gov/Structure/mmdb/";
  $mprog = "mmdbsrv.cgi";

  if ( $id ne "" ) {

    my @ids = split (',', $id);
    foreach $uid (@ids) {
      $uid =~ s/(\d+)\.\d+/$1/g;
      $url = $mbase . $mprog;

      $arg = "uid=$uid";
      $arg .= "&save=asntext&form=6&db=t&Dopt=j&Complexity=Cn3D%20Subset";

      $mmdb = do_post ($url, $arg, $tool, $email, true);

      print "$mmdb\n";
    }

    return;
  }

  test_edirect ( $dbase, $web, $key, $num, "mmdb" );

  $chunk = 100;
  for ( $start = 0; $start < $num; $start += $chunk ) {

    my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $num, $tool, $email );

    foreach $uid (@ids) {
      $url = $mbase . $mprog;

      $arg = "uid=$uid";
      $arg .= "&save=asntext&form=6&db=t&Dopt=j&Complexity=Cn3D%20Subset";

      $mmdb = do_post ($url, $arg, $tool, $email, true);

      print "$mmdb\n";

      do_sleep ();
    }
  }
}

# entfy sends e-mail

my $ntfy_help = qq{
  -email    Contact person's address
  -tool     Name of script or program

};

sub entfy {

  # ... | edirect.pl -notify

  clearflags ();

  MyGetOptions(
    $ntfy_help,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "enotify $version\n";
    print $ntfy_help;
    return;
  }

  ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    die "ERROR in notify input: $err\n\n";
  }

  test_edirect ( $dbase, $web, $key, $num, "notify" );

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }
  if ( $email eq "" ) {
    $email = get_email ();
  }
  if ( $email eq "" ) {
    die "Email value not found in notify input\n";
  }

  binmode STDOUT, ":utf8";

  if ( $num > 0 ) {
    $chunk = 100;
    for ( $start = 0; $start < $num; $start += $chunk ) {

      my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $num, $tool, $email );

      foreach $uid (@ids) {
        $txt = "echo \"https://www.ncbi.nlm.nih.gov/$dbase/$uid\n\"";
        $str = "mail -s \"A new $dbase record is in Entrez\" $email";
        system "$txt | $str";
      }

      do_sleep ();
    }
  }
}

# epost uploads UIDs or accessions

sub post_chunk {

  my $dbsx = shift (@_);
  my $webx = shift (@_);
  my $keyx = shift (@_);
  my $tulx = shift (@_);
  my $emlx = shift (@_);
  my $uids = shift (@_);
  my $qryx = shift (@_);

  $url = $base;
  $arg = "";
  $wb = $webx;

  if ( $uids ne "" ) {

    $url .= $epost;

    $arg = "db=$dbsx";
    if ( $web ne "" ) {
      $arg .= "&WebEnv=$webx";
    }
    $arg .= "&id=$uids";

  } elsif ( $qryx ne "" ) {

    $url .= $esearch;

    $qryx = map_labels ($qryx);
    $qryx = map_macros ($qryx);
    $enc = uri_escape($qryx);

    $arg = "db=$dbsx&term=$enc";
    if ( $web ne "" ) {
      $arg .= "&WebEnv=$webx";
    }
    $arg .= "&retmax=0&usehistory=y";
    if ( $rldate > 0 ) {
      $arg .= "&reldate=$rldate";
      if ( $dttype eq "" ) {
        $dttype = "PDAT";
      }
    }
    if ( $mndate ne "" and $mxdate ne "" ) {
      $arg .= "&mindate=$mndate&maxdate=$mxdate";
      if ( $dttype eq "" ) {
        $dttype = "PDAT";
      }
    }
    if ( $dttype ne "" ) {
      $arg .= "&datetype=$dttype";
    }
  }

  $output = do_post ($url, $arg, $tulx, $emlx, true);

  $webx = "";
  $keyx = "";
  $err = "";

  $webx = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
  $keyx = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
  $err = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);

  if ( $err ne "" ) {
    write_edirect ( "", "", "", "", "", $err, "", "" );
    close (STDOUT);
    die "ERROR in post output: $err\nURL: $arg\n\n";
  }

  if ( $webx eq "" ) {
    die "WebEnv value not found in post output\n";
  }
  if ( $keyx eq "" ) {
    die "QueryKey value not found in post output\n";
  }

  if ( $wb ne "" and $webx ne $wb ) {
    $err = "WebEnv mismatch in post output - WebEnv1 $wb, WebEnv2 $webx";
    write_edirect ( "", $wb, "", "", "", $err, "", "" );
    close (STDOUT);
    die "WebEnv value changed in post output - WebEnv1 $wb\nWebEnv2 $webx\n";
  }

  return $webx, $keyx;
}

my $post_help = qq{
  -db        Database name
  -id        Unique identifier(s) or accession number(s)
  -format    uid or acc
  -input     Read from file instead of stdin
  -label     Alias for query step

};

sub epost {

  # ... | edirect.pl -post -db nucleotide -format uid | ...

  clearflags ();

  MyGetOptions(
    $post_help,
    "db=s" => \$db,
    "id=s" => \$id,
    "format=s" => \$field,
    "input=s" => \$input,
    "label=s" => \$lbl,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "epost $version\n";
    print $post_help;
    return;
  }

  # convert spaces between UIDs to commas

  $id =~ s/ /,/g;
  $id =~ s/,,/,/g;

  if ( -t STDIN and not @ARGV ) {
    if ( $id ne "" ) {
      if ( $id =~ /[,0-9.]/ ) {
        $just_num = true;
      }
    }
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    die "ERROR in post input: $err\n\n";
  }

  if ( $dbase eq "" ) {
    $dbase = $db;
  }

  $dbase = lc($dbase);

  if ( $dbase eq "" ) {
    die "Must supply -db database on command line\n";
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  binmode STDOUT, ":utf8";

  if ( $field eq "" ) {
    $field = "UID";
  }

  # read data from input file instead of piping from stdin
  if ( $input ne "" ) {
    if (open (my $FILE_IN, $input)) {
      $has_num = false;
      $all_num = true;
      while ( $thisline = <$FILE_IN> ) {
        $thisline =~ s/\r//;
        $thisline =~ s/\n//;
        $thisline =~ s/^\s+//;
        $thisline =~ s/\s+$//;

        if ( $thisline =~ /^(\d+)$/ ) {
          push (@rest, $1);
          $has_num = true;
        } elsif ( $thisline =~ /^(.+)$/ ) {
          push (@rest, $1);
          $all_num = false;
        }
      }
      close ($FILE_IN);
      if ( $has_num && $all_num ) {
        $just_num = true;
      }
    } else {
      print STDERR "Unable to open input file '$input'\n";
    }
  }

  my $combo = "";
  my $pfx = "";
  my $loops = 0;

  my $accession_mode = false;
  if ( $field eq "ACCN" or $field eq "accn" or $field eq "ACC" or $field eq "acc" ) {
    $accession_mode = true;
  }

  if ( ! $just_num ) {
    if ( $dbase eq "nucleotide" or
         $dbase eq "nuccore" or
         $dbase eq "est" or
         $dbase eq "gss" or
         $dbase eq "protein" ) {
      $accession_mode = true;
      $field = "acc";
    }
  }

  if ( $field eq "UID" or $field eq "uid" ) {

    if ( $id ne "" ) {

      ( $web, $key ) = post_chunk ( $dbase, $web, $key, $tool, $email, $id, "" );

      $combo .= $pfx . "#" . $key;
      $pfx = " OR ";
      $loops++;

    } else {

      while ( @rest ) {
        my @chunk = splice(@rest, 0, 2000);

        $ids = join (',', @chunk);

        # newline to comma conversion for piped data

        $ids =~ s/\n/,/g;

        if ( $ids =~ /[a-zA-Z]/ ) {
          if ( $dbase eq "nucleotide" or
               $dbase eq "nuccore" or
               $dbase eq "est" or
               $dbase eq "gss" or
               $dbase eq "protein" ) {
            $accession_mode = true;
            $field = "acc";
          } else {
            die "Non-numeric value found in post input:$ids\n\n";
          }
        }

        ( $web, $key ) = post_chunk ( $dbase, $web, $key, $tool, $email, $ids, "" );

        $combo .= $pfx . "#" . $key;
        $pfx = " OR ";
        $loops++;

        do_sleep ();
      }
    }

  } else {

    if ( $id ne "" ) {

      my @chunk = split (',', $id);

      if ( $accession_mode ) {
        $query = join (' [ACCN] OR ', @chunk);
        $query .= " [ACCN]";
        $query =~ s/\./_/g;
      } else {
        $query = join (' OR ', @chunk);
        $query .= " [$field]";
      }

      if ( $accession_mode and $dbase eq "assembly" ) {
        $query =~ s/\[ACCN\]/[ASAC]/g;
      }

      ( $web, $key ) = post_chunk ( $dbase, $web, $key, $tool, $email, "", $query );

      $combo .= $pfx . "#" . $key;
      $pfx = " OR ";
      $loops++;

    } else {

      while ( @rest ) {
        my @chunk = splice(@rest, 0, 2000);

        if ( $accession_mode ) {
          $query = join (' [ACCN] OR ', @chunk);
        } else {
          $query = join (' OR ', @chunk);
        }

        if ( $query eq "" ) {
          die "Must pipe data into stdin\n";
        }

        if ( $accession_mode ) {
        } elsif ( ! $just_num ) {
          die "Non-numeric value found in post input\n";
        }

        if ( $accession_mode ) {
          $query .= " [ACCN]";
          $query =~ s/\./_/g;
        } else {
          $query .= " [$field]";
        }

        if ( $accession_mode and $dbase eq "assembly" ) {
          $query =~ s/\[ACCN\]/[ASAC]/g;
        }

        ( $web, $key ) = post_chunk ( $dbase, $web, $key, $tool, $email, "", $query );

        $combo .= $pfx . "#" . $key;
        $pfx = " OR ";
        $loops++;

        do_sleep ();
      }
    }
  }

  if ( $combo eq "" ) {
    die "Failure of post to find data to load\n";
  }

  if ( $loops > 1 ) {
    $url = $base . $esearch;

    $enc = uri_escape($combo);
    $arg = "db=$dbase&term=$enc&WebEnv=$web";
    $arg .= "&retmax=0&usehistory=y";

    $output = do_post ($url, $arg, $tool, $email, true);

    $web = "";
    $key = "";
    $err = "";

    $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
    $key = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
    $err = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);
  }

  ( $num, $key ) = get_count ( $dbase, $web, $key, $tool, $email );

  if ( $num eq "" ) {
    die "Count value not found in post output - WebEnv1 $web\n";
  }

  if ( $lbl ne "" and $key ne "" ) {
    $labels{"$lbl"} = "$key";
  }

  write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
}

# espel performs an ESpell search

my $spell_help = qq{
  -db       Database name
  -query    Query string

};

sub espel {

  # ... | edirect.pl -spell -db pubmed -query "asthmaa OR alergies" | ...

  clearflags ();

  MyGetOptions(
    $spell_help,
    "db=s" => \$db,
    "query=s" => \$query,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "espell $version\n";
    print $spell_help;
    return;
  }

  # convert spaces between UIDs to commas

  $id =~ s/ /,/g;
  $id =~ s/,,/,/g;

  if ( -t STDIN and not @ARGV ) {
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    die "ERROR in spell input: $err\n\n";
  }

  if ( $dbase eq "" ) {
    $dbase = $db;
  }

  $dbase = lc($dbase);

  if ( $dbase eq "" ) {
    die "Must supply -db database on command line\n";
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  binmode STDOUT, ":utf8";

  if ( $query eq "" ) {
    die "Must supply -query term expression on command line\n";
  }

  $url = $base . $espell;

  $enc = uri_escape($query);
  $arg = "db=$dbase&term=$enc";

  $data = do_post ($url, $arg, $tool, $email, true);

  Encode::_utf8_on($data);

  print "$data";
}

# ecitmtch performs an ECitMatch search

my $citmatch_help = qq{
  -journal    Journal Title
  -year       Year
  -volume     Volume
  -page       First Page
  -author     Author Name

};

sub ecitmtch {

  # ... | edirect.pl -citmatch -journal "proc natl acad sci u s a" -year 2005 ...

  clearflags ();

  MyGetOptions(
    $citmatch_help,
    "journal=s" => \$journal,
    "year=s" => \$year,
    "volume=s" => \$volume,
    "page=s" => \$page,
    "author=s" => \$author,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "ecitmatch $version\n";
    print $citmatch_help;
    return;
  }

  # convert spaces between UIDs to commas

  $id =~ s/ /,/g;
  $id =~ s/,,/,/g;

  if ( -t STDIN and not @ARGV ) {
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    die "ERROR in citation match input: $err\n\n";
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  binmode STDOUT, ":utf8";

  $url = $base . $ecitmat;

  $query = "";

  if ( $journal ne "" ) {
    $query .= $journal;
  }
  $query .= "|";
  if ( $year ne "" ) {
    $query .= $year;
  }
  $query .= "|";
  if ( $volume ne "" ) {
    $query .= $volume;
  }
  $query .= "|";
  if ( $page ne "" ) {
    $query .= $page;
  }
  $query .= "|";
  if ( $author ne "" ) {
    $query .= $author;
  }
  $query .= "||";

  $enc = uri_escape($query);
  $arg = "db=pubmed&retmode=xml&bdata=$enc";

  $data = do_post ($url, $arg, $tool, $email, true);

  Encode::_utf8_on($data);

  if ( $data =~ "NOT_FOUND" ) {
    return;
  }

  if ( $data =~ /.+\|AMBIGUOUS (.+)$/ ) {
    my $my_uids = $1;
    my @ids = split (',', $my_uids);

    foreach $uid (@ids) {
      print "$uid\n";
    }

    return;
  }

  if ( $data =~ /.+\|(\d+)$/ ) {
    my $my_uid = $1;
    print "$my_uid\n";
  }
}

# eprxy reads a file of query proxies, can also pipe from stdin

my $prxy_help = qq{
  -alias    File of aliases
  -pipe     Read aliases from stdin

};

sub eprxy {

  # ... | edirect.pl -proxy -alias file_name | ...

  clearflags ();

  MyGetOptions(
    $prxy_help,
    "pipe" => \$pipe,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "eproxy $version\n";
    print $prxy_help;
    return;
  }

  if ( -t STDIN and not @ARGV ) {
  } elsif ( $pipe ) {
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $pipe ) {
    while ( defined($thisline = <STDIN>) ) {
      $thisline =~ s/\r//;
      $thisline =~ s/\n//;
      $thisline =~ s/ +/ /g;
      $thisline =~ s/> </></g;
      if ( $thisline =~ /(.+)\t(.+)/ ) {
        $ky = $1;
        $vl = $2;
        $vl =~ s/\"//g;
        $macros{"$ky"} = "$vl";
      }
    }
  }

  write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
}

# esrch performs a new EUtils search, but can read a previous web environment value

my $srch_help = qq{
Query Specification

  -db          Database name
  -query       Query string

Document Order

  -sort        Result presentation order

Date Constraint

  -days        Number of days in the past
  -datetype    Date field abbreviation
  -mindate     Start of date range
  -maxdate     End of date range

Limit by Field

  -field       Query words individually in field
  -pairs       Query overlapping word pairs

Spell Check

  -spell       Correct misspellings in query

Miscellaneous Arguments

  -label       Alias for query step

Sort Order Examples

  -db            -sort
  ___            _____

  gene
                 Chromosome
                 Gene Weight
                 Name
                 Relevance

  geoprofiles
                 Default Order
                 Deviation
                 Mean Value
                 Outliers
                 Subgroup Effect

  pubmed
                 First Author
                 Journal
                 Last Author
                 Pub Date
                 Recently Added
                 Relevance
                 Title

  (sequences)
                 Accession
                 Date Modified
                 Date Released
                 Default Order
                 Organism Name
                 Taxonomy ID

  snp
                 Chromosome Base Position
                 Default Order
                 Heterozygosity
                 Organism
                 SNP_ID
                 Success Rate

};

sub remove_punctuation {

  my $qury = shift (@_);

  $qury =~ s/[^a-zA-Z0-9]/ /g;
  $qury =~ s/ +/ /g;

  return $qury;
}

sub remove_stop_words {

  my $qury = shift (@_);

  # split to protect against regular expression artifacts
  $qury =~ s/[^a-zA-Z0-9]/ /g;
  $qury =~ s/ +/ /g;

  my @words = split (' ', $qury);
  my $kept = "";
  my $pfx = "";

  foreach $term (@words) {

    my $trm = lc($term);
    $trm = "#$trm#";

    if ($stop_words !~ /$trm/) {
      $kept .= "$pfx$term";
      $pfx = " ";
    }
  }

  if ( $kept ne "" ) {
    $qury = $kept;
  }

  return $qury;
}

sub field_each_word {

  my $fld = shift (@_);
  my $qury = shift (@_);

  $qury =~ s/,/ /g;
  $qury =~ s/ +/ /g;

  my @words = split (' ', $qury);
  $qury = "";
  my $pfx = "";

  foreach $term (@words) {
    $qury .= "$pfx$term [$fld]";
    $pfx = " AND ";
  }

  return $qury;
}

sub merge_each_word {

  my $fld = shift (@_);
  my $qury = shift (@_);

  $qury =~ s/,/ /g;
  $qury =~ s/ +/ /g;

  my @words = split (' ', $qury);
  $qury = "";
  my $pfx = "";

  foreach $term (@words) {
    $qury .= "$pfx$term [$fld]";
    $pfx = " OR ";
  }

  return $qury;
}

sub field_each_pair {

  my $fld = shift (@_);
  my $qury = shift (@_);

  my @words = split (' ', $qury);
  $qury = "";
  my $pfx = "";

  my $prev = "";
  foreach $term (@words) {
    my $trm = lc($term);
    $trm = "#$trm#";
    if ($stop_words =~ /$trm/) {
      $prev = "";
      $term = "";
    }
    if ( $prev ne "" ) {
      $qury .= "$pfx\"$prev $term\" [$fld]";
      $pfx = " AND ";
    }
    $prev = $term;
  }

  return $qury;
}

sub esrch {

  # ... | edirect.pl -search -db nucleotide -query "M6506* [ACCN] OR 1322283 [UID]" -days 365 | ...

  clearflags ();

  MyGetOptions(
    $srch_help,
    "db=s" => \$db,
    "query=s" => \$query,
    "sort=s" => \$sort,
    "days=i" => \$rldate,
    "mindate=s" => \$mndate,
    "maxdate=s" => \$mxdate,
    "datetype=s" => \$dttype,
    "label=s" => \$lbl,
    "pub=s" => \$pub,
    "feature=s" => \$feature,
    "location=s" => \$location,
    "molecule=s" => \$molecule,
    "organism=s" => \$organism,
    "source=s" => \$source,
    "status=s" => \$status,
    "type=s" => \$gtype,
    "clean" => \$clean,
    "field=s" => \$field,
    "word" => \$word,
    "drop" => \$drop,
    "trim" => \$trim,
    "trunc" => \$trunc,
    "spell" => \$spell,
    "split=s" => \$split,
    "merge=s" => \$meadow,
    "pairs=s" => \$pair,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "esearch $version\n";
    print $srch_help;
    return;
  }

  if ( -t STDIN and not @ARGV ) {
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    die "ERROR in search input: $err\n\n";
  }

  if ( $db ne "" ) {
    $dbase = $db;
  }

  $dbase = lc($dbase);

  if ( $dbase eq "" ) {
    die "Must supply -db database on command line\n";
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  binmode STDOUT, ":utf8";

  # support all efilter shortcut flags in esearch (undocumented)
  $query = process_extras ( $query, $pub, $feature, $location, $molecule, $organism, $source, $status, $gtype );

  if ( $query eq "" ) {
    die "Must supply -query search expression on command line\n";
  }

  $url = $base . $esearch;

  $query = map_labels ($query);
  $query = map_macros ($query);

  # multi-step query cleaning (undocumented)
  if ( $clean ) {
    if ( $query =~ /^(.*)\(.+\)(.*)$/ ) {
      $query = "$1 $2";
    }
    if ( $query =~ /^ +(.+)$/ ) {
      $query = $1;
    }
    if ( $query =~ /^(.+) +$/ ) {
      $query = $1;
    }
    $query =~ s/ +/ /g;
    $query = remove_stop_words ($query);
    $query = spell_check_query ($dbase, $query);
  }

  # remove punctuation from query (undocumented)
  if ( $word ) {
    $query = remove_punctuation ($query);
  }

  # drop stop words from query (undocumented)
  if ( $drop ) {
    $query = remove_stop_words ($query);
  }

  # trim words within parentheses (undocumented)
  if ( $trim ) {
    if ( $query =~ /^(.*)\(.+\)(.*)$/ ) {
      $query = "$1 $2";
    }
  }

  # truncate words at first parenthesis (undocumented)
  if ( $trunc ) {
    if ( $query =~ /^(.+)\(.*$/ ) {
      $query = $1;
    }
  }

  # remove leading, trailing, and multiple spaces
  if ( $query =~ /^ +(.+)$/ ) {
    $query = $1;
  }
  if ( $query =~ /^(.+) +$/ ) {
    $query = $1;
  }
  $query =~ s/ +/ /g;

  # spell check each query word
  if ( $spell ) {
    $query = spell_check_query ($dbase, $query);
  }

  # force each query word to be separately fielded, combined with AND (undocumented)
  if ( $split ne "" ) {
    $query = field_each_word ($split, $query);
  }

  # force each query word to be separately fielded, combined with OR (undocumented)
  if ( $meadow ne "" ) {
    $query = merge_each_word ($meadow, $query);
  }

  # -field combines -drop and -split (-field TITL produces same behavior as Web PubMed)
  if ( $field ne "" ) {
    $query = remove_stop_words ($query);
    $query = field_each_word ($field, $query);
  }

  # -pairs separately fields query word pairs, breaking chain at stop words
  if ( $pair ne "" ) {
    $query = remove_punctuation ($query);
    if ( $query =~ /^ +(.+)$/ ) {
      $query = $1;
    }
    if ( $query =~ /^(.+) +$/ ) {
      $query = $1;
    }
    $query =~ s/ +/ /g;
    $query = field_each_pair ($pair, $query);
  }

  $enc = uri_escape($query);
  $arg = "db=$dbase&term=$enc";
  if ( $web ne "" ) {
    $arg .= "&WebEnv=$web";
  }

  if ( $sort ne "" ) {
    if ( $sort eq "Relevance" ) {
      $sort = "relevance";
    }
    $arg .= "&sort=$sort";
  }
  $arg .= "&retmax=0&usehistory=y";

  if ( $rldate > 0 ) {
    $arg .= "&reldate=$rldate";
    if ( $dttype eq "" ) {
      $dttype = "PDAT";
    }
  }
  if ( $mndate ne "" and $mxdate ne "" ) {
    $arg .= "&mindate=$mndate&maxdate=$mxdate";
    if ( $dttype eq "" ) {
      $dttype = "PDAT";
    }
  }
  if ( $dttype ne "" ) {
    $arg .= "&datetype=$dttype";
  }

  $wb = $web;

  $output = do_post ($url, $arg, $tool, $email, true);

  $web = "";
  $key = "";
  $num = "";
  $err = "";
  my $trn = "";

  $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
  $key = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
  $num = $1 if ($output =~ /<Count>(\S+)<\/Count>/);
  $err = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);
  $trn = $1 if ($output =~ /<QueryTranslation>(.+?)<\/QueryTranslation>/i);

  if ( $err ne "" ) {
    write_edirect ( "", "", "", "", "", $err, "", "" );
    close (STDOUT);
    die "ERROR in search output: $err\nURL: $arg\n\n";
  }

  if ( $web eq "" ) {
    die "WebEnv value not found in search output - WebEnv1 $wb\n";
  }
  if ( $key eq "" ) {
    die "QueryKey value not found in search output - WebEnv1 $wb\n";
  }

  if ( $wb ne "" and $web ne $wb ) {
    $err = "WebEnv mismatch in search output - WebEnv1 $wb, WebEnv2 $web";
    write_edirect ( "", $wb, "", "", "", $err, "", "" );
    close (STDOUT);
    die "WebEnv value changed in search output - WebEnv1 $wb\nWebEnv2 $web\n";
  }

  if ( $num eq "" ) {
    die "Count value not found in search output - WebEnv1 $web\n";
  }

  if ( $lbl ne "" and $key ne "" ) {
    $labels{"$lbl"} = "$key";
  }

  write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );

  if ( $verbose ) {
    my $seconds = "";
    my $end_time = Time::HiRes::time();
    my $elapsed = $end_time - $begin_time;
    if ( $elapsed > 0.0005 ) {
      if ( $elapsed =~ /(\d+)\.(\d\d\d)\d+/ ) {
        $seconds = "$1.$2";
      }
    }
    if ( $seconds ne "" ) {
      print STDERR "Elapsed time is $seconds seconds\n";
    }
  }

  if ( $log ) {
    if ( $trn ne "" ) {
      print STDERR "$trn\n";
    }
  }
}

#  eaddr returns the current user's e-mail address

sub eaddr {
  my $addr = get_email ();
  print "$addr\n";
}

#  etest is an unadvertised function for development

sub etest {
#  $addr = get_email ();
#  print "e-mail:  $addr\n";
}

# main block dispatches control to appropriate subroutine

if ( scalar @ARGV > 0 and $ARGV[0] eq "-version" ) {
  print "$version\n";
} elsif ( $fnc eq "-search" ) {
  esrch ();
} elsif ( $fnc eq "-link" ) {
  elink ();
} elsif ( $fnc eq "-filter" ) {
  efilt ();
} elsif ( $fnc eq "-summary" ) {
  eftch ();
} elsif ( $fnc eq "-fetch" ) {
  eftch ();
} elsif ( $fnc eq "-info" ) {
  einfo ();
} elsif ( $fnc eq "-post" ) {
  epost ();
} elsif ( $fnc eq "-spell" ) {
  espel ();
} elsif ( $fnc eq "-citmatch" ) {
  ecitmtch ();
} elsif ( $fnc eq "-proxy" ) {
  eprxy ();
} elsif ( $fnc eq "-contact" ) {
  ecntc ();
} elsif ( $fnc eq "-notify" ) {
  entfy ();
} elsif ( $fnc eq "-address" ) {
  eaddr ();
} elsif ( $fnc eq "-test" ) {
  etest ();
} else {
  die "Function name '$fnc' not recognized\n";
}

# close input and output files

close (STDIN);
close (STDOUT);
close (STDERR);
