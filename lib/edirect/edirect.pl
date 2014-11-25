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

# usage - edirect.pl -function arguments

use Data::Dumper;
use Encode;
use Getopt::Long;
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

$base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";

$efetch   = "efetch.fcgi";
$einfo    = "einfo.fcgi";
$elink    = "elink.fcgi";
$epost    = "epost.fcgi";
$esearch  = "esearch.fcgi";
$esummary = "esummary.fcgi";

# EDirect version number

$version = "1.00";

# utility subroutines

sub clearflags {
  @rest = ();
  %labels = ();
  %macros = ();
  $alias = "";
  $basx = "";
  $batch = false;
  $cmd = "";
  $complexity = 0;
  $db = "";
  $dbase = "";
  $dbs = "";
  $dbto = "";
  $debug = false;
  $dttype = "";
  $emaddr = "";
  $email = "";
  $err = "";
  $extrafeat = -1;
  $field = "";
  $holding = "";
  $http = "";
  $id = "";
  $just_num = false;
  $key = "";
  $lbl = "";
  $log = false;
  $max = 0;
  $min = 0;
  $mndate = "";
  $mode = "";
  $mxdate = "";
  $name = "";
  $neighbor = false;
  $num = "";
  $output = "";
  $pipe = false;
  $query = "";
  $related = false;
  $rldate = 0;
  $seq_start = 0;
  $seq_stop = 0;
  $silent = false;
  $stp = "";
  $strand = "";
  $tool = "";
  $tuul = "";
  $type = "";
  $verbose = false;
  $web = "";
}

# gets a live UID for any database

sub get_zero_uid {

  my $db = shift (@_);

  my $val = "";

  %zeroUidHash = (
    'assembly'         =>  '443538',
    'bioproject'       =>  '174115',
    'biosample'        =>  '1001166',
    'biosystems'       =>  '105909',
    'blastdbinfo'      =>  '1023214',
    'books'            =>  '1371014',
    'cdd'              =>  '189363',
    'clone'            =>  '18646800',
    'dbvar'            =>  '6173073',
    'epigenomics'      =>  '13103',
    'gap'              =>  '872875',
    'gds'              =>  '200022309',
    'gencoll'          =>  '398148',
    'gene'             =>  '3667',
    'genome'           =>  '52',
    'genomeprj'        =>  '72363',
    'geoprofiles'      =>  '16029743',
    'homologene'       =>  '510',
    'medgen'           =>  '162753',
    'mesh'             =>  '68007328',
    'ncbisearch'       =>  '15032',
    'nlmcatalog'       =>  '0404511',
    'nuccore'          =>  '1322283',
    'nucest'           =>  '338968657',
    'nucgss'           =>  '168803471',
    'nucleotide'       =>  '1322283',
    'omia'             =>  '1920',
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
    'toolkit'          =>  '242011',
    'unigene'          =>  '1132160',
    'unists'           =>  '40791'
  ); 

  if ( defined $zeroUidHash{$db} ) {
    $val = $zeroUidHash{$db};
  }

  return $val;
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
      if ( $1 == "130" ) {
        $base = "http://eutils.be-md.ncbi.nlm.nih.gov/entrez/eutils/";
      } elsif ( $1 == "165" ) {
        $base = "http://eutils.st-va.ncbi.nlm.nih.gov/entrez/eutils/";
      }
    }
    return;
  }

  if ( $basx =~ /\(#/ ) {
    $basx = map_macros ($basx);
  }

  if ( $basx !~ /^http:\/\// ) {
    $basx = "http://" . $basx;
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

  $keyx = "";
  $numx = "";
  $errx = "";

  $output = get ($url);

  if ( $output eq "" ) {
    print STDERR "No get_count output returned from '$url'\n";
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

  $data = get ($url);

  if ( $data eq "" ) {
    print STDERR "No get_uids output returned from '$url'\n";
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

    if ( $rslt eq "" ) {
      print STDERR "No do_get output returned from '$urlx'\n";
    }

    if ( $debug ) {
      print STDERR "$rslt\n";
    }

    return \$rslt;
  }

  $usragnt = new LWP::UserAgent (timeout => 300);

  $req = new HTTP::Request POST => "$urlx";
  $req->content_type('application/x-www-form-urlencoded');
  $req->content("$argx");

  $res = $usragnt->request ( $req );

  if ( $res->is_success) {
    $rslt = $res->content_ref;
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

# subroutines for each -function

# ecntc prepares the requested tool and email arguments for an EUtils pipe

sub ecntc {

  # ... | edirect.pl -contact -email darwin@beagle.edu -tool edirect_test | ...

  clearflags ();

  GetOptions (
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

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

sub efilt {

  # ... | edirect.pl -filter -query "bacteria [ORGN]" -days 365 | ...

  clearflags ();

  GetOptions (
    "query=s" => \$query,
    "days=i" => \$rldate,
    "mindate=s" => \$mndate,
    "maxdate=s" => \$mxdate,
    "datetype=s" => \$dttype,
    "label=s" => \$lbl,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $query eq "" && $rldate < 1 and $mndate eq "" and $mxdate eq "" ) {
    die "Must supply -query or -days or -mindate and -maxdate arguments on command line\n";
  }

  ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    die "ERROR in filt input: $err\n\n";
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  binmode STDOUT, ":utf8";

  if ( $dbase ne "" and $web ne "" and $key eq "" and $num eq "0" ) {
    write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
    # close (STDOUT);
    # die "QueryKey value not found in filter input\n";
    return;
  }

  test_edirect ( $dbase, $web, $key, $num, "filter" );

  $url = $base . $esearch;

  $arg = "db=$dbase&query_key=$key&WebEnv=$web";
  $arg .= "&retmax=0&usehistory=y";
  if ( $query ne "" ) {
    $query = map_labels ($query);
    $query = map_macros ($query);
    $enc = uri_escape($query);
    $arg .= "&term=$enc";
  }
  if ( $rldate > 0 ) {
    $arg .= "&reldate=$rldate";
    if ( $dttype ne "" ) {
      $dttype = "PDAT";
    }
  }
  if ( $mndate ne "" and $mxdate ne "" ) {
    $arg .= "&mindate=$mndate&maxdate=$mxdate";
    if ( $dttype ne "" ) {
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

  $output = "";

  for ( $tries = 0; $tries < 3 && $web eq ""; $tries++) {
    $output = do_post ($url, $arg, $tool, $email, true);

    $err = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);
    if ( $err eq "" ) {
      $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
      $key = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
      $num = $1 if ($output =~ /<Count>(\S+)<\/Count>/);
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
    die "ERROR in filt output: $err\nURL: $arg\nResult: $output\n\n";
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
}

# efetch -format docsum calls esmry to retrieve document summaries

sub esmry {

  my $dbase = shift (@_);
  my $web = shift (@_);
  my $key = shift (@_);
  my $num = shift (@_);
  my $id = shift (@_);
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

  if ( $dbase ne "" and $id ne "" ) {

    if ( $id eq "0" ) {

      # id "0" returns a live UID for any database

      $id = get_zero_uid ($dbase);
    }

    $url = $base . $esummary;

    $arg = "db=$dbase&id=$id";
    $arg .= "&version=2.0";

    $data = do_post ($url, $arg, $tool, $email, true);

    if ($data =~ /<eSummaryResult>/i and $data =~ /<ERROR>(.+?)<\/ERROR>/i) {
      $err = $1;
      if ( $err ne "" ) {
        if ( ! $silent ) {
          print STDERR "ERROR in esummary: $err\n";
        }
      }
    } else {
      $data =~ s/<DocumentSummary uid=\"(\d+)\">/<DocumentSummary><Id>$1<\/Id>/g;

      Encode::_utf8_on($data);

      print "$data";
    }

    return;
  }

  if ( $dbase ne "" and $web ne "" and $key eq "" and $num eq "0" ) {
    # die "QueryKey value not found in summary input\n";
    return;
  }

  test_edirect ( $dbase, $web, $key, $num, "summary" );

  $stpminusone = $stp - 1;

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

      $data =~ s/<DocumentSummary uid=\"(\d+)\">/<DocumentSummary><Id>$1<\/Id>/g;

      Encode::_utf8_on($data);

      print "$data";
    }

    sleep 1;
  }
}

# eftch can read all arguments from the command line or participate in an EUtils pipe

sub eftch {

  # ... | edirect.pl -fetch -format gp | ...

  clearflags ();

  GetOptions (
    "db=s" => \$db,
    "id=s" => \$id,
    "format=s" => \$type,
    "mode=s" => \$mode,
    "seq_start=i" => \$seq_start,
    "seq_stop=i" => \$seq_stop,
    "strand=s" => \$strand,
    "complexity=i" => \$complexity,
    "extrafeat=i" => \$extrafeat,
    "start=i" => \$min,
    "stop=i" => \$max,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "pipe" => \$pipe,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  # convert spaces between UIDs to commas

  $id =~ s/ /,/g;
  $id =~ s/,,/,/g;

  # "-format xml" is a shortcut for "-format full -mode xml"

  if ( $type eq "xml" and $mode eq "" ) {
    $type = "full";
    $mode = "xml";
  }

  if ( -t STDIN and not @ARGV ) {
  } elsif ( $db ne "" and $id ne "" ) {
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    die "ERROR in fetch input: $err\n\n";
  }

  if ( $dbase eq "" ) {
    $dbase = $db;
  }

  if ( $type eq "" and $dbase ne "" ) {
    if ( get_zero_uid ($dbase) eq "" ) {
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

    esmry ( $dbase, $web, $key, $num, $id, $min, $max, $tool, $email,
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

    # use larger chunk for UID format
    $chunk = 1000;
    for ( $start = $min; $start < $max; $start += $chunk ) {

      my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $max, $tool, $email );
      foreach $uid (@ids) {
        print "$uid\n";
      }
    }

    return;
  }

  if ( $dbase ne "" and ( $type eq "URL" or $type eq "url" ) ) {

    if ( $id ne "" ) {

      my @ids = split (',', $id);
      foreach $uid (@ids) {
        print "http://www.ncbi.nlm.nih.gov/$dbase/$uid\n";
      }

      return;
    }

    if ( $web eq "" ) {
      die "WebEnv value not found in fetch input\n";
    }

    if ( $pipe ) {
      write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
    }

    # use larger chunk for URL format
    $chunk = 1000;
    for ( $start = $min; $start < $max; $start += $chunk ) {

      my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $max, $tool, $email );
      foreach $uid (@ids) {
        print "http://www.ncbi.nlm.nih.gov/$dbase/$uid\n";
      }
    }

    return;
  }

  if ( $dbase ne "" and $id ne "" ) {

    if ( $id eq "0" ) {

      # id "0" returns a live UID for any database

      $id = get_zero_uid ($dbase);
    }

    $url = $base . $efetch;

    $arg = "db=$dbase&id=$id";
    $arg .= "&rettype=$type&retmode=$mode";
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

    print $$data;

    return;
  }

  if ( $dbase ne "" and $web ne "" and $key eq "" and $num eq "0" ) {
    # die "QueryKey value not found in fetch input\n";
    return;
  }

  test_edirect ( $dbase, $web, $key, $num, "fetch" );

  $stpminusone = $stp - 1;

  # use small chunk because fetched records could be quite large
  $chunk = 100;
  for ( $start = $min; $start < $max; $start += $chunk ) {
    $url = $base . $efetch;

    $chkx = $chunk;
    if ( $start + $chkx > $max ) {
      $chkx = $max - $start;
    }

    $arg = "db=$dbase&query_key=$key&WebEnv=$web";
    $arg .= "&rettype=$type&retmode=$mode";
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

      print $$data;
    }

    sleep 1;
  }
}

# einfo obtains names of databases or names of fields and links per database

sub einfo {

  # ... | edirect.pl -info -db pubmed | ...

  clearflags ();

  GetOptions (
    "db=s" => \$db,
    "dbs" => \$dbs,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

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

  if ($output =~ /IsTrunc/ or $output =~ /IsRang/) {
    # print STDERR "New server is out with IsTrunc/IsRang - update executable\n";
  }

  if ( $output eq "" ) {
    print STDERR "No einfo output returned from '$url'\n";
  }

  if ( $debug ) {
    print STDERR "$output\n";
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
    die "ERROR in link output: $err\nWebEnv: $wb\nURL: $arg\nResult: $output\n\n";
  }

  if ( $web eq "" ) {
    die "WebEnv value not found in link output - WebEnv1 $wb\n";
  }
  if ( $key eq "" ) {
    write_edirect ( $dbto, $wb, $key, "0", $stp, $err, $tool, $email );
    # close (STDOUT);
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
  }

  $dbase = $dbto;

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

sub elink {

  # ... | edirect.pl -link -target structure | ...

  clearflags ();

  GetOptions (
    "db=s" => \$db,
    "id=s" => \$id,
    "target=s" => \$dbto,
    "name=s" => \$name,
    "related" => \$related,
    "neighbor" => \$neighbor,
    "cmd=s" => \$cmd,
    "mode=s" => \$mode,
    "batch" => \$batch,
    "holding=s" => \$holding,
    "label=s" => \$lbl,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

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
    die "ERROR in link input: $err\n\n";
  }

  if ( $dbase eq "" ) {
    $dbase = $db;
  }

  if ( $dbto eq "" ) {
    if ( $cmd eq "acheck" or
         $cmd eq "ncheck" or
         $cmd eq "lcheck" or
         $cmd eq "llinks" or
         $cmd eq "llinkslib" or
         $cmd eq "prlinks" ) {
      $dbto = $dbase;
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

  if ( $dbase ne "" and $id ne "" ) {

    # process db and id command-line arguments instead of getting from history

    $url = $base . $elink;

    $arg = "dbfrom=$dbase&db=$dbto&cmd=$cmd&linkname=$name";
    if ( $mode ne "" ) {
      $arg .= "&retmode=$mode";
    }
    $arg .= "&id=$id";

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

    $chunk = 200;
    for ( $start = 0; $start < $num; $start += $chunk ) {

      my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $num, $tool, $email );

      $url = $base . $elink;

      $arg = "dbfrom=$dbase&db=$dbto&cmd=$cmd&linkname=$name";
      if ( $mode ne "" ) {
        $arg .= "&retmode=$mode";
      }
      $arg .= "&id=";
      $arg .= join ('&id=', @ids);

      $data = do_post ($url, $arg, $tool, $email, true);

      print "$data";
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

  $arg = "dbfrom=$dbase&db=$dbto&query_key=$key&WebEnv=$web";
  $arg .= "&cmd=$cmd&linkname=$name";
  if ( $mode ne "" ) {
    $arg .= "&retmode=$mode";
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

  $mbase = "http://www.ncbi.nlm.nih.gov/Structure/mmdb/";
  $mprog = "mmdbsrv.cgi";

  if ( $id ne "" ) {

    my @ids = split (',', $id);
    foreach $uid (@ids) {
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
    }
  }
}

# entfy sends e-mail

sub entfy {

  # ... | edirect.pl -notify

  clearflags ();

  GetOptions (
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

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
        $txt = "echo \"http://www.ncbi.nlm.nih.gov/$dbase/$uid\n\"";
        $str = "mail -s \"A new $dbase record is in Entrez\" $email";
        system "$txt | $str";
      }
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
      if ( $dttype ne "" ) {
        $dttype = "PDAT";
      }
    }
    if ( $mndate ne "" and $mxdate ne "" ) {
      $arg .= "&mindate=$mndate&maxdate=$mxdate";
      if ( $dttype ne "" ) {
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

sub epost {

  # ... | edirect.pl -post -db nucleotide -format uid | ...

  clearflags ();

  GetOptions (
    "db=s" => \$db,
    "id=s" => \$id,
    "format=s" => \$field,
    "label=s" => \$lbl,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

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
    die "ERROR in post input: $err\n\n";
  }

  if ( $dbase eq "" ) {
    $dbase = $db;
  }

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

  my $combo = "";
  my $pfx = "";
  my $loops = 0;

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

        ( $web, $key ) = post_chunk ( $dbase, $web, $key, $tool, $email, $ids, "" );

        $combo .= $pfx . "#" . $key;
        $pfx = " OR ";
        $loops++;
      }
    }

  } else {

    while ( @rest ) {
      my @chunk = splice(@rest, 0, 2000);

      $query = join (' OR ', @chunk);

      if ( $query eq "" ) {
        die "Must pipe data into stdin\n";
      }

      if ( $field eq "ACCN" or $field eq "accn" or $field eq "ACC" or $field eq "acc" ) {
        if ( $dbase eq "nucleotide" or $dbase eq "nuccore" or $dbase eq "est" or
             $dbase eq "gss" or $dbase eq "protein" ) {
          $query =~ s/\./_/g;
          $field = "ACCN";
        }
      } elsif ( ! $just_num ) {
        die "Non-numeric value found in post input\n";
      }

      $query .= " [$field]";

      ( $web, $key ) = post_chunk ( $dbase, $web, $key, $tool, $email, "", $query );

      $combo .= $pfx . "#" . $key;
      $pfx = " OR ";
      $loops++;
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

# eprxy reads a file of query proxies, can also pipe from stdin

sub eprxy {

  # ... | edirect.pl -proxy -file file_name | ...

  clearflags ();

  GetOptions (
    "pipe" => \$pipe,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

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

sub esrch {

  # ... | edirect.pl -search -db nucleotide -query "M6506* [ACCN] OR 1322283 [UID]" -days 365 | ...

  clearflags ();

  GetOptions (
    "db=s" => \$db,
    "query=s" => \$query,
    "days=i" => \$rldate,
    "mindate=s" => \$mndate,
    "maxdate=s" => \$mxdate,
    "datetype=s" => \$dttype,
    "label=s" => \$lbl,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "log" => \$log,
    "http=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

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
    die "Must supply -query search expression on command line\n";
  }

  $url = $base . $esearch;

  $query = map_labels ($query);
  $query = map_macros ($query);
  $enc = uri_escape($query);
  $arg = "db=$dbase&term=$enc";
  if ( $web ne "" ) {
    $arg .= "&WebEnv=$web";
  }
  $arg .= "&retmax=0&usehistory=y";
  if ( $rldate > 0 ) {
    $arg .= "&reldate=$rldate";
    if ( $dttype ne "" ) {
      $dttype = "PDAT";
    }
  }
  if ( $mndate ne "" and $mxdate ne "" ) {
    $arg .= "&mindate=$mndate&maxdate=$mxdate";
    if ( $dttype ne "" ) {
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

  $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
  $key = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
  $num = $1 if ($output =~ /<Count>(\S+)<\/Count>/);
  $err = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);

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
} elsif ( $fnc eq "-proxy" ) {
  eprxy ();
} elsif ( $fnc eq "-contact" ) {
  ecntc ();
} elsif ( $fnc eq "-notify" ) {
  entfy ();
} elsif ( $fnc eq "-test" ) {
  etest ();
} else {
  die "Function name '$fnc' not recognized\n";
}

# close input and output files

close (STDIN);
close (STDOUT);
close (STDERR);
