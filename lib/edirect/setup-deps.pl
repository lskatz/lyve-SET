#!/usr/bin/env perl
use warnings;
use strict;
use CPAN::MyConfig;
use CPAN;
use CPAN::HandleConfig;
use CPAN::Shell;

my $already_configured_cpan = 0;
my $root;

BEGIN {
    sub CheckAvailability
    {
        my $code = "require $_[0]";
        if (@_ > 1) {
            $code .= "; $_[1]";
        }
        eval $code;
        if ($@) {
            if ($already_configured_cpan) {
                print STDERR "Missing $_[0]; CPAN already initialized.\n";
            } else {
                print STDERR "Missing $_[0]; initializing CPAN.\n";
                CPAN::HandleConfig->load(autoconfig => 1, auto_pick => 1,
                                         doit => 1);
                CPAN::Shell::setup_output;
                CPAN::Index->reload;
                $already_configured_cpan = 1;
            }
            return 0;
        } else {
            print STDERR "Found $_[0].\n";
            return 1;
        }
    }

    alarm(3600);

    $root = $INC[0];
    if ($root !~ s:/_cpan$::) {
        die "Please invoke this script via .../setup.sh."
    }

    if ( ! exists $CPAN::Config{cpan_home} ) {
        $CPAN::Config = {
            'auto_commit'  => 1,
            'cpan_home'    => "$root/_cpan",
            'ftp_passive'  => 1,
            'install_help' => 'manual',
            'urllist'      => [@CPAN::Defaultsites]
        };
    }
    # $CPAN::DEBUG ||= $CPAN::DEBUG{'FTP'};

    if ( !CheckAvailability('local::lib',
                            'die unless $local::lib::VERSION >= 2') 
        &&  ! -d "$root/aux/lib/perl5/local" ) {
        my $ll = CPAN::Shell->expandany('local::lib');
        $ll->get;
        system('mkdir', '-p', "$root/aux/lib/perl5/local");
        system('cp', $ll->distribution->dir . "/lib/local/lib.pm",
               "$root/aux/lib/perl5/local/lib.pm");
    }
}

use lib "$root/aux/lib/perl5";
use local::lib("$root/aux", '--no-create');
my @lwp_deps = qw(Encode::Locale File::Listing
                  HTML::Parser HTML::Tagset HTML::Tree
                  HTTP:Cookies HTTP::Date HTTP::Message HTTP::Negotiate
                  IO::Socket::SSL LWP::MediaTypes LWP::Protocol::https
                  Net::HTTP URI WWW::RobotRules Mozilla::CA);
for my $module (@lwp_deps, 'Time::HiRes') {
    if ( ! CheckAvailability($module) ) {
        CPAN::Shell->install($module);
    }
}
if ( ! CheckAvailability('LWP') ) {
    CPAN::Shell->install('Bundle::LWP');
}
