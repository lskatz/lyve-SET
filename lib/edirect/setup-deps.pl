#!/usr/bin/perl -w
use strict;
use CPAN::MyConfig;
use CPAN;
use CPAN::HandleConfig;
use CPAN::Shell;

my $root;

BEGIN {
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
    CPAN::HandleConfig->load(autoconfig => 1, auto_pick => 1, doit => 1);
    CPAN::Shell::setup_output;
    CPAN::Index->reload;

    my $ll = CPAN::Shell->expandany('local::lib');
    if ( ! $ll->inst_file  &&  ! -d "$root/aux/lib/perl5/local" ) {
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
                  LWP::MediaTypes LWP::Protocol::HTTPS
                  Net::HTTP URI WWW::RobotRules);
for my $module (@lwp_deps, 'Time::HiRes') {
    if ( ! CPAN::Shell->expandany($module)->inst_file ) {
        CPAN::Shell->install($module);
    }
}
if ( ! CPAN::Shell->expandany('LWP')->inst_file ) {
    CPAN::Shell->install('Bundle::LWP');
}
