# NOTE: Derived from blib/lib/Config/Simple.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Config::Simple;

#line 1418 "blib/lib/Config/Simple.pm (autosplit into blib/lib/auto/Config/Simple/verbose.al)"
sub verbose {
  DEBUG or return;
  carp "****[$0]: " .  join ("", @_);
}

# end of Config::Simple::verbose
1;
