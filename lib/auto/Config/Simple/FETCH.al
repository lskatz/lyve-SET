# NOTE: Derived from blib/lib/Config/Simple.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Config::Simple;

#line 1440 "blib/lib/Config/Simple.pm (autosplit into blib/lib/auto/Config/Simple/FETCH.al)"
sub FETCH {
  my $self = shift;

  return $self->param(@_);
}

# end of Config::Simple::FETCH
1;
