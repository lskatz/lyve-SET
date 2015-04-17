# NOTE: Derived from blib/lib/Config/Simple.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Config::Simple;

#line 1447 "blib/lib/Config/Simple.pm (autosplit into blib/lib/auto/Config/Simple/STORE.al)"
sub STORE {
  my $self = shift;

  return $self->param(@_);
}

# end of Config::Simple::STORE
1;
