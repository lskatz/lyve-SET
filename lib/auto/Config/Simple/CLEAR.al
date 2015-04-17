# NOTE: Derived from blib/lib/Config/Simple.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Config::Simple;

#line 1462 "blib/lib/Config/Simple.pm (autosplit into blib/lib/auto/Config/Simple/CLEAR.al)"
sub CLEAR {
  my $self = shift;
  map { $self->delete($_) } $self->param();
}

# end of Config::Simple::CLEAR
1;
