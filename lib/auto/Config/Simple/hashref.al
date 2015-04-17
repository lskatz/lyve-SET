# NOTE: Derived from blib/lib/Config/Simple.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Config::Simple;

#line 1513 "blib/lib/Config/Simple.pm (autosplit into blib/lib/auto/Config/Simple/hashref.al)"
sub hashref {
  my $self = shift;

  return scalar( $self->vars() );
}

# end of Config::Simple::hashref
1;
