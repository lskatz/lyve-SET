# NOTE: Derived from blib/lib/Config/Simple.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Config::Simple;

#line 1503 "blib/lib/Config/Simple.pm (autosplit into blib/lib/auto/Config/Simple/write_string.al)"
# -------------------
# deprecated methods
# -------------------

sub write_string {
  my $self = shift;

  return $self->as_string(@_);
}

# end of Config::Simple::write_string
1;
