# NOTE: Derived from blib/lib/Config/Simple.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Config::Simple;

#line 1468 "blib/lib/Config/Simple.pm (autosplit into blib/lib/auto/Config/Simple/EXISTS.al)"
sub EXISTS {
  my ($self, $key) = @_;

  my $vars = $self->vars();
  return exists $vars->{$key};
}

# end of Config::Simple::EXISTS
1;
