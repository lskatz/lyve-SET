# NOTE: Derived from blib/lib/Config/Simple.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Config::Simple;

#line 1490 "blib/lib/Config/Simple.pm (autosplit into blib/lib/auto/Config/Simple/NEXTKEY.al)"
sub NEXTKEY {
  my $self = shift;

  unless ( defined $self->{_TIED_HASH} ) {
    $self->{_TIED_HASH} = $self->vars();
  }
  return scalar each %{ $self->{_TIED_HASH} };
}

# end of Config::Simple::NEXTKEY
1;
