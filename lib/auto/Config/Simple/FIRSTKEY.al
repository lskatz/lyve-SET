# NOTE: Derived from blib/lib/Config/Simple.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Config::Simple;

#line 1477 "blib/lib/Config/Simple.pm (autosplit into blib/lib/auto/Config/Simple/FIRSTKEY.al)"
sub FIRSTKEY {
  my $self = shift;

  # we make sure that tied hash is created ONLY if the program
  # needs to use this functionality.
  unless ( defined $self->{_TIED_HASH} ) {    
    $self->{_TIED_HASH} = $self->vars();
  }
  my $temp = keys %{ $self->{_TIED_HASH} };
  return scalar each %{ $self->{_TIED_HASH} };
}

# end of Config::Simple::FIRSTKEY
1;
