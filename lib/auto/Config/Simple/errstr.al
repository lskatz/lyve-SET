# NOTE: Derived from blib/lib/Config/Simple.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Config::Simple;

#line 1525 "blib/lib/Config/Simple.pm (autosplit into blib/lib/auto/Config/Simple/errstr.al)"
sub errstr {
  my $self = shift;
  return $self->error(@_);
}

# end of Config::Simple::errstr
1;
