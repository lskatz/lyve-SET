# NOTE: Derived from blib/lib/Config/Simple.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Config::Simple;

#line 1403 "blib/lib/Config/Simple.pm (autosplit into blib/lib/auto/Config/Simple/dump.al)"
sub dump {
  my ($self, $file, $indent) = @_;

  require Data::Dumper;
  my $d = new Data::Dumper([$self], [ref $self]);
  $d->Indent($indent||2);
  if ( defined $file ) {
    sysopen(FH, $file, O_WRONLY|O_CREAT|O_TRUNC, 0666) or die $!;
    print FH $d->Dump();
    CORE::close(FH) or die $!;
  }
  return $d->Dump();
}

# end of Config::Simple::dump
1;
