# NOTE: Derived from blib/lib/Config/Simple.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Config::Simple;

#line 1365 "blib/lib/Config/Simple.pm (autosplit into blib/lib/auto/Config/Simple/import_from.al)"
# imports names from a file. Compare with import_names.
sub import_from {
  my ($class, $file, $arg) = @_;

  if ( ref($class) ) {
    croak "import_from() is not an object method.";
  }
  # this is a hash support
  if ( defined($arg) && (ref($arg) eq 'HASH') ) {
    my $cfg = $class->new($file) or return;
    map { $arg->{$_} = $cfg->param($_) } $cfg->param();
    return $cfg;
  }
  # following is the original version of our import_from():
  unless ( defined $arg ) {
    $arg = (caller)[0];
  }  
  my $cfg = $class->new($file) or return;
  $cfg->import_names($arg);
  return $cfg;
}

# end of Config::Simple::import_from
1;
