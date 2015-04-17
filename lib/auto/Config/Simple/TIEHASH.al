# NOTE: Derived from blib/lib/Config/Simple.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Config::Simple;

#line 1426 "blib/lib/Config/Simple.pm (autosplit into blib/lib/auto/Config/Simple/TIEHASH.al)"
#------------------
# tie() interface
#------------------

sub TIEHASH {
  my ($class, $file, $args) = @_;

  unless ( defined $file ) {
    croak "Usage: tie \%config, 'Config::Simple', \$filename";
  }  
  return $class->new($file);
}

# end of Config::Simple::TIEHASH
1;
