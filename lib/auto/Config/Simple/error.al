# NOTE: Derived from blib/lib/Config/Simple.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Config::Simple;

#line 1390 "blib/lib/Config/Simple.pm (autosplit into blib/lib/auto/Config/Simple/error.al)"
sub error {
  my ($self, $msg) = @_;

  if ( $msg ) {
    $errstr = $msg;
  }
  return $errstr;
}

# end of Config::Simple::error
1;
