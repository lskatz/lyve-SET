# NOTE: Derived from blib/lib/Config/Simple.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Config::Simple;

#line 1343 "blib/lib/Config/Simple.pm (autosplit into blib/lib/auto/Config/Simple/import_names.al)"
# imports names into the caller's namespace as global variables.
# I'm not sure how secure this method is. Hopefully someone will
# take a look at it for me
sub import_names {
  my ($self, $namespace) = @_;

  unless ( defined $namespace ) {    
    $namespace = (caller)[0];
  }
  if ( $namespace eq 'Config::Simple') {
    croak "You cannot import into 'Config::Simple' package";
  }
  my %vars = $self->vars();
  no strict 'refs';
  while ( my ($k, $v) = each %vars ) {
    $k =~ s/\W/_/g;
    ${$namespace . '::' . uc($k)} = $v;
  }
}

# end of Config::Simple::import_names
1;
