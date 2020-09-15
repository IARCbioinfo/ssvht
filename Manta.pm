package Manta;

=head1 NAME
Manta

=head1 DESCRIPTION

This object perform specific operations on the output of Manta caller

=head2 Available methods


=cut

use strict;
use Data::Dumper;

sub new{
  my ($packagename) = @_;
  my $self = {};
  bless ($self, $packagename);
  return ($self);
}


#add missing values to a set of DELLY2 SVs
sub normalize_sv{
   my $self=shift;
   my $svobj=shift;
   foreach my $item($svobj->{entries}){
     print Dumper($item);
   }
}

1;
