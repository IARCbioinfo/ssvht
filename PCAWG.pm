package PCAWG;

=head1 NAME
PCAWG

=head1 DESCRIPTION

This object read the PCAWG calls from transformed bedpe files

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

#add missing values to a set of PCAWG SVs
sub normalize_sv{
   my $self=shift;
   my $svobj=shift;
   #PCAWG items
   foreach my $item (@{$svobj->{entries}}){
     $item->{info}->{SVMETHOD}="PCAWG";
     my $len=$item->{info}->{SVLEN};
     #the lenght of the SV goes to positive
     $item->{info}->{SVLEN}=abs($len);
     #my $type=$item->{info}->{SVTYPE};
     #PCAWG CIPOS and CIEND intervals are equal 0
      #$item->{info}->{CIPOS}="100,100";
      #$item->{info}->{CIEND}="100,100";
   }

}

1;
