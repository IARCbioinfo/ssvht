package GNOMAD;

=head1 NAME
GNOMAD

=head1 DESCRIPTION

This object read the GNOMAD calls from NCBI

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

#add missing values to a set of GNOMAD SVs
sub normalize_sv{
   my $self=shift;
   my $svobj=shift;
   #GNOMAD don't have predictions for
   foreach my $item (@{$svobj->{entries}}){
     $item->{info}->{SVMETHOD}="GNOMAD";
     my $len=$item->{info}->{SVLEN};
     #the lenght of the SV goes to positive
     $item->{info}->{SVLEN}=abs($len);
     my $type=$item->{info}->{SVTYPE};
     #GNOMAD CIPOS and CIEND intervals
      $item->{info}->{CIPOS}="100,100";
      $item->{info}->{CIEND}="100,100";
      my $chr=$item->{CHROM};
      #todo: check this on the fly
      $chr="chr".$chr;
      $item->{CHROM}=$chr;
   }

}

1;
