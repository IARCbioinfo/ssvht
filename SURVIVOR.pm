package SURVIVOR;

=head1 NAME
SURVIVOR

=head1 DESCRIPTION

This object perform specific operations on the output of SURVIVOR caller

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

#add missing values to a set of SURVIVOR SVs
sub normalize_sv{
   my $self=shift;
   my $svobj=shift;
   foreach my $item (@{$svobj->{entries}}){
     $item->{info}->{SVMETHOD}="SURVIVOR";
     my $len=$item->{info}->{SVLEN};
     #the lenght of the SV goes to positive
     $item->{info}->{SVLEN}=abs($len);
     my $type=$item->{info}->{SVTYPE};
     #SURVIVOR pon
     if($type eq "TRA" or $type eq "BND"){
       $item->{info}->{POS2}=$item->{info}->{END};
       #print Dumper($item);
     }
   }
}

1;
