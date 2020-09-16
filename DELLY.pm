package DELLY;

=head1 NAME
DELLY

=head1 DESCRIPTION

This object perform specific operations on the output of delly2 caller

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
   foreach my $item (@{$svobj->{entries}}){

     $item->{info}->{SVMETHOD}="DELLY";
     #we ask for PE support
     if(!defined $item->{info}->{PE}){
       $item->{info}->{PE}=0;
     }
     if(!defined $item->{info}->{SR}){
       $item->{info}->{SR}=0;
     }
     #total support of the SV
     $item->{info}->{PE_SR}=$item->{info}->{PE}+$item->{info}->{SR};
     #we compute the lengh of the SV
     my $type=$item->{info}->{SVTYPE};
     if($type eq "DEL" or $type eq "DUP" or $type eq "INV"){
            my $l=abs($item->{POS}-$item->{info}->{END});
            $item->{info}->{SVLEN}=$l;
     }

     if($type eq "BND"){
             $item->{info}->{SVLEN}=0;
     }

     #print join("\t",$item->{ID},$item->{info}->{PE},$item->{info}->{SR},$type,$item->{info}->{SVLEN})."\n";
     #print Dumper($item);
   }
}

1;
