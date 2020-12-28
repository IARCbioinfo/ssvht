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
     #print Dumper($item);
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
     #insertions has this varible already
     #<ID=SVLEN,Number=1,Type=Integer,Description="Insertion length for SVTYPE=INS.">
     my $type=$item->{info}->{SVTYPE};
     if($type eq "DEL" or $type eq "DUP" or $type eq "INV"){
            my $l=abs($item->{POS}-$item->{info}->{END});
            $item->{info}->{SVLEN}=$l;
     }
     #we get genovars
    # FORMAT	GT	Genotype
#FORMAT	GL	Log10-scaled genotype likelihoods for RR
#FORMAT	GQ	Genotype Quality
#FORMAT	FT	Per-sample genotype filter
#FORMAT	RC	Raw high-quality read counts or base counts for the SV
#FORMAT	RCL	Raw high-quality read counts or base counts for the left control regio
#FORMAT	RCR	Raw high-quality read counts or base counts for the right control regi
#FORMAT	CN	Read-depth based copy-number estimate for autosomal sites
#FORMAT	DR	# high-quality reference pairs
#FORMAT	DV	# high-quality variant pairs
#FORMAT	RR	# high-quality reference junction reads
#FORMAT	RV	# high-quality variant junction reads
     #print Dumper($item->{geno});
     #print join(" ",$item->{geno}->[0]->{DV}[0],$item->{geno}->[0]->{RV}[0])."\n";
     #for tumor-only samples there is only one genome
     my $sv_sup=$item->{geno}->[0]->{DV}[0]+$item->{geno}->[0]->{RV}[0];
     my $ref_sup=$item->{geno}->[0]->{DR}[0]+$item->{geno}->[0]->{RR}[0];
      #we add total supporting reads for ref/alt
      $item->{info}->{RFS}=$sv_sup+$ref_sup;
      if($item->{info}->{RFS} == 0){
        $item->{info}->{RFS}=1;
      }
      #we add AF
      $item->{info}->{RAF}=($sv_sup)/($item->{info}->{RFS});


     #we change the BND type of DELLY to TRA = TRANSLOCATION
     ###ALT=<ID=BND,Description="Translocation">
     if($type eq "BND"){
             $item->{info}->{SVLEN}=0;
             $item->{info}->{SVTYPE}="TRA";#change BND for translocation
     }

     #print join("\t",$item->{ID},$item->{info}->{PE},$item->{info}->{SR},$type,$item->{info}->{SVLEN})."\n";
     #print Dumper($item);
   }
}

1;
