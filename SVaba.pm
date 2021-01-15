package SVaba;

=head1 NAME
SVaba

=head1 DESCRIPTION

This object perform specific operations on the output of SVaba caller

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
   foreach my $item(@{$svobj->{entries}}){

     $item->{info}->{SVMETHOD}="SVaba";
##FORMAT=<ID=LR,Number=1,Type=Float,Description="Log-odds that this variant is REF vs AF=0.5">
##FORMAT=<ID=SL,Number=1,Type=Float,Description="Alignment-quality Scaled log-odds, where LO is LO * (MAPQ - 2*NM)/60">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of discordant-supported reads for this variant">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of spanning reads for this variants">
##FORMAT=<ID=PL,Number=.,Type=Float,Description="Normalized likelihood of the current genotype">
##FORMAT=<ID=GQ,Number=1,Type=String,Description="Genotype quality (currently not supported. Always 0)">
##FORMAT=<ID=LO,Number=1,Type=Float,Description="Log-odds that this variant is real vs artifact">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth of coverage: Number of reads covering site.">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allele depth: Number of reads supporting the variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Most likely genotype">
my $sv_sup=0;
my $ref_sup=0;
#we ask for PE support
if(!defined $item->{geno}[0]->{DR}){
  $item->{info}->{PE}=0;
}else{
  #we add all the PE support for REF and ALT
  #$item->{info}->{PE}=$item->{geno}[0]->{PR}[0]+$item->{geno}[0]->{PR}[1];
  $item->{info}->{PE}=$item->{geno}[0]->{DR}[0];
  $sv_sup+=$item->{geno}[0]->{DR}[0];
  #$ref_sup+=$item->{geno}[0]->{PR}[0];
}
if(!defined $item->{geno}->[0]->{SR}){
  $item->{info}->{SR}=0;
}else{
  #$item->{info}->{SR}=$item->{geno}[0]->{SR}[0]+$item->{geno}[0]->{SR}[1];
 $item->{info}->{SR}=$item->{geno}[0]->{SR}[0];
  $sv_sup+=$item->{geno}[0]->{SR}[0];
  #$ref_sup+=$item->{geno}[0]->{SR}[0];
}
#total support of the SV for ALT
$item->{info}->{PE_SR}=$item->{geno}[0]->{AD}[0];
#we compute the lengh of the SV

#we add total supporting reads for ref/alt
$item->{info}->{RFS}=$item->{geno}[0]->{DP}[0];
if($item->{info}->{RFS} == 0){
  $item->{info}->{RFS}=1;
}

#we add AF
$item->{info}->{RAF}=($item->{info}->{PE_SR})/($item->{info}->{RFS});
#We add some tags related to SVlen
my $type=$item->{info}->{SVTYPE};
if($type eq "DEL" or $type eq "DUP" or $type eq "INV" or $type eq "INS"){
       #my $l=abs($item->{info}->{SPAN});
       $item->{info}->{SVLEN}=$item->{info}->{SPAN};
       my ($chr_pos)=_parse_trans_name($item->{ALT});
       $item->{info}->{END}=$chr_pos->{POS2};
       $item->{info}->{SVLEN}=abs($item->{POS}-$item->{info}->{END});
}

 if($type eq "BND"){
   $item->{info}->{SVLEN}=0;
   $item->{info}->{SVTYPE}="TRA";#change BND for translocation
   my ($chr_pos)=_parse_trans_name($item->{ALT});
   $item->{info}->{CHR2}=$chr_pos->{CHR2};
   $item->{info}->{POS2}=$chr_pos->{POS2};
 }
#we add CIPOS and CIEND
$item->{info}->{CIPOS}="-100,100";
$item->{info}->{CIEND}="-100,100";
#print Dumper($item);

   }
}


#function that parse the orientation of translocations and inversions
sub _parse_trans_name{
  my ($value)=(@_);
  my $val=();
  #notation from VCFv4.3
  # s t[p[ piece extending to the right of p is joined after t
  # s t]p] reverse comp piece extending left of p is joined after t
  # s ]p]t piece extending to the left of p is joined before t
  # s [p[t reverse comp piece extending right of p is joined before t
  if($value =~/(\w+):(\d+)/){
    #print $1." ".$2."\n";
    $val->{CHR2}=$1;
    $val->{POS2}=$2;
  }
  return $val;
}

1;
