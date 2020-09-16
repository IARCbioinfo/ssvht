package SVE;

=head1 NAME
SV

=head1 DESCRIPTION

This object perform several operation on SV lines

=head2 Available methods


=cut

use strict;
use Data::Dumper;

sub new{
  my ($packagename, $line, $lgeno) = @_;
  chomp $line;
  my $self = {el => $line,lg=>$lgeno};
  bless ($self, $packagename);
  $self->_create_entry();
  return ($self);
}


sub _create_entry{
      my $self=shift;
      #my $lgeno=shift;

      #parse subparts of the entry
      my $tags=$self->_get_tags();

      my $geno;
      if($self->{lg}){
        $geno=$self->_get_geno();
      }
      $self->{info}=$tags;
      $self->{geno}=$geno;
      ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  B00JAM1 B00JAM2

      my @data=split("\t",$self->{el});
      $self->{CHROM}=$data[0];
      $self->{POS}=$data[1];
      $self->{ID}=$data[2];
      $self->{REF}=$data[3];
      $self->{ALT}=$data[4];
      $self->{QUAL}=$data[5];
      $self->{FILTER}=$data[6];
      #$self->{INFO}=$data[0];
      #my $sv=();

      #$sv->{brp1_chr}=$data[0];
      #$sv->{brp1_pos}=$data[1];

      #$sv->{brp2_chr}=$tags->{CHR2};
      #$sv->{brp2_pos}=$tags->{END};
      #print Dumper(@data,$tags,$sv, $geno);
}

sub check_mandatory_tags{




  my @man=["SVTYPE","SVLEN","PR","SR","CHR2","END"];
  foreach my $v (@man){
    }


}

#parse the genotype information of the tool
sub _get_geno{
  my $self=shift;
  my $geno=();
  #chomp $entry;
  my @d=split("\t",$self->{el});
  my @gtags=split(":",$d[8]);
  for(my $i=9; $i<=$#d; $i++){
    my @d2=split(":",$d[$i]);
    my $tags=();
    for(my $j=0; $j<=$#gtags; $j++){
          my @gv=split(",",$d2[$j]);
          push(@{$tags->{$gtags[$j]}},$_) foreach (@gv);

    }
    push(@{$geno},$tags);
  }

  #we return the set of genotypes
  return $geno;
}

#parse the tags of the tool
sub _get_tags{
    my $self=shift;
    my $tags=();
    #chomp $entry;
    my @d=split("\t",$self->{el});
    my @vals=split(";",$d[7]);
    foreach my $t(@vals){
      my ($tag,$v)=split("=",$t);
          if(!defined $v){
            $v=1;
          }
          $tags->{$tag}=$v;
    }

    return $tags;
}

1;
