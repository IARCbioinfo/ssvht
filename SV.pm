package SV;


=head1 NAME

SV

=head1 DESCRIPTION

This object perform several operation on SV  using VCF files

=head2 Available methods


=cut

use strict;
use Data::Dumper;
use SVE;
use DELLY;
use Manta;
use SVaba;

sub new{
  my ($packagename, $vcf) = @_;
  my $self = {vcffile => $vcf};
  bless ($self, $packagename);
  return ($self);
}


sub _get_autosomes_X_numbers{
  my @chrs=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X");
  return @chrs;
}

sub _get_autosomes_X_chr{
  my @chrs=("chr1","chr2","chr3","chr4","chr5","chr6","chr7",
  "chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
  "chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX");
  return @chrs;
}

sub _get_hash_autosomes_X_chr{
    my $self= shift;
    my $hchr=();
    foreach my $c($self->_get_autosomes_X_chr()){
        $hchr->{$c}=1;
    }
    return $hchr;
}

sub _parse_ctg_entry{
      my $self=shift;
      my $entry=shift;
      chomp $entry;
      my $s=index($entry,"<");
      my $e=index($entry,">");
      #print join(" ",$s,$e,$entry,substr($entry,$s+1,abs($s+1-$e)))."\n";
      my $ht=();
      foreach my $en(split(",",substr($entry,$s+1,abs($s+1-$e)))){
        my ($t,$v)=split("=",$en);
        #print join(" ",$t,$v);
          $ht->{$t}=$v;
      }
      #print join(" ",$s,$e,$entry,substr($entry,$s+3,abs($s+3-$e)))."\n";
      #print $e;
      #print Dumper($ht);
      return $ht;
}



sub norm_svs{
      my $self=shift;
      my $hchr=$self->_get_hash_autosomes_X_chr();
      open(VCF, $self->{vcffile}) or die "cannot open VCF file\n";
      while (my $line=<VCF>) {
        if ($line =~m/^#/){
             next;
        }else{
          #we have an entry
          my $entry= new SVE($line);
          #next if($entry->{CHR1} eq)
          #print Dumper($entry);
          push (@{$self->{entries}},$entry);
        }
      }

      #function that guest the caller
      my ($caller)=$self->_guess_caller();
      #add missing values to the entries and merge complex ones
      $caller->normalize_sv($self);
}


#function that guest the caller currently expect SVaba, Manta and Delly2
sub _guess_caller{
  my $self=shift;
  my $caller="unknow";

  foreach my $item(@{$self->{entries}}){
     if($item->{info}->{SVMETHOD} =~m/DELLY/){
       #print "Delly2\n";
       $caller="DELLY2";
     }elsif($item->{ID} =~m/Manta/){
       #print "Manta\n";
       $caller="Manta";
     }elsif(defined $item->{info}->{SCTG}){
       #print "SVaba\n";
       $caller="SVaba";
     }else{
       print "Unknow caller\n";
     }
  }

  my $tool=();
  if($caller eq "DELLY2"){
    $tool = new DELLY();
    #$tool->normalize_sv($self);
  }elsif($caller eq "Manta"){
    $tool = new Manta();
    #$tool->normalize_sv($self);
  }elsif($caller eq "SVaba"){
    $tool = new SVaba();
  }

  return $tool;
}



sub filter_by_region{
    my $self=shift;
    my $prefix=shift;
    my $hchr=$self->_get_hash_autosomes_X_chr();
    #print Dumper($hchr);
    #filter in header
    ##contig=<ID=chr1,length=248956422>
    open(VCF, $self->{vcffile}) or die "cannot open VCF file\n";
    my $p =substr($self->{vcffile},0,length($self->{vcffile})-3)."autosomes_plus_X.vcf";
    open(OUT, ">".$p) or die "cannot create the output file\n";

    while (my $line=<VCF>) {
      if ($line =~m/^#/){
        #print $line;
          if($line=~m/#contig=/){
              my ($e)=$self->_parse_ctg_entry($line);
              print OUT $line if(defined $hchr->{$e->{ID}});
          }else{
            print OUT $line;
          }
      }else{
        #is an entry
        my @tmp=split("\t",$line);
        print OUT $line if(defined $hchr->{$tmp[0]});
      }
    }
    close(OUT);
}
#for the moment insertion and deletions only
sub add_type_entry{
  my $self=shift;
  my $entry=shift;
  my $m=shift;

  chomp $entry;
  my @d=split("\t",$entry);

 #chrX    138269888       .       C       CAATGTCATCAGTTAAGGCAGGAACAGGCCATTTTCACTTCTTTTGTGGTGG    60      .       .       GT      1/1
 #PRECISE;SVMETHOD=PAFTS;SVTYPE=INS;SVLEN=325
 my $info="";
 #deletion
 if(length($d[3]) > length($d[4])){
     $info="PRECISE;SVMETHOD=$m;SVTYPE=DEL;SVLEN=-".length($d[3]);
 }else{
   #insertion
   $info="PRECISE;SVMETHOD=$m;SVTYPE=INS;SVLEN=".length($d[4]);
 }
return join("\t",@d[0 .. 6],$info,@d[8 .. $#d]);

}


sub add_sv_type{
  my $self=shift;
  my $m=shift;
  if(length($m) == 0 ){
    $m="Denovo";
  }

  open(VCF, $self->{vcffile}) or die "cannot open VCF file\n";
  my $p =substr($self->{vcffile},0,length($self->{vcffile})-3)."svtypes.vcf";
  open(OUT, ">".$p) or die "cannot create the output file\n";
  my $types=<<EOF;
##fileformat=VCFv4.3
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=INVDUP,Description="InvertedDUP with unknown boundaries">
##ALT=<ID=TRA,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
##INFO=<ID=UNRESOLVED,Number=0,Type=Flag,Description="An insertion that is longer than the read and thus we cannot predict the full size.">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
EOF
my $info=1;
  while (my $line=<VCF>) {
    if ($line =~m/^#/){
        if($line=~m/^##INFO/ and $info == 1){
           print OUT $types;
           print OUT $line;
           $info++;
        }else{
          print OUT $line;
        }
    }else{
        #we are in an entry
        print OUT $self->add_type_entry($line,$m)."\n";
    }
  }
    close(OUT);

}


sub _get_tags{
    my $self=shift;
    my $entry=shift;
    my $tags=();
    chomp $entry;
    my @d=split("\t",$entry);
    my @vals=split(";",$d[7]);
    foreach my $t(@vals){
      my ($tag,$v)=split("=",$t);
          $tags->{$tag}=$v;
    }
    return $tags;
}


sub keep_indels_svs{
      my $self=shift;
      open(VCF, $self->{vcffile}) or die "cannot open VCF file\n";
      my $p =substr($self->{vcffile},0,length($self->{vcffile})-3)."sv_insdel.vcf";
      open(OUT, ">".$p) or die "cannot create the output file\n";

      while (my $line=<VCF>) {
        if ($line =~m/^#/){
              print OUT $line;
        }else{
            #we are in an entry
            my $tags=$self->_get_tags($line);
            # print Dumper($tags);
            if($tags->{SVTYPE} eq "INS" or $tags->{SVTYPE} eq "DEL"){
                if(abs($tags->{SVLEN}) >= 50 ){
                  print OUT $line;
                }
            }
        }
      }
        close(OUT);
}


1; #EOM