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
use SURVIVOR;
use GNOMAD;
use PCAWG;

sub new{
  my ($packagename, $vcf) = @_;
  my $self = {vcffile => $vcf};
  bless ($self, $packagename);
  return ($self);
}



sub _get_autosomes_X_numbers{
  my @chrs=("1","2","3","4","5","6","7",
            "8","9","10","11","12","13",
            "14","15","16","17","18","19",
            "20","21","22","X");
  return @chrs;
}

sub _get_autosomes_X_chr{
  my @chrs=("chr1","chr2","chr3","chr4","chr5","chr6","chr7",
            "chr8","chr9","chr10","chr11","chr12","chr13",
            "chr14","chr15","chr16","chr17","chr18","chr19",
            "chr20","chr21","chr22","chrX");
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
      my $load_geno=shift;
      my $hchr=$self->_get_hash_autosomes_X_chr();
      open(VCF, $self->{vcffile}) or die "cannot open VCF file\n";
      while (my $line=<VCF>) {
        if ($line =~m/^#/){
             next;
        }else{
          #we have an entry
          my $entry= new SVE($line,$load_geno);
          push (@{$self->{entries}},$entry);
        }
      }

      #function that guest the caller
      my ($caller)=$self->_guess_caller();
      #add missing values to the entries and merge complex ones
      $caller->normalize_sv($self);
      #foreach my $item (@{$self->{entries}}){
      #  print join("\t",$item->{ID},$item->{info}->{PE},$item->{info}->{SR},$item->{info}->{SVTYPE},$item->{info}->{SVLEN})."\n";
      #}
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
     }elsif($item->{info}->{SVMETHOD} =~m/SURVIVOR/){
        $caller="SURVIVOR";
     }elsif(defined $item->{info}->{DBVARID}){
         $caller="GNOMAD";
     }elsif($item->{info}->{SVMETHOD} =~m/PCAWG/){
        $caller="PCAWG";
     }else{
       print "Unknow caller\n";
     }
  }

  #we
  my $tool=();
  if($caller eq "DELLY2"){
    $tool = new DELLY();
    #$tool->normalize_sv($self);
  }elsif($caller eq "Manta"){
    $tool = new Manta();
    #$tool->normalize_sv($self);
  }elsif($caller eq "SVaba"){
    $tool = new SVaba();
  }elsif($caller eq "SURVIVOR"){
    $tool = new SURVIVOR();
  }elsif($caller eq "GNOMAD"){
    $tool = new GNOMAD();
  }elsif($caller eq "PCAWG"){
    $tool = new PCAWG();
  }

  return $tool;
}


sub basic_filters{
    my $self=shift;
    my $rsup=shift;
    my $min_len=shift;
    my $max_len=shift;

    my $hchr=$self->_get_hash_autosomes_X_chr();
    my $remove_by_length=0;
    my $remove_by_chr=0;
    my $remove_by_readsupport=0;
    my $total_variants=0;
    my $alive=0;#pass filter variants
    my $total=0;#total variants

    #remove variants outside
    foreach my $item (@{$self->{entries}}){
      #print Dumper($item);
      $total++;
      my $type=$item->{info}->{SVTYPE};
      #mark the variant as alive
      $item->{info}->{ALIVE}=1;

      if($type ne "BND"){
            if(!defined $hchr->{$item->{CHROM}}){
                    $item->{info}->{ALIVE}=0;
                    $remove_by_chr++;
            }

            if($item->{info}->{SVLEN} < $min_len or $item->{info}->{SVLEN} > $max_len){
                  $remove_by_length++;
                  $item->{info}->{ALIVE}=0;
            }
      }else{
            if(!defined $hchr->{$item->{CHROM}} or !defined $hchr->{$item->{info}->{CHR2}}){
                $item->{info}->{ALIVE}=0;
                $remove_by_chr++;
            }
      }
      #we remove the variants if is supported by less than < $rsup reads
      if($item->{info}->{PE_SR} < $rsup){
          $item->{info}->{ALIVE}=0;
          $remove_by_readsupport++;
      }
      $alive++ if($item->{info}->{ALIVE} > 0);
      #print join("\t",$item->{ID},$item->{CHROM},$item->{POS},$item->{info}->{PE},$item->{info}->{SR},$item->{info}->{SVTYPE},$item->{info}->{SVLEN},$item->{info}->{ALIVE})."\n";
    }

  print "############## BASIC FILTERS START #################\n";
  print "Total SV       : $total\n";
  print "Total Alive SV : $alive\n";
  print "Filter by CHR  : $remove_by_chr\n";
  print "Filter by Length[$min_len,$max_len]: $remove_by_length\n";
  print "Filter by support[<$rsup]: $remove_by_readsupport\n";
  print "############## BASIC FILTERS END #################\n";
}

#print the listed colums to
sub print_matrix{
  my $self=shift;
  my $prefix=shift;

  my @cols=("PON","PON_SUPP","PON_TYPE","PON_IDS","PON_BC1","PON_BC2",
            "GNOMAD","GNOMAD_AC","GNOMAD_TYPE","GNOMAD_IDS","GNOMAD_BC1","GNOMAD_BC2",
            "PCAWG","PCAWG_SUP","PCAWG_TYPE","PCAWG_IDS","PCAWG_BC1","PCAWG_BC2",
            "SOMATIC","SOMATIC_TYPE","SOMATIC_IDS");

  open(FILE,">".$prefix.".txt") or die "cannot open $prefix.txt file\n";
  print FILE join(" ","ID","CHROM","POS","TYPE",@cols)."\n";

  foreach my $item (@{$self->{entries}}){
       next if($item->{info}->{ALIVE} == 0);
       my @tmp=();
        foreach my $v (@cols){
              if(defined $item->{info}->{$v}){
                push(@tmp,$item->{info}->{$v});
              }else{
                push(@tmp,0);
                print "$v not defined on SV info\n";
              }
        }
        print FILE join(" ",$item->{ID},$item->{CHROM},$item->{POS},$item->{info}->{SVTYPE},@tmp)."\n";
  }

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
