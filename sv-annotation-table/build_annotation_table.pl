
###############################################################################
# Author: Alex Di Genova
# Laboratory: ERABLE/INRIA
# Copyright (c)
# year: 2017
###############################################################################
use Data::Dumper;
use Getopt::Std;
use Set::IntervalTree;
use strict;

sub usage {
   print "$0 usage : -a <gtf file>  -f <fragile_sites> -r <rna_fusions> -b <vcf>  -s <tumorID>\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:s:f:r:", \%opts );
if ( !defined $opts{a}  ) {
   usage;
}

my ($frs)=load_fragile_sites($opts{f});
#print Dumper($frs);
my $t_frs=_build_interval_tree_bed($frs);
#we load the fusions of the sample
my ($f_table,$fbed,$fhash)=load_mRNA_fusions($opts{r},$opts{s});

#function that create the exons/introns with annotations to CDS/UTRs etc
my ($genes,$exons,$introns)=load_exons_genes($opts{a});
my $t_exons=_build_interval_tree_bed($exons);
my $t_introns=_build_interval_tree_bed($introns);
#my $t_fus=_build_interval_tree_bed($fbed);


# we iterate the VCF file to generate the table
open(VCF, $opts{b}) or die "cannot open VCF file\n";

print join("\t","TumorID","SV_ID","CHROM","POS","CIPOS",
                      "SVTYPE","PES",
                      "STRANDS","SVLEN",
                      "CHR2","END","CIEND",
                      "SUPP", "CALLERS",
                      "E1.name","E1.CDS","E1.UTR","E1.chr","E1.start","E1.stop","E1.strand","E1.exon_number","E1.gene_id","E1.gene_name","E1.gene_type",
                      "I1.chr","I1.start","I1.stop","I1.strand","I1.number","I1.gene_id","I1.gene_name",
                      "F1.chr","F1.start","F1.stop","F1.name","F1.Inducer","F1.frequency",
                      "E2.name","E2.CDS","E2.UTR","E2.chr","E2.start","E2.stop","E2.strand","E2.exon_number","E2.gene_id","E2.gene_name","E2.gene_type",
                      "I2.chr","I2.start","I2.stop","I2.strand","I2.number","I2.gene_id","I2.gene_name",
                      "F2.chr","F2.start","F2.stop","F2.name","F2.Inducer","F2.frequency",
                      "GF.name","GF.confidence","GF.read_support","GF.breakpoint1","GF.breakpoint2","GF.closest_genomic_breakpoint1","GF.closest_genomic_breakpoint2","GF.coverage1","GF.coverage2","GF.type")."\n";
#print join("\t",)

while(my $line=<VCF>){
  next if($line=~m/^#/);
  chomp $line;
  my $item=();
  $item->{el}=$line;
  ($item)=_create_entry($item);
  my $type=$item->{info}->{SVTYPE};
  #SURVIVOR pon
  if($type eq "TRA" or $type eq "BND"){
  $item->{info}->{POS2}=$item->{info}->{END};
  #print Dumper($item);
  }

  my $n_exons_b1=0;
  my $n_exons_b2=0;
  my $n_introns_b1=0;
  my $n_introns_b2=0;
  # we match the exons
  my ($exon_o1,$exon_o2)=overlap_feats($item, $t_exons,1);
  # we match the introns
  my ($intron_o1,$intron_o2)=overlap_feats($item, $t_introns,1);
  #we match the fragile sites
  my ($fra_o1,$fra_o2)=overlap_feats($item,$t_frs,1);
  # we match the mRNA fusions
  #my ($rna_o1,$rna_o2)=overlap_feats($item,$t_fus,10000);


  if(scalar(@$exon_o1)>0){
      my ($fe1)=filter_exons_results($exon_o1,$exons);
     $n_exons_b1=scalar(@{$fe1});
  }
  if(scalar(@$exon_o2)>0){
      my ($fe2)=filter_exons_results($exon_o2,$exons);
      $n_exons_b2=scalar(@{$fe2});
  }

  if(scalar(@$exon_o1)==0){
          my ($i1)=filter_introns_results($intron_o1,$introns,$genes);
          #print Dumper($i1);
          $n_introns_b1= defined($i1) ? scalar(@$i1):0;
  }

  if(scalar(@$exon_o2)==0){
          my ($i2)=filter_introns_results($intron_o2,$introns,$genes);
          $n_introns_b2= defined($i2) ? scalar(@$i2):0;
  }

  my @fra_s1=();
  if(scalar(@{$fra_o1}) > 0){
        my ($fra1)=filter_fragile_results($fra_o1,$frs);
        push(@fra_s1,$fra1->{chr},$fra1->{start},$fra1->{stop},$fra1->{name},$fra1->{tags}->{Inducer},$fra1->{tags}->{frequency});
  }else{
    @fra_s1=("NA") x 6;
  }
  my @fra_s2=();
  if(scalar(@{$fra_o2}) > 0){
        my ($fra2)=filter_fragile_results($fra_o2,$frs);
        push(@fra_s2,$fra2->{chr},$fra2->{start},$fra2->{stop},$fra2->{name},$fra2->{tags}->{Inducer},$fra2->{tags}->{frequency});
  }else{
    @fra_s2=("NA") x 6;
  }

  # we print all combinations of exon1--exon2
  if($n_exons_b1 > 0 && $n_exons_b2 > 0){
    my ($fe1)=filter_exons_results($exon_o1,$exons);
    my ($fe2)=filter_exons_results($exon_o2,$exons);
    foreach my $e1 (@{$fe1}){
      my $re1=print_exon($e1);
      foreach my $e2 (@{$fe2}){
          my $re2=print_exon($e2);
          my $p_fus=check_fusion($e1->{tags}->{gene_name},$e2->{tags}->{gene_name},$fhash,$f_table);
          print join("\t",$opts{s},$item->{ID},$item->{CHROM},$item->{POS},$item->{info}->{CIPOS},
                                $item->{info}->{SVTYPE},$item->{info}->{PES},
                                $item->{info}->{STRANDS},$item->{info}->{SVLEN},
                                $item->{info}->{CHR2},$item->{info}->{END},$item->{info}->{CIEND},
                                $item->{info}->{SUPP}, $item->{info}->{CALLERS},
                                @$re1,("NA") x 7,@fra_s1,
                                @$re2,("NA") x 7,@fra_s2,
                                @$p_fus)."\n";
        }

    }


  }elsif($n_exons_b1 > 0 && $n_introns_b2 > 0){
    #exon1 -- intron2
    my ($fe1)=filter_exons_results($exon_o1,$exons);
    my ($fi2)=filter_introns_results($intron_o2,$introns,$genes);
    foreach my $e1 (@{$fe1}){
      my $re1=print_exon($e1);
      foreach my $i2 (@{$fi2}){
          my $ri2=print_intron($i2);
          my $p_fus=check_fusion($e1->{tags}->{gene_name},$i2->{gene_name},$fhash,$f_table);
          print join("\t",$opts{s},$item->{ID},$item->{CHROM},$item->{POS},$item->{info}->{CIPOS},
                                $item->{info}->{SVTYPE},$item->{info}->{PES},
                                $item->{info}->{STRANDS},$item->{info}->{SVLEN},
                                $item->{info}->{CHR2},$item->{info}->{END},$item->{info}->{CIEND},
                                $item->{info}->{SUPP}, $item->{info}->{CALLERS},
                                @$re1,("NA") x 7,@fra_s1,
                                ("NA") x 11,@$ri2,@fra_s2,
                                @$p_fus)."\n";
        }

    }


    #intron1 -- #exon2
  }elsif($n_exons_b2 > 0 && $n_introns_b1 > 0){
    #intron1 -- intron2
    my ($fi1)=filter_introns_results($intron_o1,$introns,$genes);
    my ($fe2)=filter_exons_results($exon_o2,$exons);
    foreach my $i1 (@{$fi1}){
      my $ri1=print_intron($i1);
      foreach my $e2 (@{$fe2}){
          my $re2=print_exon($e2);
          my $p_fus=check_fusion($i1->{gene_name},$e2->{tags}->{gene_name},$fhash,$f_table);
          print join("\t",$opts{s},$item->{ID},$item->{CHROM},$item->{POS},$item->{info}->{CIPOS},
                                $item->{info}->{SVTYPE},$item->{info}->{PES},
                                $item->{info}->{STRANDS},$item->{info}->{SVLEN},
                                $item->{info}->{CHR2},$item->{info}->{END},$item->{info}->{CIEND},
                                $item->{info}->{SUPP}, $item->{info}->{CALLERS},
                                ("NA") x 11,@$ri1,@fra_s1,
                                @$re2,("NA") x 7,@fra_s2,
                                @$p_fus)."\n";
        }

    }


  }elsif($n_introns_b1 > 0 && $n_introns_b2 > 0){
    my ($fi1)=filter_introns_results($intron_o1,$introns,$genes);
    my ($fi2)=filter_introns_results($intron_o2,$introns,$genes);
    foreach my $i1 (@{$fi1}){
      my $ri1=print_intron($i1);
      foreach my $i2 (@{$fi2}){
        my $ri2=print_intron($i2);
        my $p_fus=check_fusion($i1->{gene_name},$i2->{gene_name},$fhash,$f_table);
        #we print the results
        print join("\t",$opts{s},$item->{ID},$item->{CHROM},$item->{POS},$item->{info}->{CIPOS},
                              $item->{info}->{SVTYPE},$item->{info}->{PES},
                              $item->{info}->{STRANDS},$item->{info}->{SVLEN},
                              $item->{info}->{CHR2},$item->{info}->{END},$item->{info}->{CIEND},
                              $item->{info}->{SUPP}, $item->{info}->{CALLERS},
                              ("NA") x 11,@$ri1,@fra_s1,
                              ("NA") x 11,@$ri2,@fra_s2,
                              @$p_fus)."\n";
      }
   }
  }elsif($n_exons_b1 > 0){
    my ($fe1)=filter_exons_results($exon_o1,$exons);
    foreach my $e1 (@{$fe1}){
      my $re1=print_exon($e1);
      print join("\t",$opts{s},$item->{ID},$item->{CHROM},$item->{POS},$item->{info}->{CIPOS},
                      $item->{info}->{SVTYPE},$item->{info}->{PES},
                      $item->{info}->{STRANDS},$item->{info}->{SVLEN},
                      $item->{info}->{CHR2},$item->{info}->{END},$item->{info}->{CIEND},
                      $item->{info}->{SUPP}, $item->{info}->{CALLERS},
                      @$re1,("NA") x 7,@fra_s1,
                      ("NA") x 11,("NA") x 7,@fra_s2,
                      ("NA") x 10)."\n";
    }


  }elsif($n_exons_b2 > 0){
    my ($fe2)=filter_exons_results($exon_o2,$exons);
  foreach my $e2 (@{$fe2}){
      my $re2=print_exon($e2);
      print join("\t",$opts{s},$item->{ID},$item->{CHROM},$item->{POS},$item->{info}->{CIPOS},
                      $item->{info}->{SVTYPE},$item->{info}->{PES},
                      $item->{info}->{STRANDS},$item->{info}->{SVLEN},
                      $item->{info}->{CHR2},$item->{info}->{END},$item->{info}->{CIEND},
                      $item->{info}->{SUPP}, $item->{info}->{CALLERS},
                      ("NA") x 11,("NA") x 7,@fra_s1,
                      @$re2,("NA") x 7,@fra_s2,
                      ("NA") x 10)."\n";
  }

  }elsif($n_introns_b1 > 0){
    my ($fi1)=filter_introns_results($intron_o1,$introns,$genes);
    foreach my $i1 (@{$fi1}){
      my $ri1=print_intron($i1);
        #we print the results
        print join("\t",$opts{s},$item->{ID},$item->{CHROM},$item->{POS},$item->{info}->{CIPOS},
                              $item->{info}->{SVTYPE},$item->{info}->{PES},
                              $item->{info}->{STRANDS},$item->{info}->{SVLEN},
                              $item->{info}->{CHR2},$item->{info}->{END},$item->{info}->{CIEND},
                              $item->{info}->{SUPP}, $item->{info}->{CALLERS},
                              ("NA") x 11,@$ri1,@fra_s1,
                              ("NA") x 11,("NA") x 7,@fra_s2,
                              ("NA") x 10)."\n";
   }



  }elsif($n_introns_b2 > 0){
    my ($fi2)=filter_introns_results($intron_o2,$introns,$genes);
      foreach my $i2 (@{$fi2}){
        my $ri2=print_intron($i2);
        print join("\t",$opts{s},$item->{ID},$item->{CHROM},$item->{POS},$item->{info}->{CIPOS},
                              $item->{info}->{SVTYPE},$item->{info}->{PES},
                              $item->{info}->{STRANDS},$item->{info}->{SVLEN},
                              $item->{info}->{CHR2},$item->{info}->{END},$item->{info}->{CIEND},
                              $item->{info}->{SUPP}, $item->{info}->{CALLERS},
                              ("NA") x 11,("NA") x 7,@fra_s1,
                              ("NA") x 11,@$ri2,@fra_s2,
                              ("NA") x 10)."\n";
    }
  }else{
    print join("\t",$opts{s},$item->{ID},$item->{CHROM},$item->{POS},$item->{info}->{CIPOS},
                    $item->{info}->{SVTYPE},$item->{info}->{PES},
                    $item->{info}->{STRANDS},$item->{info}->{SVLEN},
                    $item->{info}->{CHR2},$item->{info}->{END},$item->{info}->{CIEND},
                    $item->{info}->{SUPP}, $item->{info}->{CALLERS},
                    ("NA") x 11,("NA") x 7,@fra_s1,
                    ("NA") x 11,("NA") x 7,@fra_s2,
                    ("NA") x 10)."\n";
  }

  #print join("\t",$opts{s},$item->{ID},$n_exons_b1,$n_exons_b2,$n_introns_b1,$n_introns_b2)."\n";

}

sub print_intron{
  my ($i)=@_;
  my $re1=();
  push(@$re1,$i->{chr},$i->{start},$i->{stop},$i->{strand},$i->{number},$i->{gene_id},$i->{gene_name});
  return $re1;
}

sub check_fusion{
  my ($g1,$g2,$fhash,$f_table)=@_;
  my $f_index=-1;
  if(defined $fhash->{join("__",$g1,$g2)}){
      $f_index=$fhash->{join("__",$g1,$g2)};
  }elsif(defined $fhash->{join("__",$g2,$g1)}){
      $f_index=$fhash->{join("__",$g2,$g1)};
  }

  my $p_fus=();
  if($f_index >=0){
    my $f=@$f_table[$f_index];
    push(@$p_fus,$f->{name},$f->{confidence},$f->{read_support},$f->{breakpoint1},$f->{breakpoint2},
                ,$f->{closest_genomic_breakpoint1},$f->{closest_genomic_breakpoint2},$f->{coverage1},
                $f->{coverage2},$f->{type});
  }else{
    @$p_fus=("NA") x 10;
  }
  return $p_fus;
}


sub print_exon{
  my ($e)=@_;
  my $re1=();
  push(@$re1,$e->{name},$e->{CDS},$e->{UTR},
       $e->{chr},$e->{start},$e->{stop},$e->{strand},
       $e->{tags}->{exon_number},$e->{tags}->{gene_id},$e->{tags}->{gene_name},$e->{tags}->{gene_type});
   return $re1;
}


sub filter_fragile_results{
  my ($res, $db)=@_;
  my @hash=();
  my $index=0;
  my $i=0;
  foreach my $r (@$res){
        my $re=@$db[$r->{index}];
        #print Dumper($re);
        push(@hash,$re);
        if($re->{tags}->{frequency} =~m/common/i){
          $index=$i;
        }
        $i++;
  }
  #we return the most frequence fragile site
  if($i > 1){
    return($hash[$index]);
  }else{
    return ($hash[0]);
  }

}
#we filter the introns
sub filter_introns_results{
  my ($res, $db)=@_;
  my $hash=();
  foreach my $r (@$res){
        my $re=@$db[$r->{index}];
        #print Dumper($re);
        push(@$hash,$re);
  }
  return $hash;
}
#we filter the exons from the same gene
sub filter_exons_results{
    my ($res, $db)=@_;
    my $hash=();
    foreach my $r (@$res){
        my $re=@$db[$r->{index}];
        #we store all the exons from the same gene
        if(!defined $hash->{$re->{tags}->{gene_id}}){
          $hash->{$re->{tags}->{gene_id}}=$re;
        }elsif($re->{CDS} >=1 and $hash->{$re->{tags}->{gene_id}}->{CDS} == 0){
          $hash->{$re->{tags}->{gene_id}}=$re;
        }elsif($re->{UTR} >=1 and $hash->{$re->{tags}->{gene_id}}->{UTR} == 0){
          $hash->{$re->{tags}->{gene_id}}=$re;
        }
    }
    #print Dumper($hash);
    my $f=();
    foreach my $re (keys %{$hash}){
        push (@$f,$hash->{$re});
    }
    return $f;
}


#We load the fusions from arriba and
sub load_mRNA_fusions{
      my ($file,$s)=@_;
	open(FILE,$file) or die "cannot open file $file\n";
	#we read the header
	my $h=<FILE>;
	chomp $h;
	my @values=split("\t",$h);
	my $index=0;
	my $fbed=();
	my $fusions=();
  my $fhash=();
	while(my $line=<FILE>){
	    chomp $line;
            my @data=split /\t/,$line;
	    next if($data[1] ne $s);
	    my $tmp=();
	    for(my $i=0;$i<$#values; $i++){
		$tmp->{$values[$i]}=$data[$i];
	    }
	    $tmp->{read_support}=$tmp->{split_reads1}+$tmp->{split_reads2}+$tmp->{discordant_mates};
	    $tmp->{index}=$index;
	    #we drops some features
	    $tmp->{fusion_transcript}="";
	    $tmp->{peptide_sequence}="";
	    #$tmp->{name}=join("__",$tmp->{})
	    $tmp->{name}=join("__",$tmp->{gene1},$tmp->{gene2});
	    my $tmp_b1=();
	       my ($a,$b)=split(":",$tmp->{breakpoint1});
		  $tmp_b1->{chr}="chr".$a;
		  $tmp_b1->{start}=$b;
		  $tmp_b1->{stop}=$b+1;
		  $tmp_b1->{index}=$index;
		  #$tmp_b1->{name}=join("__",$tmp->{gene1},$tmp->{gene2});
		  $tmp_b1->{name}=$tmp->{gene1};
	        #we split the breakpoint2
		($a,$b)=split(":",$tmp->{breakpoint2});
	    my $tmp_b2=();
		  $tmp_b2->{chr}="chr".$a;
		  $tmp_b2->{start}=$b;
		  $tmp_b2->{stop}=$b+1;
		  $tmp_b2->{index}=$index;
		  #$tmp_b2->{name}=join("__",$tmp->{gene1},$tmp->{gene2});
		  $tmp_b2->{name}=$tmp->{gene2};
      #we save a hash with gene1 and gene2 to make overlap by gene names rather than fusions
      if(!defined $fhash->{join("__",$tmp->{gene1},$tmp->{gene2})}){
          $fhash->{join("__",$tmp->{gene1},$tmp->{gene2})}=$index;
      }
		push(@$fusions, $tmp);
		push(@$fbed,$tmp_b1);
		push(@$fbed,$tmp_b2);
	  $index++;
	}
  #return the fusions, the bed and a hash with fhash
	return ($fusions,$fbed,$fhash);
}


sub load_fragile_sites{
    my ($sites)=@_;
    open(FILE,$sites) or die "cannot open $sites\n";
    my $frs=();
    my $i=0;
    #we load the fragile sites
    while(my $line=<FILE>){
      chomp $line;
      my @d=split("\t",$line);
      my $tmp=();
      ($tmp->{chr},$tmp->{start},$tmp->{stop},
      $tmp->{strand},$tmp->{type}) = ($d[0],$d[3],$d[4],$d[6],$d[2]);
      foreach my $t(split(";",$d[8])){
        my ($v,$v2)=split("=",$t);
        $tmp->{tags}->{$v}=$v2;
      }
      $tmp->{name}=$tmp->{tags}->{name};
      if(!defined $tmp->{tags}->{Inducer}){
        $tmp->{tags}->{Inducer}="NA";
      }
      $tmp->{index}=$i;
      if($tmp->{chr} eq "chrx"){
        $tmp->{chr}="chrX";
      }
      $i++;
      #there is a bad line in the file
      #if($tmp->{start} > $tmp->{stop}){
        #print $line."\n";
        #print Dumper(@d);
        #print $tmp->{start}." ".$tmp->{stop}."\n";
      #}
      push(@{$frs},$tmp);
    }
    return ($frs);
}


sub get_results{
    my ($res, $db)=@_;
    my $hash=();
    foreach my $r (@$res){
        my $re=@$db[$r->{index}];
        #print join("\t",$r->{name},$r->{index},$r->{chr},$re->{name},$re->{start},$re->{stop},$re->{tags}->{gene_type})."\n"
        $hash->{IDS}->{$re->{name}}++;
        if($re->{type} eq "gene"){
          $hash->{GN}->{$re->{tags}->{"gene_name"}}++;
          $hash->{GT}->{$re->{tags}->{"gene_type"}}++;
        }
        #print Dumper($re);
    }
    #print Dumper($hash);
    return $hash;
}



sub overlap_feats{
   my ($item, $db, $delta)=@_;
   # we build the breakpoints
   my ($b1,$b2)=_get_breakpoints($item);
   #we add an additional delta def is 1
   my $rbk1=$db->fetch($b1->{start}-$delta,$b1->{stop}+$delta);
   my $rbk2=$db->fetch($b2->{start}-$delta,$b2->{stop}+$delta);
   #we filter results by chromosome
  my ($fr1,$fr2)=_filter_results_by_chr($item,$rbk1,$rbk2);

}




#filter results by chromosome, INV,INS,DEL,DUP,TRA
sub _filter_by_chr{
    my ($chr,$rbk1)=@_;
    #we filter the list  of breapoints by chr
    my $f1=[];#array of filtered breakpoints1
    foreach my $r(@{$rbk1}){
        if($r->{chr} eq $chr){
             push(@{$f1},$r);
          }
      }
      return $f1;
}


#filter breakpoint results by chr, cause of how the interval tree is build
sub _filter_results_by_chr{
  my ($item,$rbk1,$rbk2)=@_;
  my $f1=[];#array of filtered breakpoints1
  my $f2=[]; #array of filtered breakpoint2
  my $type=$item->{info}->{SVTYPE};
  my $chr=$item->{CHROM};#chr for breakpoint 1;

  ($f1)=_filter_by_chr($chr,$rbk1);
  #we change chr in case of a BND or a TRAN for breakpoint 2
  if($type eq "BND" or $type eq "TRA"){
       $chr=$item->{info}->{CHR2};
  }
  #we filter the list of breakpoint2 by chr
  ($f2)=_filter_by_chr($chr,$rbk2);
  #list of filtered breakpoints
  return ($f1,$f2);
}



#funtions to read an entry from a VCF file
sub _create_entry{
      my $item=shift;
      my $tags=_get_tags($item);
      $item->{info}=$tags;
      #$item->{geno}=$geno;
      my @data=split("\t",$item->{el});
      $item->{CHROM}=$data[0];
      $item->{POS}=$data[1];
      $item->{ID}=$data[2];
      $item->{REF}=$data[3];
      $item->{ALT}=$data[4];
      $item->{QUAL}=$data[5];
      $item->{FILTER}=$data[6];

      return $item;
}



#parse the tags of the tool
sub _get_tags{
    my $item=shift;
    my $tags=();
    #chomp $entry;
    my @d=split("\t",$item->{el});
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




#funtion that take into acount the type of variant to build the breakpoint of each SV
sub _get_breakpoints{
      my $item=shift;
      #breakpoint1
      my ($b1s,$b1e)=split(",",$item->{info}->{CIPOS});
      $b1s++ if($b1s ==0);
      $b1e++ if($b1e ==0);
      #start/stop for breakpoint 2
      my ($b2s,$b2e)=split(",",$item->{info}->{CIEND});
      $b2s++ if($b2s ==0);
      $b2e++ if($b2e ==0);
      #get the type of the current SV
      my $type=$item->{info}->{SVTYPE};
      my $brk1=();
      my $brk2=();
      #breakpoint 1 boundaries
      $brk1->{start}=$item->{POS}-abs($b1s);
      $brk1->{stop}=$item->{POS}+abs($b1e);
      if($type ne "BND" and $type ne "TRA"){
          #breakpoint 2 boundaries
          $brk2->{start}=$item->{info}->{END}-abs($b2s);
          $brk2->{stop}=$item->{info}->{END}+abs($b2e);
      }else{
        #breakpoint 2 boundaries
        $brk2->{start}=$item->{info}->{POS2}-abs($b2s);
        $brk2->{stop}=$item->{info}->{POS2}+abs($b2e);
      }
      return ($brk1,$brk2);
}


#load the exons from the GTF files marking whether an exon is CDS, UTR and its associated gene
sub load_exons_genes{
    my ($file,$type)=@_;
    open(GTF, "gzip -dc $file |") or die "cannot open file\n";
    my $exons=(); #exons
    my $genes=(); #gene features
    my $exon2genes=(); #exon 2 gene
    my $introns=();

    while (my $line=<GTF>) {
            next if($line=~m/^#/);
            chomp $line;
            my @d=split("\t",$line);
            #we get the features and tags for each value
            my $e=();
            ($e->{chr},$e->{start},$e->{stop},
            $e->{strand},$e->{type}) = ($d[0],$d[3],$d[4],$d[6],$d[2]);
            #we parse the tags of the file
            $e->{tags}=parse_tags($d[8]);
            $e->{name}=$e->{tags}->{join("_",$d[2],"id")};
            if($d[2] eq "CDS" or $d[2] eq "UTR"){
              $e->{name}=$e->{tags}->{join("_","exon","id")};
            }
            #we skyp non exonic genes
            if($d[2] eq "gene"){
              $genes->{$e->{name}}=$e; #we save the gene information
            }elsif($d[2] eq "exon"){
              if(!defined $exons->{$e->{name}}){
                $e->{CDS}=0;
                $e->{UTR}=0;
                $exons->{$e->{name}}=$e;

              }
            }elsif($d[2] eq "CDS"){
              if(defined $exons->{$e->{name}}){
                $exons->{$e->{name}}->{CDS}++;
              }else{
                print STDERR $e->{name}. "is not present on exon db\n";
              }
            }elsif($d[2] eq "UTR"){
              if(defined $exons->{$e->{name}}){
                $exons->{$e->{name}}->{UTR}++;
              }else{
                print STDERR $e->{name}. "is not present on exon db\n";
              }
            }
    }
    #print Dumper($exons);
    #we sort the exons by chromosome, gene and position
    my @sexons= sort {$exons->{$a}->{chr} cmp $exons->{$b}->{chr} ||
                      $exons->{$a}->{tags}->{gene_id} cmp $exons->{$b}->{tags}->{gene_id} ||
                      $exons->{$a}->{start} <=> $exons->{$b}->{start} } keys %{$exons};

    #print Dumper(@sexons);
    my $index_exon=0;
    my $index_intron=0;
    my $sorted_exons=();
    my $gname=$exons->{$sexons[0]}->{tags}->{gene_id};
    my $uniq_exons=();
    foreach my $n(@sexons){
      my $e=$exons->{$n};
      my $g=$genes->{$e->{tags}->{gene_id}};#we get the gene
      my $gl=abs($g->{start}-$g->{stop});
      if($gname eq $g->{tags}->{gene_id}){
        push(@$uniq_exons,$e);
      }else{
        #print join(" ",$gname,$g->{tags}->{gene_id});

        my ($intron)=build_introns($uniq_exons,$genes);
        if(defined $intron){
            foreach my $in(@{$intron}){
              $in->{index}=$index_intron;
              $index_intron++;
              push(@$introns,$in);
            }
        }else{
          #genes with no introns
          #print Dumper($intron);
          #print STDERR $gname." do not have introns\n";
        }
        $gname=$g->{tags}->{gene_id};
        $uniq_exons=();
        push(@$uniq_exons,$e);
      }
      #my $g=$genes->{$e->{tags}->{gene_id}};#we get the gene
      #print join(" ",$e->{chr},$e->{name},$e->{tags}->{gene_id},$e->{CDS},$e->{UTR},$e->{start},$e->{stop},$e->{strand},
      #$g->{tags}->{gene_name},$g->{start},$g->{stop},$g->{strand},$g->{tags}->{gene_type})."\n";
      $e->{index}=$index_exon;
      $index_exon++;

      push(@$sorted_exons,$e);
    }


    #we add the index to exons and introns
    #we create the introns for each gene
    #print Dumper($introns);
    return ($genes,$sorted_exons,$introns);

    #exit 0;
    #we return the coordiantes for each gene

}

sub build_introns{
    my ($uniq_e,$genes)=@_;
    my $g=$genes->{@$uniq_e[0]->{tags}->{gene_id}};
    my $l=abs($g->{start}-$g->{stop});
    #we init the array with N
    my @cigar=("N") x $l;
    #print Dumper($uniq_e);
    #print Dumper(@g_space);
    #we mark all the exon sequence on the gene regions
    foreach my $e(@$uniq_e){
      my $se=$e->{start}-$g->{start};
      my $ee=abs($e->{start}-$e->{stop});
      #print join(" ",$se,$ee,$l)."\n";
      for(my $i=$se;$i<=$se+$ee;$i++){
        $cigar[$i]="E";
      }
    }
    #we build the introns start/stop sequences
    my $first=-1;
    my $last=0;
    #we assume that the gene start by a E
    if($cigar[0] eq "N"){
      print STDERR $g->{tags}->{gene_id}." Start with N\n";
    }else{
      #we have observed that is due to exons wiht start=stop
      $cigar[0]="E";
    }
    #We check some conditions
    if($cigar[$#cigar] eq "N"){
      print STDERR $g->{tags}->{gene_id}." End with N\n";
    }
    my $introns=();
    my $n=1;
    for(my $i=1; $i <$l ; $i++){
        if($cigar[$i] eq "E" and $cigar[$i-1] eq "N"){
          my $tmp=();
            $tmp->{start}=$g->{start}+$first+1;
            #case of micro-introns of lenght 1 e.g ENE
            if($last == 0){
                $tmp->{stop}=$g->{start}+$first+1;
            }else{
                $tmp->{stop}=$g->{start}+$first+$last;
           }
            $tmp->{gene_id}=$g->{tags}->{gene_id};
            $tmp->{chr}=$g->{chr};
            $tmp->{strand}=$g->{strand};
            $tmp->{number}=$n;
            $tmp->{name}=join("_",$g->{tags}->{gene_id},$n);
            $tmp->{gene_name}=$g->{tags}->{gene_name};
            push(@{$introns},$tmp);
            $n++;
            #debug micro-introns
          # if($tmp->{start} > $tmp->{stop}){
          #      print join("",@cigar)."\n";
          #      foreach(@{$introns}){
          #          print join(" ",$_->{start},$_->{stop},$i,$first,$last,$l,$_->{g_id})."\n";
          #     }
          # }
            #print join(" ","INTRON",$first,$last)."\n";
        $first=-1; $last=0;
      }elsif($cigar[$i] eq "N" and $first== -1){
          $first=$i;
        }elsif($cigar[$i] eq "N"){
          $last++;
        }
    }
    #print join("",@cigar)."\n";
    #print Dumper($introns);
    return ($introns);
}

#load different types of features from GTF file
sub load_type{
        my ($file,$type)=@_;
  open(GTF, "gzip -dc $file |") or die "cannot open file\n";
  #exonic regions of protein coding genes
  my $exons=();
  my $i=0;
  my $t_skyp=0;
  #my $t="";
  while (my $line=<GTF>) {
          next if($line=~m/^#/);
          chomp $line;
          my @d=split("\t",$line);
          #we skyp non exonic genes
          next if($d[2] ne $type);
          #print Dumper(@d);
          my $e=();
          ($e->{chr},$e->{start},$e->{stop},
          $e->{strand},$e->{type}) = ($d[0],$d[3],$d[4],$d[6],$d[2]);
          #we parse the tags of the file
          $e->{tags}=parse_tags($d[8]);
          $e->{name}=$e->{tags}->{join("_",$type,"id")};
          if($type eq "CDS" or $type eq "UTR"){
            $e->{name}=$e->{tags}->{join("_","exon","id")};
          }
          #just a warning
          if(abs($e->{start}-$e->{stop}) == 0){
            #print STDERR $e->{name}." skip exon, start and stop are equal\n";
            $t_skyp++;
            #$t=$e->{type};
          }else{
          $e->{index}=$i;
          push(@{$exons},$e);
          $i++;
          }
          #sprint Dumper($e);
  }
  close(GTF);
  #if($t_skyp){
    print STDERR "A total of $t_skyp $type were skyped [start == stop]\n";
  #}
  return $exons;
}



#build the intervaltree from a bedfile
sub  _build_interval_tree_bed{
my ($bed_file)=@_;
#print Dumper($bed_file);
#we create the interval tree for storing the break points of each SV
my $tree = Set::IntervalTree->new;
foreach my $e (@{$bed_file}){
    if(abs($e->{start}-$e->{stop})==0){
        #print STDERR  Dumper($e);
        print STDERR "start equall to stop\n";
        #$line." start equal to stop\n";
      }else{
        #print join(" ",$e->{name},$e->{start},$e->{stop})."\n";
        $tree->insert({chr=>$e->{chr},name=>$e->{name},index=>$e->{index}},$e->{start},$e->{stop});
      }
 }
 #$tree->print();
 #print Dumper($tree);
  return ($tree);
}

#we parse the tags
sub parse_tags{
    my ($t)=@_;
    my @vv=split (";",$t);
    my $tags=();
    #print Dumper(@vv);
    foreach my $d (@vv){
        my ($n,$v)=split(" ",$d);
        $v=~s/\"//g;
        #print join(" ",$n,$v)."\n";
        $tags->{$n}=$v;
    }
    return $tags;
}
