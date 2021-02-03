
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
   print "$0 usage : -a <gtf file>  -b <vcf>  -s <tumorID>\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:s:", \%opts );
if ( !defined $opts{a}  ) {
   usage;
}

#load annotation data and store it in the interval tree
my $exons=load_type($opts{a},"exon");
my $t_exon=_build_interval_tree_bed($exons);
my $cds=load_type($opts{a},"CDS");
my $t_cds=_build_interval_tree_bed($cds);
my $genes=load_type($opts{a},"gene");
my $t_gene=_build_interval_tree_bed($genes);


=bla
foreach my $e (@{$cds}[10 .. 20]){
    my $results = $t_cds->fetch($e->{start},$e->{stop});
    #print Dumper($results);
    foreach my $r (@$results){
        my $re=@$cds[$r->{index}];
        print join("\t",$r->{name},$r->{index},$r->{chr},$re->{name},$re->{start},$re->{stop},$re->{tags}->{gene_type},$e->{start},$e->{stop})."\n"
    }
    print "\n";
}
=cut
#we load the VCF file and match the data
#we can load a list or a single file
open(VCF, $opts{b}) or die "cannot open VCF file\n";

print join("\t","TumorID","CHROM","POS",
                      "SVTYPE","PES",
                      "STRANDS","SVLEN",
                      "CHR2","END",
                      "SUPP","CALLERS",
                      "ENS_BRK1", "GN_BRK1", "GT_BRK1", "CDS_BRK1", "EXON_BRK1","CDS_BRK1_D100", "EXON_BRK1_D100",
                      "ENS_BRK2", "GN_BRK2", "GT_BRK1","CDS_BRK2", "EXON_BRK2","CDS_BRK2_D100", "EXON_BRK2_D100",
                      "HIT_GENE","HIT_CDS","HIT_EXON","HIT_CDS_D100","HIT_EXON_D100")."\n";

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
      #this get the breakpoints considering the CIPOS and CIEND post for SVs
      #my ($b1,$b2)=_get_breakpoints($item);
      #delta = 1 for CDS and exons
      my ($cdso1,$cdso2)=overlap_feats($item, $t_cds,1);
      my ($eo1,$eo2)=overlap_feats($item, $t_exon,1);
      #delta = 100
      my ($cdso1_100,$cdso2_100)=overlap_feats($item, $t_cds,100);
      my ($eo1_100,$eo2_100)=overlap_feats($item, $t_exon,100);
      # delta = 1 for genes
      my ($go1,$go2)=overlap_feats($item, $t_gene,1);

      #we parse the results
      #my $coding=0;
      #my $gene=0;
      #print Dumper($eo1,$eo2);
      #print Dumper($eo1_100,$eo2_100);
      #print Dumper($go1,$go2);
      #ENS_BRK1 GN_BRK1 CDS_BRK1 EXON_BRK1 CDS_BRK1_D100 EXON_BRK1_D100 ENS_BRK2 GN_BRK2 CDS_BRK2 EXON_BRK2 CDS_BRK2_D100 EXON_BRK2_D100
      my @v1=();
      my @v2=();
      my $h_cds=0;
      my $h_g=0;
      my $h_e=0;
      my $h_cds_100=0;
      my $h_e_100=0;

      #we check the variables to complete the table
      if(scalar(@$go1) > 0){
          #print Dumper($go1);
          my $hash=get_results($go1,$genes);
          push(@v1,join(";",keys %{$hash->{IDS}}),join(";",keys %{$hash->{GN}}),join(";",keys %{$hash->{GT}}));
          $h_g=1;
      }else{
        push(@v1,"NA","NA","NA");
      }

      if(scalar(@$go2) > 0){
          #print Dumper($go1);
          my $hash=get_results($go2,$genes);
          push(@v2,join(";",keys %{$hash->{IDS}}),join(";",keys %{$hash->{GN}}),join(";",keys %{$hash->{GT}}));
          #push(@v2,join(";",keys %$hash->{IDS}),join(";",keys %$hash->{GN}),join(";",keys %$hash->{GT}));
          $h_g=1;
      }else{
        push(@v2,"NA","NA","NA");
      }

      #CDSs
      if(scalar(@$cdso1) > 0){
          #print Dumper($go1);
          my $hash=get_results($cdso1,$cds);
          push(@v1,join(";",keys %{$hash->{IDS}}));
          $h_cds=1;
      }else{
        push(@v1,"NA");
      }

      if(scalar(@$cdso2) > 0){
          #print Dumper($go1);
          my $hash=get_results($cdso2,$cds);
          push(@v2,join(";",keys %{$hash->{IDS}}));
          $h_cds=1;
      }else{
        push(@v2,"NA");
      }

      #Exons
      if(scalar(@$eo1) > 0){
          #print Dumper($go1);
          my $hash=get_results($eo1,$exons);
          push(@v1,join(";",keys %{$hash->{IDS}}));
          $h_e=1;
      }else{
        push(@v1,"NA");
      }

      if(scalar(@$eo2) > 0){
          #print Dumper($go1);
          my $hash=get_results($eo2,$exons);
          push(@v2,join(";",keys %{$hash->{IDS}}));
          $h_e=1;
      }else{
        push(@v2,"NA");
      }


      #CDSs
      if(scalar(@$cdso1_100) > 0){
          #print Dumper($go1);
          my $hash=get_results($cdso1_100,$cds);
          push(@v1,join(";",keys %{$hash->{IDS}}));
          $h_cds_100=1;
      }else{
        push(@v1,"NA");
      }

      if(scalar(@$cdso2_100) > 0){
          #print Dumper($go1);
          my $hash=get_results($cdso2_100,$cds);
          push(@v2,join(";",keys %{$hash->{IDS}}));
          $h_cds_100=1;
      }else{
        push(@v2,"NA");
      }

      #Exons 100
      if(scalar(@$eo1_100) > 0){
          #print Dumper($go1);
          my $hash=get_results($eo1_100,$exons);
          push(@v1,join(";",keys %{$hash->{IDS}}));
          $h_e_100=1;
      }else{
        push(@v1,"NA");
      }

      if(scalar(@$eo2_100) > 0){
          #print Dumper($go1);
          my $hash=get_results($eo2_100,$exons);
          push(@v2,join(";",keys %{$hash->{IDS}}));
          $h_e_100=1;
      }else{
        push(@v2,"NA");
      }

        print join("\t",$opts{s},$item->{ID},$item->{CHROM},$item->{POS},
                              $item->{info}->{SVTYPE},$item->{info}->{PES},
                              $item->{info}->{STRANDS},$item->{SVLEN},
                              $item->{info}->{CHR2},$item->{info}->{END},
                              $item->{info}->{SUPP}, $item->{info}->{CALLERS},@v1,@v2,$h_g,$h_cds,$h_e,$h_cds_100,$h_e_100)."\n";

      #print Dumper($item);
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





#print Dumper($t_exon);
#$t_exon->print();
=bla
foreach my $e (@{$exons}[10 .. 20]){
    #$t_exon->
    my $results = $t_exon->fetch($e->{start},$e->{stop});
    #print Dumper($results);
    foreach my $r (@$results){
        my $re=@$exons[$r->{index}];
        print join("\t",$r->{name},$r->{index},$r->{chr},$re->{name},$re->{start},$re->{stop},$re->{tags}->{gene_type},$e->{start},$e->{stop})."\n"
    }
    print "\n";
}
print "GENES\n";
foreach my $e (@{$genes}[10 .. 20]){
    #$t_exon->
    my $results = $t_gene->fetch($e->{start},$e->{stop});
    #print Dumper($results);
    foreach my $r (@$results){
        my $re=@$genes[$r->{index}];
        print join("\t",$r->{name},$r->{index},$r->{chr},$re->{name},$re->{start},$re->{stop},$re->{tags}->{gene_type},$e->{start},$e->{stop})."\n"
    }
    print "\n";
}
=cut

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
          if($type eq "CDS"){
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
