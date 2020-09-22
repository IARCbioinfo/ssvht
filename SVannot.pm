package SVannot;

=head1 NAME
SVannot

=head1 DESCRIPTION

This object use an interlvalTree to overlap and annoted a given set of SVs with an SV database.
The SVdatabase might be a panel of normals or a external SV resource like GNOMAD.

=head2 Available methods

=cut

use strict;
use Set::IntervalTree;
use Data::Dumper;

sub new{
  my ($packagename) = @_;
  my $self = {};
  bless ($self, $packagename);
  return ($self);
}

#add missing values to a set of SURVIVOR SVs
sub annot_pon_sv{
   my $self=shift;
   my $pon=shift; #MERGE from SURVIVOR
   my $target=shift; #SOMATICS calls

   my $delta=1000;#delta for breakpoint clustering

   my $tree =$self->_build_interval_tree($pon);
   #we compute the overlaps with the given target calls
   foreach my $item (@{$target->{entries}}){
     my ($brk1,$brk2)=_get_breakpoints($item);
     #add delta to breakpoints
     my $rb1a=$brk1->{start}-$delta/2;
     my $rb1b=$brk1->{stop}+$delta/2;
     #fetch results from the tree
     my $rbk1=$tree->fetch($rb1a,$rb1b);

     print join(" ",$item->{CHROM},$item->{ID},$rb1a,$rb1b,"breakpoint1",abs($rb1b-$rb1a),scalar(@$rbk1))."\n";
     #print Dumper(@$rbk1);
     my $rb2a=$brk2->{start}-$delta/2;
     my $rb2b=$brk2->{stop}+$delta/2;
     my $rbk2=$tree->fetch($rb2a,$rb2b);
     print join(" ",$item->{CHROM},$item->{ID},$rb2a,$rb2b,"breakpoint2",abs($rb2b-$rb2a),scalar(@$rbk2))."\n";

     #DataDumper avoid the printing of equal objects
     #print Dumper($rbk1);
     #print Dumper($rbk2);
     #we filter the breakpoints considering the intervaltree construction, which ignore the chromosome
     my ($fr1,$fr2)=_filter_results_by_chr($item,$rbk1,$rbk2);
     _select_results($item,$fr1,$fr2);
     
   }

}



sub _select_results{
   my ($item,$rbk1,$rbk2)=@_;
   #print join(" ", "count from fucntion",scalar(@$rbk1),scalar(@$rbk2))."\n";
   my $ac=scalar(@$rbk1);
   my $bc=scalar(@$rbk2);
   #print Dumper($rbk1,$rbk2);
   if($ac==0 and $bc==0){
     $item->{info}->{PON}=0;#there is no a matching SV on the panel of normals
     $item->{info}->{PON_IDS}=0;#ids of each matching pon
     $item->{info}->{PON_TYPE}=0;#type of each matching pon
     $item->{info}->{PON_SUPP}=0;#number of genomes supporting the PON
     #print Dumper($item);
     print join(" ","NO-MATCH",$item->{ID},0,0,0,0)."\n";
     return;
   }

   #we have a posible match
   #$tree->insert({chr=>$item->{CHROM},brk=>1,index=>$i,type=>$item->{info}->{SVTYPE},
   #               id=>$item->{ID}},$brk1->{start},$brk1->{stop});

   my $jbk=();
   #first breakpoint
   foreach my $bp(@{$rbk1}){
        $jbk->{$bp->{id}}->{$bp->{brk}}++;
   }
   #second breakpoint
   foreach my $bp(@{$rbk2}){
        $jbk->{$bp->{id}}->{$bp->{brk}}++;
   }

   print Dumper($jbk);
}

sub _filter_by_chr{
    my ($chr,$rbk1)=@_;
    #we filter the list  of breapoin1 by chr
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

  my ($f1)=_filter_by_chr($chr,$rbk1);
  #we change chr in case of a BND or a TRAN for breakpoint 2
  if($type eq "BND" and $type eq "TRA"){
       $chr=$item->{info}->{CHR2};
  }
  #we filter the list of breakpoint2 by chr
  my ($f2)=_filter_by_chr($chr,$rbk2);
  #list of filtered breakpoints
  return ($f1,$f2);
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

#funtion that populate an intervalTree using the GivenSV
sub _build_interval_tree{
  my $self=shift;
  my $pon=shift; #MERGE from SURVIVOR

  #we create the interval tree for storing the break points of each SV
  my $tree = Set::IntervalTree->new;
  my $i=0;

  foreach my $item (@{$pon->{entries}}){
    my ($brk1,$brk2)=_get_breakpoints($item);
    my $type=$item->{info}->{SVTYPE};

    if(abs($brk1->{start}-$brk1->{stop}) ==0 or abs($brk2->{start}-$brk2->{stop}) == 0 ){
      print join("\t",$item->{CHROM},$item->{info}->{SVTYPE},$item->{ID},$item->{POS},
           join("-",$brk1->{start},$brk1->{stop}),$item->{info}->{END},
            join("-",$brk2->{start},$brk2->{stop}))."\n";
            print Dumper($item);
            next;
    }
    #insert breakpoint 1 into the intervalTree
    $tree->insert({chr=>$item->{CHROM},brk=>1,index=>$i,type=>$item->{info}->{SVTYPE},
                   id=>$item->{ID}},$brk1->{start},$brk1->{stop});

    #anything that is not a Translocation
    if($type ne "BND" and $type ne "TRA"){
         #print join("\t",$item->{CHROM},$item->{info}->{SVTYPE},$item->{ID},$item->{POS},
          #     join("-",$item->{POS}-abs($b1s),$item->{POS}+abs($b1e)),$item->{info}->{END},
          #     join("-",$item->{info}->{END}-abs($b2s),$item->{info}->{END}+abs($b2e)))."\n";
          #we add the second breakpoint
          $tree->insert({chr=>$item->{CHROM},brk=>2,index=>$i,type=>$item->{info}->{SVTYPE},
                        id=>$item->{ID}},$brk2->{start},$brk2->{stop});
    }else{
       #  print join("\t",$item->{CHROM},$item->{info}->{SVTYPE},$item->{ID},$item->{POS},
       #        join("-",$item->{POS}-$b1s,$item->{POS}+$b1e),$item->{info}->{CHR2},$item->{info}->{POS2},
       #        join("-",$item->{info}->{POS2}-$b2s,$item->{info}->{POS2}+$b2e))."\n";
          $tree->insert({chr=>$item->{info}->{CHR2},brk=>2,index=>$i,type=>$item->{info}->{SVTYPE},
                          id=>$item->{ID}},$brk2->{start},$brk2->{stop});
    }
    $i++;#is an index for fast access to the structure
  }
  return $tree;
}

1;
