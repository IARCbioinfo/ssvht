package SVannot;

=head1 NAME
SVannot

=head1 DESCRIPTION

This object use an interlvalTree to overlap and annoted a give SVobject

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
     print join(" ",$item->{CHROM},$item->{ID},$brk1->{start}-$delta,$brk1->{stop}+$delta)."\n";
     my $rbk1=$tree->fetch($brk1->{start}-$delta,$brk1->{stop}-$delta+$delta);
     print Dumper($rbk1);
     print join(" ",$item->{CHROM},$item->{ID},$brk2->{start}-$delta,$brk2->{stop}+$delta)."\n";
     my $rbk2=$tree->fetch($brk2->{start}-$delta,$brk2->{stop}-$delta+$delta);
     print Dumper($rbk2);
   }

}

#funtion that take into acount the type of variant to build the breakpoint of each SV
sub _get_breakpoints{
      my $item=shift;
      #breakpoint1
      my ($b1s,$b1e)=split(",",$item->{info}->{CIPOS});
      #start/stop for breakpoint 2
      my ($b2s,$b2e)=split(",",$item->{info}->{CIEND});
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


sub _build_interval_tree{
  my $self=shift;
  my $pon=shift; #MERGE from SURVIVOR

  #we create the interval tree for storing the break points of each SV
  my $tree = Set::IntervalTree->new;
  my $i=0;
  foreach my $item (@{$pon->{entries}}){
    #start/stop for breakpoint 1
    #my ($b1s,$b1e)=split(",",$item->{info}->{CIPOS});
    #start/stop for breakpoint 2
    #my ($b2s,$b2e)=split(",",$item->{info}->{CIEND});
    my ($brk1,$brk2)=_get_breakpoints($item);
    my $type=$item->{info}->{SVTYPE};

    #insert breakpoint 1 into the intervalTree
    $tree->insert({chr=>$item->{CHROM},brk=>1,index=>$i,type=>$item->{info}->{SVTYPE},
                id=>$item->{ID}},$brk1->{start},$brk1->{stop});

    #anything that is not a Translocation
    if($type ne "BND" and $type ne "TRA"){
         #print join("\t",$item->{CHROM},$item->{info}->{SVTYPE},$item->{ID},$item->{POS},
          #     join("-",$item->{POS}-abs($b1s),$item->{POS}+abs($b1e)),$item->{info}->{END},
          #     join("-",$item->{info}->{END}-abs($b2s),$item->{info}->{END}+abs($b2e)))."\n";
         #we add the first breakpoint
         #$tree->insert({chr=>$item->{CHROM},brk=>1,index=>$i,type=>$item->{info}->{SVTYPE},
          #           id=>$item->{ID}},$item->{POS}-abs($b1s),$item->{POS}+abs($b1e));
         #we add the second breakpoint
         #$tree->insert({chr=>$item->{CHROM},brk=>2,index=>$i,type=>$item->{info}->{SVTYPE},
          #             id=>$item->{ID}},$item->{info}->{END}-abs($b2s),$item->{info}->{END}+abs($b2e));

          #$tree->insert({chr=>$item->{CHROM},brk=>1,index=>$i,type=>$item->{info}->{SVTYPE},
          #            id=>$item->{ID}},$item->{POS}-abs($b1s),$item->{POS}+abs($b1e));
          #we add the second breakpoint
          $tree->insert({chr=>$item->{CHROM},brk=>2,index=>$i,type=>$item->{info}->{SVTYPE},
                        id=>$item->{ID}},$brk2->{start},$brk2->{stop});

    }else{
       #  print join("\t",$item->{CHROM},$item->{info}->{SVTYPE},$item->{ID},$item->{POS},
       #        join("-",$item->{POS}-$b1s,$item->{POS}+$b1e),$item->{info}->{CHR2},$item->{info}->{POS2},
       #        join("-",$item->{info}->{POS2}-$b2s,$item->{info}->{POS2}+$b2e))."\n";
          #$tree->insert({chr=>$item->{CHROM},brk=>1,index=>$i,type=>$item->{info}->{SVTYPE},
          #           id=>$item->{ID}},$item->{POS}-abs($b1s),$item->{POS}+abs($b1e));
          #we add the second breakpoint
          #$tree->insert({chr=>$item->{info}->{CHR2},brk=>2,index=>$i,type=>$item->{info}->{SVTYPE},
          #             id=>$item->{ID}},$item->{info}->{POS2}-abs($b2s),$item->{info}->{POS2}+abs($b2e));
               $tree->insert({chr=>$item->{info}->{CHR2},brk=>2,index=>$i,type=>$item->{info}->{SVTYPE},
                              id=>$item->{ID}},$brk2->{start},$brk2->{stop});
    }
    $i++;#is an index for fast access to the structure
  }

  return $tree;
}



1;
