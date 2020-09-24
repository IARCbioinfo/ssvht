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
   my $type=shift; #consider type of variant
   my $delta=shift;#delta around bearkpoints

   #my $delta=1000;#delta for breakpoint clustering

   my $tree =$self->_build_interval_tree($pon);
   my $total_vars=0;
   my $total_annotations=0;
   #we compute the overlaps with the given target calls
   foreach my $item (@{$target->{entries}}){
     next if($item->{info}->{ALIVE} == 0);
     my ($brk1,$brk2)=_get_breakpoints($item);
     $total_vars++;
     #add delta to breakpoints
     my $rb1a=$brk1->{start}-$delta/2;
     my $rb1b=$brk1->{stop}+$delta/2;
     #fetch results from the tree
     my $rbk1=$tree->fetch($rb1a,$rb1b);
    #print Dumper(@$rbk1);
     my $rb2a=$brk2->{start}-$delta/2;
     my $rb2b=$brk2->{stop}+$delta/2;
     #fecth the result from the tree for the second breakpoint
     my $rbk2=$tree->fetch($rb2a,$rb2b);
    #DataDumper avoid the printing of equal objects
     #print Dumper($rbk1);
     #print Dumper($rbk2);
     #we filter the breakpoints considering the intervaltree construction, which ignore the chromosome
     my ($fr1,$fr2)=_filter_results_by_chr($item,$rbk1,$rbk2);
     #we filter by SVTYPE
     if($type == 1){
        $fr1=_filter_by_svtype($item->{info}->{SVTYPE},$fr1);
        $fr2=_filter_by_svtype($item->{info}->{SVTYPE},$fr2);
     }

     #print join(" ",$item->{CHROM},$item->{ID},$rb1a,$rb1b,"breakpoint1",abs($rb1b-$rb1a),scalar(@$rbk1), $item->{info}->{SVLEN})."\n";
     #print join(" ",$item->{CHROM},$item->{ID},$rb2a,$rb2b,"breakpoint2",abs($rb2b-$rb2a),scalar(@$rbk2), $item->{info}->{SVLEN})."\n";
     #we rank selected results
     my ($matches)=_select_matchs($item,$fr1,$fr2);

     my $nmatches=scalar(@{$matches});
     #the variant do not match PON
     if($nmatches == 0){
       $item->{info}->{PON}=0;#there is no a matching SV on the panel of normals
       $item->{info}->{PON_IDS}=0;#ids of each matching pon
       $item->{info}->{PON_TYPE}=0;#type of each matching pon
       $item->{info}->{PON_SUPP}=0;#number of genomes supporting the PON
     }else{
        my $tmp=();

        foreach my $r (@{$matches}){
            my $p=@{$pon->{entries}}[$r->{index}];
            push(@{$tmp->{PON_SUPP}},$p->{info}->{SUPP});
            push(@{$tmp->{PON_TYPE}},$p->{info}->{SVTYPE});
            push(@{$tmp->{PON_IDS}},$r->{name});
            #print Dumper($r);
            #print Dumper($p);
        }
        $total_annotations++;
        #we fill the info
        $item->{info}->{PON}=$nmatches;#there is $nmatches matching the PON
        $item->{info}->{PON_SUPP}=join(",",@{$tmp->{PON_SUPP}});
        $item->{info}->{PON_TYPE}=join(",",@{$tmp->{PON_TYPE}});
        $item->{info}->{PON_IDS}=join(",",@{$tmp->{PON_IDS}});
     }
     #print join(" ",$item->{CHROM},$item->{ID}, $item->{info}->{PON},
     #                $item->{info}->{PON_SUPP},$item->{info}->{PON_TYPE},
     #                $item->{info}->{PON_IDS})."\n";
     #print join(" ",$item->{CHROM},$item->{ID},scalar(@{$matches}),$nmatches)."\n";
   }
   print "Total SVs : $total_vars\nTotal Annotated PON : $total_annotations\nPON annotated (\%) : ".$total_annotations/$total_vars."\n";
}


#funtion that select the right results
sub _select_matchs{
   my ($item,$rbk1,$rbk2)=@_;

   #$tree->insert({chr=>$item->{CHROM},brk=>1,index=>$i,type=>$item->{info}->{SVTYPE},
   #               id=>$item->{ID}},$brk1->{start},$brk1->{stop});

   my $jbk=();
   my $both_bp_match=[];
   #first breakpoint list

   foreach my $bp(@{$rbk1}){
        $jbk->{$bp->{index}}->{bpd}->{$bp->{brk}}++; #check the number of overlapping
        $jbk->{$bp->{index}}->{bpq}->{1}++;
        $jbk->{$bp->{index}}->{name}=$bp->{id}; #the name of the SVs might be the same for the PON
        $jbk->{$bp->{index}}->{type}=$bp->{type};
        $jbk->{$bp->{index}}->{chr}=$bp->{chr};
        $jbk->{$bp->{index}}->{len}=$bp->{len};
        $jbk->{$bp->{index}}->{index}=$bp->{index};
        #print join(" ",$bp->{chr},$bp->{brk},$bp->{index},$bp->{type},$bp->{id})."\n";
   }
   #second breakpoint list
   foreach my $bp(@{$rbk2}){
        #$jbk->{$bp->{id}}->{$bp->{brk}}++;
        $jbk->{$bp->{index}}->{bpd}->{$bp->{brk}}++;
        $jbk->{$bp->{index}}->{bpq}->{2}++;
        $jbk->{$bp->{index}}->{name}=$bp->{id};
        $jbk->{$bp->{index}}->{type}=$bp->{type};
        $jbk->{$bp->{index}}->{chr}=$bp->{chr};
        $jbk->{$bp->{index}}->{len}=$bp->{len};
        $jbk->{$bp->{index}}->{index}=$bp->{index};
        #print join(" ",$bp->{chr},$bp->{brk},$bp->{index},$bp->{type},$bp->{id})."\n";
   }

   #print Dumper($jbk);
   #we iterate the jbk looking for hits with two breakpoints
   foreach my $id (keys %{$jbk}){
     my ($bpq)=scalar(keys %{$jbk->{$id}->{bpq}});
     my ($bpd)=scalar(keys %{$jbk->{$id}->{bpd}});
     if($bpd == 2 and $bpq == 2){
        my $r=$jbk->{$id};
        push(@{$both_bp_match},$r);
     }
   }
   #print Dumper($both_bp_match);
   return $both_bp_match;
}

#filter by type
sub _filter_by_svtype{
      my ($type,$rbk1)=@_;
      #we filter the list  of breapoints by chr
      my $f1=[];#array of filtered breakpoints1
      foreach my $r(@{$rbk1}){
        #print join(" ",$r->{type},$type)."\n";
          if($r->{type} eq $type){
               push(@{$f1},$r);
            }
        }
        return $f1;
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

  my ($f1)=_filter_by_chr($chr,$rbk1);
  #we change chr in case of a BND or a TRAN for breakpoint 2
  if($type eq "BND" or $type eq "TRA"){
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
                   id=>$item->{ID},len=>$item->{info}->{SVLEN}},$brk1->{start},$brk1->{stop});

    #anything that is not a Translocation
    if($type ne "BND" and $type ne "TRA"){
          #we add the second breakpoint
          $tree->insert({chr=>$item->{CHROM},brk=>2,index=>$i,type=>$item->{info}->{SVTYPE},
                        id=>$item->{ID},len=>$item->{info}->{SVLEN}},$brk2->{start},$brk2->{stop});
    }else{
        #print join(" ",$item->{info}->{CHR2},$item->{info}->{SVTYPE},$brk2->{start},$brk2->{stop},$item->{info}->{POS2})."\n";
          $tree->insert({chr=>$item->{info}->{CHR2},brk=>2,index=>$i,type=>$item->{info}->{SVTYPE},
                          id=>$item->{ID},len=>$item->{info}->{SVLEN}},$brk2->{start},$brk2->{stop});
    }
    $i++;#is an index for fast access to the structure
  }
  return $tree;
}

1;
