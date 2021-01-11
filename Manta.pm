package Manta;

=head1 NAME
Manta

=head1 DESCRIPTION

This object perform specific operations on the output of Manta caller

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
   #we control the MATEID

my $events=();
my $mate_id=();
my $index=0;
my $id2pos=();
   foreach my $item (@{$svobj->{entries}}){
   #}
   #print Dumper($item);
   #print ref($item->{info})."\n";
   $item->{info}->{SVMETHOD}="MANTA";

   ##FORMAT=<ID=PR,Number=.,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">
   ##FORMAT=<ID=SR,Number=.,Type=Integer,Description="Split reads for the ref and alt alleles in the order listed, for reads where P(allele|rea
   my $sv_sup=0;
   my $ref_sup=0;
   #we ask for PE support
   if(!defined $item->{geno}[0]->{PR}){
     $item->{info}->{PE}=0;
   }else{
     #we add all the PE support for REF and ALT
     #$item->{info}->{PE}=$item->{geno}[0]->{PR}[0]+$item->{geno}[0]->{PR}[1];
     $item->{info}->{PE}=$item->{geno}[0]->{PR}[1];
     $sv_sup+=$item->{geno}[0]->{PR}[1];
     $ref_sup+=$item->{geno}[0]->{PR}[0];
   }
   if(!defined $item->{geno}->[0]{SR}){
     $item->{info}->{SR}=0;
   }else{
     #$item->{info}->{SR}=$item->{geno}[0]->{SR}[0]+$item->{geno}[0]->{SR}[1];
    $item->{info}->{SR}=$item->{geno}[0]->{SR}[1];
     $sv_sup+=$item->{geno}[0]->{SR}[1];
     $ref_sup+=$item->{geno}[0]->{SR}[0];
   }
   #total support of the SV for ALT
   $item->{info}->{PE_SR}=$item->{info}->{PE}+$item->{info}->{SR};
   #we compute the lengh of the SV
   #insertions has this varible already
   #<ID=SVLEN,Number=1,Type=Integer,Description="Insertion length for SVTYPE=INS.">
   my $type=$item->{info}->{SVTYPE};
   if($type eq "DEL" or $type eq "DUP" or $type eq "INV"){
          my $l=abs($item->{POS}-$item->{info}->{END});
          $item->{info}->{SVLEN}=$l;
   }

   #my $sv_sup=$item->{geno}->[0]->{DV}[0]+$item->{geno}->[0]->{RV}[0];
   #my $ref_sup=$item->{geno}->[0]->{DR}[0]+$item->{geno}->[0]->{RR}[0];
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
     #DELLY SVs
     #chr3    61351885        BND00008083     A       A[chr2:138988501[       2146    PASS    PRECISE;SVTYPE=BND;SVMETHOD=EMBL.DELLYv0.8.3;END=61351886;CHR2=chr2;POS2=138988501;PE=18;MAPQ=60;CT=3to5;CIPOS=-4,4;CIEND=-4,4;SRMAPQ=60;INSLEN=0;HOMLEN=5;SR=18;SRQ=1;CONSENSUS=CAATGGTACATATTTGCATATCTAAACACAGAAAAGTTACCGTAAAAATATAGAAGATAAAAATGGTAGATCTGTACAGGGCACTTGCCATGATTGGAGCTTGCAGGTTGCTGTGAGCTCCCTTCTTTCTTCTTTCTTTTTTGTTTCTAGTTGCTCTCTGATTCTCCAGTCATCCAATTTTCTCAGTGTTCAGGATACAGTAGGATACAAG;CE=1.95155;RDRATIO=1;SOMATIC       GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV   0/1:-47.6375,0,-152.027:10000:PASS:35420:78095:42675:2:68:17:56:23      0/0:0,-11.7257,-121.486:117:PASS:21762:45101:23339:2:47:1:39:0
     #chr5    3244539 BND00014259     T       ]chr3:31255793]T        2340    PASS    PRECISE;SVTYPE=BND;SVMETHOD=EMBL.DELLYv0.8.3;END=3244540;CHR2=chr3;POS2=31255793;PE=22;MAPQ=60;CT=5to3;CIPOS=-1,1;CIEND=-1,1;SRMAPQ=60;INSLEN=0;HOMLEN=0;SR=17;SRQ=0.990698;CONSENSUS=ATTCCCAAAACTCCATGTTTTAGGGGCCCACATGGCCTGCATTGCTGGCAATTTTATCTTGTTCTTACAGAAAACTTATACAGCTCATGATGCAGTTTTGTGAGTTCCACTCACCAGGAAGGCAAAAGGGAGCTGATGTGCAGAGATCATATAATGAGAGCAGAAGCAAGAGAGAGGAGTAGGAGGTGCTCGGCTTCCTTCAACAACCAGCTCTG;CE=1.99197;RDRATIO=1;SOMATIC      GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV   0/1:-52.1764,0,-103.644:10000:PASS:47679:118177:70498:2:55:22:43:24     0/0:0,-9.01764,-91.6867:90:PASS:20802:41779:20977:2:37:0:30:0
     #chr6    125746853       BND00017813     C       C[chr1:228020469[       1140    PASS    IMPRECISE;SVTYPE=BND;SVMETHOD=EMBL.DELLYv0.8.3;END=125746854;CHR2=chr1;POS2=228020469;PE=19;MAPQ=60;CT=3to5;CIPOS=-486,486;CIEND=-486,486;RDRATIO=1;SOMATIC     GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV   0/1:-89.6166,0,-347.617:10000:PASS:43175:69429:26254:2:62:19:0:0        0/0:0,-11.4391,-228:114:PASS:23723:47468:23745:2:38:0:0:0
     #chr7    83460601        BND00020264     T       T[chr4:87922736[        125     PASS    IMPRECISE;SVTYPE=BND;SVMETHOD=EMBL.DELLYv0.8.3;END=83460602;CHR2=chr4;POS2=87922736;PE=5;MAPQ=26;CT=3to5;CIPOS=-522,522;CIEND=-522,522;RDRATIO=1;SOMATIC        GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV   0/1:-4.55395,0,-181.201:46:PASS:39968:78385:38417:2:35:10:0:0   0/0:0,-7.52528,-135.3:75:PASS:16451:33639:17188:2:25:0:0:0
     #chr8    108316177       BND00023173     T       T]chr3:14425702]        2040    PASS    PRECISE;SVTYPE=BND;SVMETHOD=EMBL.DELLYv0.8.3;END=108316178;CHR2=chr3;POS2=14425702;PE=23;MAPQ=60;CT=3to3;CIPOS=-2,2;CIEND=-2,2;SRMAPQ=60;INSLEN=0;HOMLEN=1;SR=11;SRQ=0.970149;CONSENSUS=CAGAGCCCACCAGAATTATTCAAACTGGCCAGTCTTAAGCTATTTACTCTGCCATACGTTGCCTTTTCCATGGAAACCCCAAATAAAGGCTTTGGGCTAGATTTTTCTATCTATATACCCACCCCTCACCTCCTTAGAGCTTCCTTTCACTCCTGCAAAGTGTGTCACCCACTGCAGCCCCGGGTGTAATCAGCTTAT;CE=1.95665;RDRATIO=1;SOMATIC     GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV   0/1:-52.9672,0,-72.1807:10000:PASS:43554:73141:29587:2:32:24:29:24      0/0:0,-5.40994,-55.1914:54:PASS:16650:36764:20114:2:21:0:18:0
     #chr8    142500250       BND00023623     G       G]chr2:32916440]        215     PASS    PRECISE;SVTYPE=BND;SVMETHOD=EMBL.DELLYv0.8.3;END=142500251;CHR2=chr2;POS2=32916440;PE=0;MAPQ=0;CT=3to3;CIPOS=-4,4;CIEND=-4,4;SRMAPQ=60;INSLEN=0;HOMLEN=5;SR=5;SRQ=0.986842;CONSENSUS=GTTTCGCCCCCGCCCCACACCGCTCCTCCACTTCCCCCCGCGCCGCTCCTCCACCTCCCCACCCGCCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CE=0.874688;RDRATIO=1;SOMATIC      GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV   0/1:-3.41923,0,-21.3519:34:PASS:40161:70963:30802:2:0:0:11:5    0/0:0,-2.02045,-15.5132:20:PASS:24013:42049:18036:2:0:0:7:0
     #MANTA translocations (reciprocal translocation)
     #chr1    152103932       MantaBND:7015:0:4:1:0:0:1       T       [chrX:89128055[T        .       PASS    SVTYPE=BND;MATEID=MantaBND:7015:0:4:1:0:0:0;IMPRECISE;CIPOS=-160,161;EVENT=MantaBND:7015:0:4:0:0:0:0;SOMATIC;SOMATICSCORE=35;JUNCTION_SOMATICSCORE=0;BND_DEPTH=41;MATE_BND_DEPTH=16     PR      31,0    52,3
     #chr1    152104248       MantaBND:7015:0:4:0:0:0:0       G       G]chrX:89128210]        .       PASS    SVTYPE=BND;MATEID=MantaBND:7015:0:4:0:0:0:1;IMPRECISE;CIPOS=-132,132;EVENT=MantaBND:7015:0:4:0:0:0:0;SOMATIC;SOMATICSCORE=35;JUNCTION_SOMATICSCORE=0;BND_DEPTH=49;MATE_BND_DEPTH=25     PR      28,0    71,2
     #chrX    89128055        MantaBND:7015:0:4:1:0:0:0       G       [chr1:152103932[G       .       PASS    SVTYPE=BND;MATEID=MantaBND:7015:0:4:1:0:0:1;IMPRECISE;CIPOS=-181,182;EVENT=MantaBND:7015:0:4:0:0:0:0;SOMATIC;SOMATICSCORE=35;JUNCTION_SOMATICSCORE=0;BND_DEPTH=16;MATE_BND_DEPTH=41     PR      31,0    52,3
     #chrX    89128210        MantaBND:7015:0:4:0:0:0:1       T       T]chr1:152104248]       .       PASS    SVTYPE=BND;MATEID=MantaBND:7015:0:4:0:0:0:0;IMPRECISE;CIPOS=-181,181;EVENT=MantaBND:7015:0:4:0:0:0:0;SOMATIC;SOMATICSCORE=35;JUNCTION_SOMATICSCORE=0;BND_DEPTH=25;MATE_BND_DEPTH=49     PR      28,0    71,2
     #Expected format
     #chr1    152103932       MantaBND:7015:0:4:1:0:0:1       T       [chrX:89128055[T        .       PASS    SVTYPE=TRA;IMPRECISE;CIPOS=-160,161;SVMETHOD=MANTA;END=61351886;CHR2=chrX;POS2=89128055;EVENT=MantaBND:7015:0:4:0:0:0:0;SOMATIC;SOMATICSCORE=35;JUNCTION_SOMATICSCORE=0;BND_DEPTH=41;MATE_BND_DEPTH=16;     PR      31,0    52,3
           $item->{info}->{SVLEN}=0;
           $item->{info}->{SVTYPE}="TRA";#change BND for translocation
           my ($chr_pos)=_parse_trans_name($item->{ALT});
           $item->{info}->{CHR2}=$chr_pos->{CHR2};
           $item->{info}->{POS2}=$chr_pos->{POS2};
   }

   #we translocation events
   if(defined $item->{info}->{EVENT} and $item->{info}->{SVTYPE} eq "TRA"){
     push(@{$events->{$item->{info}->{EVENT}}},$index);
   }elsif(defined $item->{info}->{MATEID}){
      #is an event involving just 2 ids
     $mate_id->{$item->{ID}}=$item->{info}->{MATEID};
     $id2pos->{$item->{ID}}=$index;
    }
   $index++;
 }

 #print Dumper($item);
 #we collapse multiples BND to a single line description, we pick the BND with the highest read support.
 #print Dumper($events);
 #print Dumper($mate_id);
 #print Dumper($id2pos);
 #print join("\t",$item->{ID},$item->{info}->{PE},$item->{info}->{SR},$type,$item->{info}->{SVLEN})."\n";
 #print Dumper($item);

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
