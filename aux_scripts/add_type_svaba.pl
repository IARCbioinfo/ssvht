
###############################################################################
# Author: Alex Di Genova
# Laboratory: IARC/SCG/RCG
# Copyright (c)
# year: 2021
###############################################################################
use Data::Dumper;
use Getopt::Std;
use strict;

sub usage {
   print "$0 usage : -a  -b  -c\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:c:", \%opts );
if ( !defined $opts{a}  ) {
   usage;
}

#first passs we associate sv clusters
open(VCF, $opts{a}) or die "cannot open VCF file\n";
my $sgrup=();
while(my $line=<VCF>){
	chomp $line;
	next if($line =~m/^#/);
	my @d=split("\t",$line);
	#my ($t)=parse_tags($d[7]);
  my ($id)=split(":",$d[2]);
  push(@{$sgrup->{$id}},$line);
	#print Dumper($t);
}
close(VCF);
my ($types)=add_types($sgrup);
#print Dumper($types);

## we perform the second pass to add the types

my $header=<<EOF;
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="PE confidence interval around END">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="PE confidence interval around POS">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for POS2 coordinate in case of an inter-chromosomal translocation">
##INFO=<ID=POS2,Number=1,Type=Integer,Description="Genomic position for CHR2 in case of an inter-chromosomal translocation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Insertion length for SVTYPE=INS.">
EOF




open(VCF, $opts{a}) or die "cannot open VCF file\n";
while(my $line=<VCF>){
	chomp $line;
  if($line =~m/^#/){
    if($line =~m/^#CHROM/){
      print $header;
      print $line."\n";
    }else{
      print $line."\n";
    }
  }else{
    #print $line."\n";
    my @d=split("\t",$line);
    if(defined $types->{$d[2]}){
          #print $line."\n";
          my $t=$types->{$d[2]}->{type};
          my $l=$types->{$d[2]}->{len};
          my $e=$types->{$d[2]}->{END};
          #my $t=$types->{$d[2]}->{type};
          $d[7]=~s/SVTYPE=BND//g;
          $d[7].=join(";","SVMETHOD=SVABA","SVTYPE=$t","SVLEN=$l","END=$e","CIPOS=-300,300;CIEND=-300,300");
          if($t eq "BND"){
            my $p2=$types->{$d[2]}->{POS2};
            my $c2=$types->{$d[2]}->{CHR2};
            $d[7].=";CHR2=$c2;POS2=$p2";
          }
          #print "".$d[7]."\n";
          print join("\t",@d)."\n";
    }
  }

}

#function that add types
sub add_types{
  my ($sgrup)=@_;
  my $hash=();
  my $ins=0;
  my $del=0;
  my $inv=0;
  my $tra=0;
  my $total=0;
  #print Dumper($sgrup);
  foreach my $s (keys %{$sgrup}){
      my $part=scalar(@{$sgrup->{$s}});
      #my $p1,$p2;
      if($part == 2){
         my $p1=parse_line(@{$sgrup->{$s}}[0]);
         my $p2=parse_line(@{$sgrup->{$s}}[1]);
         my $b1=parse_brackets($p1->{ALT});
         my $b2=parse_brackets($p2->{ALT});

         my $type="BND";
         my $svlen=0;
         # p1:(\w\].+\]) :1 p2:(\w\].+\]):1
         # p1:(\[.+\[\w) :2 p2:(\[.+\[\w):2
         if(($b1 == $b2) and ($b1==1 or $b1==2)  and ($p1->{CHROM} eq $p2->{CHROM})){
           $type="INV";
           $svlen=abs($p1->{POS}-$p2->{POS})+1;
           $inv++;
         #}elsif(($b1 == 3 and $b2 == 4) || ($b1 == 4 and $b2 == 3)){
         # p1:(\w\[.+\[) :3 p2:\].+\]\w:4
       }elsif(($b1 == 3 and $b2 == 4) and ($p1->{CHROM} eq $p2->{CHROM})){
           $type="DEL";
           $svlen=abs($p1->{POS}-$p2->{POS})+1;
           $del++;

         # p1:\w].+\] 1  p2: \[.+\[\w :2
       }elsif(($b1 == 1 and $b2 == 2) and ($p1->{CHROM} eq $p2->{CHROM})){
            $type="INS";
            $svlen=abs($p1->{POS}-$p2->{POS})+1;
            $ins++;
        }else{
           $type="BND";
           $svlen=0;
           $tra++;
        }
         #CHROM
         #if($p1->{CHROM} ne $p2->{CHROM}){
          # $type="TRA";
          # $svlen=0;
          # $tra++;
         #}
         #print join(" ",$p1->{CHROM},$p1->{POS},$p1->{REF},$p1->{ALT},$p2->{CHROM},$p2->{POS},$p2->{REF},$p2->{ALT}, parse_brackets($p1->{ALT}),parse_brackets($p2->{ALT}),$type,$svlen)."\n";
         #print join(" ",$p1->{CHROM},$p1->{POS},$p1->{REF},$p1->{ALT},
         #                $p2->{CHROM},$p2->{POS},$p2->{REF},$p2->{ALT},
         #                parse_brackets($p1->{ALT}),
         #                parse_brackets($p2->{ALT}),$svlen,$type)."\n";
         #tags to add to info of
         $hash->{$p1->{ID}}->{type}=$type;
         $hash->{$p1->{ID}}->{len}=$svlen;
         $hash->{$p1->{ID}}->{CHR2}=$p2->{CHROM};
         $hash->{$p1->{ID}}->{POS2}=$p2->{POS};
         $hash->{$p1->{ID}}->{END}=$p2->{POS};
         if($type eq "BND"){
           $hash->{$p1->{ID}}->{END}=$p1->{POS}+1;
         }
      }elsif($part == 1){
        #my $p1=parse_line(@$sgrup->{$s}[0]);
        print  STDERR "skyping $s, it is composed of only one mate\n";
      }else{
        print STDERR "skyping $s, it composed of more than one mate\n";
      }

      $total++;
  }

  print STDERR join(" ","TRA=$tra","INV=$inv","DEL=$del","ins=$ins", "total=$total")."\n";

  return $hash;

}

sub parse_brackets{
   my ($b)=@_;

   if($b=~m/\w\].+\]/){
     return 1;
   }
   if($b=~m/\[.+\[\w/){
     return 2;
   }

   if($b=~m/\w\[.+\[/){
     return 3;
   }
   if($b=~m/\].+\]\w/){
     return 4;
   }

   return 5;
}

#we parse the VCF entry
sub parse_line{
     my ($line)=@_;
     my $h=();
     #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  C000KGK_alt_HW5NGCCXY.DUAL292_BQSRecalibrated.bam
     #chr1    136610  49:1    C       ]chr1:137149]C  31      PASS    DISC_MAPQ=27;EVDNC=DSCRD;IMPRECISE;MAPQ=27;MATEID=49:2;MATENM=-1;NM=-1;NUMPARTS=0;SCTG

     my @d=split("\t",$line);
     $h->{CHROM}=$d[0];
     $h->{POS}=$d[1];
     $h->{ID}=$d[2];
     $h->{REF}=$d[3];
     $h->{ALT}=$d[4];
     $h->{FORMAT}=parse_tags($d[7]);
     return $h;
}

sub parse_tags{
    my ($col)=@_;
    my $tags=();
    my @vals=split(";",$col);
	#print Dumper(@vals);
    foreach my $t(@vals){
      my ($tag,$v)=split("=",$t);
          $tags->{$tag}=$v;
    }
    return $tags;
}
