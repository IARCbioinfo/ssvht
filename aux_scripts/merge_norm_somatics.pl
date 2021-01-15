##############################################################################
# Author: Alex Di Genova
# Laboratory: IARC/SCG/RCG
# Copyright (c)
# year: 2021
###############################################################################
use Data::Dumper;
use Getopt::Std;
use strict;

sub usage {
   print "$0 usage : -a <normal> -b <somatics> -c <output>\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:c:", \%opts );
if ( !defined $opts{a}  or !defined $opts{b} ) {
   usage;
}


#first passs we associate sv clusters

open(VCFB, $opts{b}) or die "cannot open VCF file $opts{b}\n";

#we load the somatics mutations
my @somatics=();
my $s_bam="";
my $t_som=0;
while(my $line=<VCFB>){
  #next if($line=~m/^#/);
  chomp $line;
  if($line=~m/^#/){
    if($line=~m/^#CHROM/){
      my @tmp=split("\t",$line);
      $s_bam=$tmp[10];
    }
  }else{
  chomp $line;
  my @d=split("\t",$line);
  my $tmp=();
  $tmp->{CHR}=$d[0];
  $tmp->{POS}=$d[1];
  $d[2].=":SOM";
  $tmp->{line}=join("\t",@d[0 .. 8],$d[10]);
  push(@somatics,$tmp);
  $t_som++;
 }
}


my $cs=shift(@somatics);
my $t_som_p=0;

open(VCFA, $opts{a}) or die "cannot open VCF file $opts{a}\n";
while (my $line=<VCFA>) {
chomp $line;
if($line=~m/#/){
  if($line=~m/#CHROM/){
      my @h=split("\t",$line);
      print join("\t",@h[0 .. 8],$h[10])."\n";
      if($h[10] ne $s_bam){
        print STDERR "DIFF BAMS $h[10] $s_bam"."\n";
        exit 1;
      }
  }else{
    print $line."\n";
  }
}else{
  #print $line."\n";
  my @d=split("\t",$line);
  $d[2].=":GERM";
  while($d[0] eq $cs->{CHR} and $cs->{POS} <= $d[1]){
        print $cs->{line}."\n";
        $cs=shift(@somatics);
        $t_som_p++;
  }
  print join("\t",@d[0 .. 8],$d[10])."\n";
}

}


if($t_som ne $t_som_p){
  print STDERR "ERROR: not all somatics were printed [$t_som, $t_som_p]\n";
}else{
  print STDERR "All Somatics were print: Total sommatics $t_som and total printed $t_som_p\n";
}
