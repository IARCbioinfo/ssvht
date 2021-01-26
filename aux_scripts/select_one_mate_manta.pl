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
   print "$0 usage : -a <manta>\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:", \%opts );
if ( !defined $opts{a}){
   usage;
}


#first pass to check which are the events
open(VCF, $opts{a}) or die "cannot open VCF file\n";
my $keep=();
while(my $line=<VCF>){
  chomp $line;
  next if($line=~m/^#/);
  my @d=split("\t",$line);
  my $tags=parse_tags($d[7]);
  #print Dumper($tags);
  #simple entry
  if(!defined $tags->{MATEID}){
      $keep->{$d[2]}=1;
  }
  #we have to enter just one side of the mate, by defaul is the first one
  if(defined $tags->{MATEID}){
        if(!defined $keep->{$d[2]} and !defined $keep->{$tags->{MATEID}}){
          $keep->{$d[2]}=1;
        }
  }

}
close(VCF);

open(VCF, $opts{a}) or die "cannot open VCF file\n";

while(my $line=<VCF>){
  chomp $line;
  if($line=~m/^#/){
    print $line."\n";
  }else{
  my @d=split("\t",$line);
    if(defined $keep->{$d[2]}){
      print $line."\n";
    }
  }
}
close(VCF);
#print Dumper($keep);

#parse the genotype information of the tool
sub _get_geno{
  my ($line)=@_;
  #chomp $entry;
  my $geno=();
  my @d=split("\t",$line);
  my @gtags=split(":",$d[8]);
  for(my $i=9; $i<=$#d; $i++){
    my @d2=split(":",$d[$i]);
    my $tags=();
    for(my $j=0; $j<=$#gtags; $j++){
          my @gv=split(",",$d2[$j]);
          push(@{$tags->{$gtags[$j]}},$_) foreach (@gv);

    }
    push(@{$geno},$tags);
  }
  #we return the set of genotypes
  return $geno;
}

#this is the col7
sub parse_geno_id{
  my ($col)=@_;
  my @vals=split(":",$col);
  $vals[7]=~s/_/:/g;
  return $vals[7];
}

sub parse_tags{
    my ($col)=@_;
    my $tags=();
    my @vals=split(";",$col);
    foreach my $t(@vals){
      my ($tag,$v)=split("=",$t);
          $tags->{$tag}=$v;
    }
    return $tags;
}

#parse SVaba tags
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



sub print_header{
  my $ctg_h=<<EOF;
##INFO=<ID=CALLERS,Number=1,Type=String,Description="SV callers supporting the SV">
##INFO=<ID=PES,Number=1,Type=Integer,Description="Maximum number of reads supporting the SV (pair-end+split)">
##FORMAT=<ID=ID,Number=1,Type=String,Description="Variant ID from input.">
##FORMAT=<ID=PE_SR,Number=1,Type=Integer,Description="Number of reads supporting the SV (pair-end+split)">
##FORMAT=<ID=RF_SP,Number=1,Type=Float,Description="Random forest probability of SV class somatic">
##FORMAT=<ID=SVT,Number=1,Type=String,Description="Type of the SV">
##FORMAT=<ID=SVL,Number=1,Type=Integer,Description="Length of the SV">
##FORMAT=<ID=RAF,Number=1,Type=Float,Description="Allele Frequency of the SV">
##FORMAT=<ID=RFS,Number=1,Type=Integer,Description="Read coverage of SV">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
##contig=<ID=chrM,length=16569>
EOF

return $ctg_h;
}
