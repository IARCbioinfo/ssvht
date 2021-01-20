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
   print "$0 usage : -a <VCF> -b <RF-result> -s <SAMPLE> -t <caller>\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:s:t:", \%opts );
if ( !defined $opts{a}  or !defined $opts{b} or !defined $opts{s} ) {
   usage;
}


my ($sv_pass)=load_rf_calls($opts{s},$opts{b});
#print Dumper($sv_pass);

#pick useful values for random forest
my @cmodel=("PON_BC1","PON_BC2","GNOMAD_AC",
           "GNOMAD_BC1","GNOMAD_BC2",
           "PCAWG_BC1","PCAWG_BC2","COSMIC_GENE",
           ,"EXON","CONSERVATION_BC1",
           "CONSERVATION_BC2","CNV_TCN1", "CNV_TCN2",
           "CNV_CF1", "CNV_CF2", "SVLEN","RAF","RFS","Predicted.prob.1","PE_SR");


open(VCF, $opts{a}) or die "cannot open VCF file\n";
while(my $line=<VCF>){
	chomp $line;
  if($line =~m/^#/){
    if($line =~m/^#CHROM/){
      print get_header($opts{t});
      print $line."\n";
    }elsif($line=~m/##contig/){
      # we do not print the alt contigs, just main chromosomes...
    }else{
      print $line."\n";
    }
  }else{
    #print $line."\n";
    my @d=split("\t",$line);
    #print Dumper(@d);
    next if (!defined $sv_pass->{$d[2]});
    #print $line."\n";
    my ($tags)=parse_tags($d[7]);
    #my $new_tags=();
    #my $ntags="";
    my @ntags=();
    foreach my $v (@cmodel){
      if(!defined $tags->{$v}){
        #$new_tags->{$v}=$sv_pass->{$d[2]}->{$v};
        #$ntags="$v=",$sv_pass->{$d[2]}->{$v}.";";
        push(@ntags,"$v=".$sv_pass->{$d[2]}->{$v});
      }
    }
    my $n_tag=join(";",@ntags,$d[7]);
    #print Dumper($tags);
    $n_tag=~s/Predicted.prob.1/RF_SP/;
    #print Dumper($new_tags);
   $d[7]=$n_tag;
   print join("\t",@d)."\n";
  }
}



#function that load the RF calls
sub load_rf_calls{
  my ($sample, $file)=@_;

  open(RF,$file) or die "cannot open file $file\n";
  my $header=<RF>;
  chomp $header;
  my @cols=split("\t",$header);
  #print Dumper(@cols);
  my $hash=();
  while(my $line =<RF>){
      chomp $line;
      my @d=split("\t",$line);
      #print Dumper(@d);
      #we skyp calls based on the given sample
      next if($d[1] ne $sample);
      #print $line."\n";
      for(my $i=1; $i<$#d; $i++){
           $hash->{$d[2]}->{$cols[$i-1]}=$d[$i];
      }
  }
  close(RF);

  #print Dumper($hash);
  return $hash;
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



sub get_header{
  my ($caller)=@_;

  #info from the RF classifier
my $ctg_h=<<EOF;
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

my $hdelly=<<EOF;
##INFO=<ID=RF_SP,Number=1,Type=Float,Description="Random forest probability of SV class somatic">
##INFO=<ID=RAF,Number=1,Type=Float,Description="Allele Frequency of the SV">
##INFO=<ID=RFS,Number=1,Type=Integer,Description="Read coverage of SV">
##INFO=<ID=CONSERVATION_BC1,Number=1,Type=Integer,Description="CONSERVATION breakpoint 1">
##INFO=<ID=CONSERVATION_BC2,Number=1,Type=Integer,Description="CONSERVATION breakpoint 2">
##INFO=<ID=GNOMAD_BC1,Number=1,Type=Integer,Description="Number of GNOMAD SVs around breakpoint 1">
##INFO=<ID=GNOMAD_BC2,Number=1,Type=Integer,Description="Number of  GNOMAD SVs around breakpoint 2">
##INFO=<ID=GNOMAD_AC,Number=1,Type=Integer,Description="GNOMAD minor allele count">
##INFO=<ID=PCAWG_BC1,Number=1,Type=Integer,Description="Number of  PCAWG SVs around breakpoint 1">
##INFO=<ID=PCAWG_BC2,Number=1,Type=Integer,Description="Number of  PCAWG SVs around breakpoint 2">
##INFO=<ID=PON_BC1,Number=1,Type=Integer,Description="Number of PON SVs around breakpoint 1">
##INFO=<ID=PON_BC2,Number=1,Type=Integer,Description="Number of PON SVs around breakpoint 2">
##INFO=<ID=CNV_TCN1,Number=1,Type=Integer,Description="Facets total Copy number on breakpoint 1">
##INFO=<ID=CNV_TCN2,Number=1,Type=Integer,Description="Facets total Copy number on breakpoint 2">
##INFO=<ID=CNV_CF1,Number=1,Type=Float,Description="Facets cellular fraction on breakpoint 1">
##INFO=<ID=CNV_CF2,Number=1,Type=Float,Description="Facets cellular fraction on breakpoint 2">
##INFO=<ID=EXON,Number=1,Type=Integer,Description="SV hit an EXON">
##INFO=<ID=COSMIC_GENE,Number=1,Type=Integer,Description="SV hit an COSMIC gene">
##INFO=<ID=PE_SR,Number=1,Type=Integer,Description="Numbe of reads supporting the SV (pair-end+split)">
EOF

if($caller =~m/delly/){
  return $ctg_h.$hdelly;
}

my $hmanta=<<EOF;
##INFO=<ID=RF_SP,Number=1,Type=Float,Description="Random forest probability of SV class somatic">
##INFO=<ID=RAF,Number=1,Type=Float,Description="Allele Frequency of the SV">
##INFO=<ID=RFS,Number=1,Type=Integer,Description="Read coverage of SV">
##INFO=<ID=CONSERVATION_BC1,Number=1,Type=Integer,Description="CONSERVATION breakpoint 1">
##INFO=<ID=CONSERVATION_BC2,Number=1,Type=Integer,Description="CONSERVATION breakpoint 2">
##INFO=<ID=GNOMAD_BC1,Number=1,Type=Integer,Description="Number of GNOMAD SVs around breakpoint 1">
##INFO=<ID=GNOMAD_BC2,Number=1,Type=Integer,Description="Number of  GNOMAD SVs around breakpoint 2">
##INFO=<ID=GNOMAD_AC,Number=1,Type=Integer,Description="GNOMAD minor allele count">
##INFO=<ID=PCAWG_BC1,Number=1,Type=Integer,Description="Number of  PCAWG SVs around breakpoint 1">
##INFO=<ID=PCAWG_BC2,Number=1,Type=Integer,Description="Number of  PCAWG SVs around breakpoint 2">
##INFO=<ID=PON_BC1,Number=1,Type=Integer,Description="Number of PON SVs around breakpoint 1">
##INFO=<ID=PON_BC2,Number=1,Type=Integer,Description="Number of PON SVs around breakpoint 2">
##INFO=<ID=CNV_TCN1,Number=1,Type=Integer,Description="Facets total Copy number on breakpoint 1">
##INFO=<ID=CNV_TCN2,Number=1,Type=Integer,Description="Facets total Copy number on breakpoint 2">
##INFO=<ID=CNV_CF1,Number=1,Type=Float,Description="Facets cellular fraction on breakpoint 1">
##INFO=<ID=CNV_CF2,Number=1,Type=Float,Description="Facets cellular fraction on breakpoint 2">
##INFO=<ID=EXON,Number=1,Type=Integer,Description="SV hit an EXON">
##INFO=<ID=COSMIC_GENE,Number=1,Type=Integer,Description="SV hit an COSMIC gene">
##INFO=<ID=PE_SR,Number=1,Type=Integer,Description="Numbe of reads supporting the SV (pair-end+split)">
EOF

#header for manta
if($caller =~m/manta/){
  return $ctg_h.$hmanta;
}

my $hsvaba=<<EOF;
##INFO=<ID=RF_SP,Number=1,Type=Float,Description="Random forest probability of SV class somatic">
##INFO=<ID=RAF,Number=1,Type=Float,Description="Allele Frequency of the SV">
##INFO=<ID=RFS,Number=1,Type=Integer,Description="Read coverage of SV">
##INFO=<ID=CONSERVATION_BC1,Number=1,Type=Integer,Description="CONSERVATION breakpoint 1">
##INFO=<ID=CONSERVATION_BC2,Number=1,Type=Integer,Description="CONSERVATION breakpoint 2">
##INFO=<ID=GNOMAD_BC1,Number=1,Type=Integer,Description="Number of GNOMAD SVs around breakpoint 1">
##INFO=<ID=GNOMAD_BC2,Number=1,Type=Integer,Description="Number of  GNOMAD SVs around breakpoint 2">
##INFO=<ID=GNOMAD_AC,Number=1,Type=Integer,Description="GNOMAD minor allele count">
##INFO=<ID=PCAWG_BC1,Number=1,Type=Integer,Description="Number of  PCAWG SVs around breakpoint 1">
##INFO=<ID=PCAWG_BC2,Number=1,Type=Integer,Description="Number of  PCAWG SVs around breakpoint 2">
##INFO=<ID=PON_BC1,Number=1,Type=Integer,Description="Number of PON SVs around breakpoint 1">
##INFO=<ID=PON_BC2,Number=1,Type=Integer,Description="Number of PON SVs around breakpoint 2">
##INFO=<ID=CNV_TCN1,Number=1,Type=Integer,Description="Facets total Copy number on breakpoint 1">
##INFO=<ID=CNV_TCN2,Number=1,Type=Integer,Description="Facets total Copy number on breakpoint 2">
##INFO=<ID=CNV_CF1,Number=1,Type=Float,Description="Facets cellular fraction on breakpoint 1">
##INFO=<ID=CNV_CF2,Number=1,Type=Float,Description="Facets cellular fraction on breakpoint 2">
##INFO=<ID=EXON,Number=1,Type=Integer,Description="SV hit an EXON">
##INFO=<ID=COSMIC_GENE,Number=1,Type=Integer,Description="SV hit an COSMIC gene">
##INFO=<ID=PE_SR,Number=1,Type=Integer,Description="Number of reads supporting the SV (pair-end+split)">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
EOF

#header for manta
if($caller =~m/svaba/){
  return $ctg_h.$hsvaba;
}


}
