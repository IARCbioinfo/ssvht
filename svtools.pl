###############################################################################
# Author: Alex Di Genova
# Laboratory: ERABLE/INRIA
# Copyright (c)
# year: 2020
###############################################################################
use Data::Dumper;
use Getopt::Std;
use FindBin;
use lib "$FindBin::Bin";
use SV;
use SVannot;

use strict;

sub usage {
   print "$0 usage : -a <vcf_tumor>  -b <pon_vcf> -c <GNOMAD.vcf> -d <PCAWG.vcf> -e <CNV-READS> -s <Somatic.vcf> -p <prefix>\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:c:d:e:s:p:", \%opts );
if ( !defined $opts{a} or !defined $opts{p}) {
   usage;
}
#target file
my $target = new SV($opts{a});
#filter
my $ftype=1;#true means that vars are filers by type
my $fdelta=1000; #average distance for breakpoint overlap
#object to annotate SVs
my $sva=new SVannot();
#remove SVs on non-chr
#add type and lenght as well as method
#add SR and PE support to predictions
#store all in an internal array
#load the genotype information
$target->norm_svs(1);#load genotype information
#remove SVs shorter than 30 bp, matching to alternative chromosomes or with read support lower than 5
$target->basic_filters(3,30,500000000);

#load somatic variants of the sample
my $som = new SV($opts{s});
#annotate using PANEL of normals
$som->norm_svs(0);#do not load genotype information
#annotate custom PON SVs and the context
$sva->annot_Somatic_sv($som,$target,$ftype,$fdelta); #match target using the PON
$som=();#we free the memory of the variable

#load custom PON file
my $pon = new SV($opts{b});
#annotate using PANEL of normals
$pon->norm_svs(0);#do not load genotype information
#annotate custom PON SVs and the context
$sva->annot_customPON_sv($pon,$target,$ftype,$fdelta); #match target using the PON
$pon=();#we free the memory of the variable

#load GNOMAD PON file
my $gnomad= new SV($opts{c});
#do not load genotype information
$gnomad->norm_svs(0);
#annotate GNOMAD SVs and the context
$sva->annot_gnomad_sv($gnomad,$target,$ftype,$fdelta);
$gnomad=(); #we free the memory of the variable

#load PCAWG file
my $pcawg= new SV($opts{d});
$pcawg->norm_svs(0);
#annotate PCAWG SVs and the context
$sva->annot_pcawg_sv($pcawg,$target,$ftype,$fdelta);
$pcawg=();#we free the memory of the variable
#$sva->annot_breapoint_coverage();
#print the matrix to train the RF tool
$target->print_matrix($opts{p});
#annoted COSMIC SVs are of low quality and were replaced by PCAWG
#$sv->annotate_cosmic();
