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
   print "$0 usage : -a <vcf_tumor>  -b <pon_vcf>\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:", \%opts );
if ( !defined $opts{a}  ) {
   usage;
}

my $target = new SV($opts{a});
my $pon = new SV($opts{b});
#filter
my $ftype=0;#true means that vars are filers by type
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
#annotate using PANEL of normals
$pon->norm_svs(0);#do not load genotype information
$sva->annot_pon_sv($pon,$target,$ftype,$fdelta); #match target using the PON
#annotate GNOMAD
#$sv->annotate_gnomad();
#annoted COSMIC
#$sv->annotate_cosmic();
