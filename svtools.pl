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


use strict;

sub usage {
   print "$0 usage : -a <vcf>  -b <method>\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:", \%opts );
if ( !defined $opts{a}  ) {
   usage;
}

my $sv = new SV($opts{a});
#remove SVs on non-chr
#add type and lenght as well as method
#add SR and PE support to predictions
#store all in an internal array
$sv->norm_svs();
