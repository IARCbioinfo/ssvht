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
   print "$0 usage : -a <manta> -b <delly> -c <svaba> -p <prefix>\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:c:p:", \%opts );
if ( !defined $opts{a}  or !defined $opts{b} or !defined $opts{c} or !defined $opts{p}){
   usage;
}

my $SURVIVOR="SURVIVOR";

#system(@args) == 0
#    or die "system @args failed: $?";
my $s_manta=load_caller_data($opts{a});
my $s_delly=load_caller_data($opts{b});
my $s_svaba=load_caller_data($opts{c});
#print Dumper($s_manta);
system("echo \"$opts{a}\n$opts{b}\n$opts{c}\" > $opts{p}.lst") == 0 or die "cannot create file for survivor failed: $?\n";
system("$SURVIVOR merge $opts{p}.lst 1000 1 0 0 0 0 $opts{p}.survivor.vcf > $opts{p}.survivor.log") ==0 or die "cannot run surviror failed: $?\n";
#merge matched.txt 1000 1 0 0 0 0 matched.vcf
#get_survivor_result("$opts{p}.survivor.vcf");
#ID:PE_SR:RF_SP:SVTYPE:SVLEN:RAF:RFS
open(FILE,"$opts{p}.survivor.vcf") or die "cannot open file\n";
open(OUT,">$opts{p}.integration.vcf") or die "cannot create $opts{p}.integration.vcf\n";
    while(my $line=<FILE>){
     #next if ($line=~m/^#/);
     chomp $line;
     if($line =~m/^#/){
       if($line =~m/^#CHROM/){
         print OUT print_header();
         print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMANTA\tDELLY\tSVABA\n";
         #print $line."\n";
       }elsif($line=~m/##contig/ or $line=~m/##FORMAT/){
         # we do not print the alt contigs, just main chromosomes...
       }else{
         print OUT $line."\n";
       }
     }else{
     chomp $line;
     my @d=split("\t",$line);
     my $tags=parse_tags($d[7]);
     my @caller=();
     my @f_manta=();
     my @f_delly=();
     my @f_svaba=();
     my $max_pe_sr=0;
     my $trsup=0;
     # manta
     my $geno1=parse_geno_id($d[9]);
     if($geno1 ne "NaN"){
     push(@caller,"Manta");
     my $g1=$geno1;
         $g1=~s/:/_/g;
     push(@f_manta,$g1,
     $s_manta->{$geno1}->{tags}->{PE_SR},
     $s_manta->{$geno1}->{tags}->{RF_SP},
     $s_manta->{$geno1}->{tags}->{SVTYPE},
     $s_manta->{$geno1}->{tags}->{SVLEN},
     $s_manta->{$geno1}->{tags}->{RAF},
     $s_manta->{$geno1}->{tags}->{RFS});
     $trsup+=$s_manta->{$geno1}->{tags}->{PE_SR};
    $max_pe_sr =  $s_manta->{$geno1}->{tags}->{PE_SR} if($s_manta->{$geno1}->{tags}->{PE_SR} > $max_pe_sr);
     }else{
       push(@f_manta,"NaN","NaN","NaN","NaN","NaN","NaN","NaN");
       #@f_manta
     }
     # delly
     my $geno2=parse_geno_id($d[10]);
     if($geno2 ne "NaN"){
       push(@caller,"Delly");
       #push(@rsupp,$s_delly->{$geno2}->{tags}->{PE_SR});
       #push(@rtype,$s_delly->{$geno2}->{tags}->{SVTYPE});
       #push(@rlen,$s_delly->{$geno2}->{tags}->{SVLEN});
       #push(@rfp,$s_delly->{$geno2}->{tags}->{RF_SP});
       my $g2=$geno2;
           $g2=~s/:/_/g;
       push(@f_delly,$g2,
       $s_delly->{$geno2}->{tags}->{PE_SR},
       $s_delly->{$geno2}->{tags}->{RF_SP},
       $s_delly->{$geno2}->{tags}->{SVTYPE},
       $s_delly->{$geno2}->{tags}->{SVLEN},
       $s_delly->{$geno2}->{tags}->{RAF},
       $s_delly->{$geno2}->{tags}->{RFS});

       $trsup+=$s_delly->{$geno2}->{tags}->{PE_SR};
       $max_pe_sr =  $s_delly->{$geno2}->{tags}->{PE_SR} if($s_delly->{$geno2}->{tags}->{PE_SR} > $max_pe_sr);
     }else{
       #push(@rsupp,0);
       push(@f_delly,"NaN","NaN","NaN","NaN","NaN","NaN","NaN");
     }
     # svaba
     my $geno3=parse_geno_id($d[11]);

     if($geno3 ne "NaN"){
       push(@caller,"SVaba");
       #push(@rsupp,$s_svaba->{$geno3}->{tags}->{PE_SR});
       #push(@rtype,$s_svaba->{$geno3}->{tags}->{SVTYPE});
       #push(@rlen,$s_svaba->{$geno3}->{tags}->{SVLEN});
       #push(@rfp,$s_svaba->{$geno3}->{tags}->{RF_SP});
       my $g3=$geno3;
           $g3=~s/:/_/g;
       push(@f_svaba,$g3,
       $s_svaba->{$geno3}->{tags}->{PE_SR},
       $s_svaba->{$geno3}->{tags}->{RF_SP},
       $s_svaba->{$geno3}->{tags}->{SVTYPE},
       $s_svaba->{$geno3}->{tags}->{SVLEN},
       $s_svaba->{$geno3}->{tags}->{RAF},
       $s_svaba->{$geno3}->{tags}->{RFS});
       $max_pe_sr = $s_svaba->{$geno3}->{tags}->{PE_SR} if($s_svaba->{$geno3}->{tags}->{PE_SR} > $max_pe_sr);

        $trsup+=$s_svaba->{$geno3}->{tags}->{PE_SR};
     }else{
       push(@f_svaba,"NaN","NaN","NaN","NaN","NaN","NaN","NaN");
     }
     if($tags->{SUPP} > 1){
        #print join(" ",$tags->{SUPP},$d[2],$geno1,$geno2,$geno3,join(",",@caller),join(":",@rsupp),join(":",@rtype),join(":",@rlen),join(":",@rfp))."\n";
        #print join(" ",$tags->{SUPP},$d[2],join(",",@caller),"ID:PE_SR:RF_SP:SVTYPE:SVLEN:RAF:RFS",join(":",@f_manta),join(":",@f_delly),join(":",@f_svaba))."\n";
        # we use the SURVIVOR calls
        $d[7]="CALLERS=".join(",",@caller).";PES=$max_pe_sr;".$d[7];
        print OUT join("\t",@d[0..7],"ID:PE_SR:RF_SP:SVT:SVL:RAF:RFS",join(":",@f_manta),join(":",@f_delly),join(":",@f_svaba))."\n";


      }elsif($trsup >=15){
	#if($geno3 eq "NaN"){
        $d[7]="CALLERS=".join(",",@caller).";PES=$max_pe_sr;".$d[7];
        print OUT join("\t",@d[0..7],"ID:PE_SR:RF_SP:SVT:SVL:RAF:RFS",join(":",@f_manta),join(":",@f_delly),join(":",@f_svaba))."\n";
	#}elsif($trsup >= 20){
        #$d[7]="CALLERS=".join(",",@caller).";PES=$max_pe_sr;".$d[7];
        #print OUT join("\t",@d[0..7],"ID:PE_SR:RF_SP:SVT:SVL:RAF:RFS",join(":",@f_manta),join(":",@f_delly),join(":",@f_svaba))."\n";
	#}

      }
     }
  }



#callers
sub load_caller_data{
    my ($file)=@_;
    open(VCF,$file) or die "cannot open file\n";
    my $caller2id=();
    while(my $line=<VCF>){
     next if ($line=~m/^#/);
     my @d=split("\t",$line);
     my $tags=parse_tags($d[7]);
     $caller2id->{$d[2]}->{tags}=$tags;
     $caller2id->{$d[2]}->{line}=$line;
   }
   close(VCF);
   return $caller2id;
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
