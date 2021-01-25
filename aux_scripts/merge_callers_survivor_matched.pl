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
my $s_svaba=load_caller_data($opts{c},"SVaba");
my $s_delly=load_caller_data($opts{b},"Delly");
my $s_manta=load_caller_data($opts{a},"Manta");



#print Dumper($s_manta);
system("echo \"$opts{a}\n$opts{b}\n$opts{c}\" > $opts{p}.lst") == 0 or die "cannot create file for survivor failed: $?\n";
system("$SURVIVOR merge $opts{p}.lst 1000 1 0 0 0 50 $opts{p}.survivor.vcf > $opts{p}.survivor.log") ==0 or die "cannot run surviror failed: $?\n";
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
     if($tags->{SUPP} > 1 ){
        #print join(" ",$tags->{SUPP},$d[2],$geno1,$geno2,$geno3,join(",",@caller),join(":",@rsupp),join(":",@rtype),join(":",@rlen),join(":",@rfp))."\n";
        #print join(" ",$tags->{SUPP},$d[2],join(",",@caller),"ID:PE_SR:RF_SP:SVTYPE:SVLEN:RAF:RFS",join(":",@f_manta),join(":",@f_delly),join(":",@f_svaba))."\n";
        # we use the SURVIVOR calls
        if($max_pe_sr >=5){
        $d[7]="CALLERS=".join(",",@caller).";PES=$max_pe_sr;".$d[7];
        print OUT join("\t",@d[0..7],"ID:PE_SR:RF_SP:SVT:SVL:RAF:RFS",join(":",@f_manta),join(":",@f_delly),join(":",@f_svaba))."\n";
       }


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
    my ($file,$caller)=@_;
    open(VCF,$file) or die "cannot open file\n";
    my $caller2id=();
    while(my $line=<VCF>){
     next if ($line=~m/^#/);
     chomp $line;
     my @d=split("\t",$line);
     my $tags=parse_tags($d[7]);
     $caller2id->{$d[2]}->{tags}=$tags;
     $caller2id->{$d[2]}->{line}=$line;
     #we complete the read support, AF and DP
     my $sv_sup=0;#read support for the call in the tumor genome
     my $ref_sup=0;#reference support for the call
     my $geno=_get_geno($line);
     #print Dumper($geno);
     if($caller eq "Manta"){
       #pair-end support
       if(defined $geno->[1]->{PR}){
         $ref_sup+=$geno->[1]->{PR}[0];
         $sv_sup+=$geno->[1]->{PR}[1];
         #split-read suport
       }
       if(defined $geno->[1]->{SR}){
         $ref_sup+=$geno->[1]->{SR}[0];
         $sv_sup+=$geno->[1]->{SR}[1];
       }

       #Manta defines the SVLEN
       $caller2id->{$d[2]}->{tags}->{PE_SR}=$sv_sup;
     }elsif($caller eq "Delly"){
        #GENO
        ##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference pairs">
        ##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant pairs">
        ##FORMAT=<ID=RR,Number=1,Type=Integer,Description="# high-quality reference junction reads">
        ##FORMAT=<ID=RV,Number=1,Type=Integer,Description="# high-quality variant junction reads">
        #for delly in matched data the tumor BAM is raw 0 to check
        $sv_sup=$geno->[0]->{DV}[0]+$geno->[0]->{RV}[0];
        $ref_sup=$geno->[0]->{DR}[0]+$geno->[0]->{RR}[0];
        #print join("\t","PE=".$tags->{PE},"SR=".$tags->{SR},$sv_sup,$ref_sup)."\n";
        $caller2id->{$d[2]}->{tags}->{PE_SR}=$tags->{PE}+$tags->{SR};
        #only defines the SVLEN for insertions
        my $type=$tags->{SVTYPE};

        if($type eq "DEL" or $type eq "DUP" or $type eq "INV"){
               my $l=abs($d[1]-$tags->{END});
                $caller2id->{$d[2]}->{tags}->{SVLEN}=$l;
        }

     }elsif($caller eq "SVaba"){
       #total support of the SV for ALT
       $sv_sup=$geno->[1]->{AD}[0];
       #we compute the lengh of the SV

       #we add total supporting reads for ref/alt
       $ref_sup=$geno->[1]->{DP}[0];
       $ref_sup=1 if($ref_sup ==0);
       $caller2id->{$d[2]}->{tags}->{PE_SR}=$sv_sup;
       $caller2id->{$d[2]}->{tags}->{RAF}=$sv_sup/$ref_sup;
       $caller2id->{$d[2]}->{tags}->{RFS}=$ref_sup;
      my $type=$tags->{SVTYPE};
      if($type eq "DEL" or $type eq "DUP" or $type eq "INV" or $type eq "INS"){
               my ($chr_pos)=_parse_trans_name($d[4]);
               $tags->{END}=$chr_pos->{POS2};
               $tags->{SVLEN}=abs($d[1]-$tags->{END});
       }
       #print Dumper($geno);


     }

     my $rfs=$sv_sup+$ref_sup;
     if($rfs == 0){
       $rfs=1;
     }
     #we add AF
     my $raf=($sv_sup)/($rfs);

     #print join("\t",$ref_sup,$sv_sup, "RAF=$raf","RFS=$rfs","Manta")."\n";
     $caller2id->{$d[2]}->{tags}->{RF_SP}=1;
     #RAF and RFS are computed when not defined
     $caller2id->{$d[2]}->{tags}->{RAF}=$raf if(!defined $caller2id->{$d[2]}->{tags}->{RAF});
     $caller2id->{$d[2]}->{tags}->{RFS}=$rfs if(!defined $caller2id->{$d[2]}->{tags}->{RFS});
     $caller2id->{$d[2]}->{tags}->{SVTYPE}="BND"  if(!defined $caller2id->{$d[2]}->{tags}->{SVTYPE});
     $caller2id->{$d[2]}->{tags}->{SVLEN}=0 if(!defined $caller2id->{$d[2]}->{tags}->{SVLEN});

     #print Dumper($caller2id->{$d[2]}->{tags});


   }
   close(VCF);
   return $caller2id;
}


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
