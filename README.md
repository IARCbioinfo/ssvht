# SSVht
A set of scripts to assist the calling of somatic structural variants from short reads using a random forest classifier.

## STEP 0: Collect Reference Databases

### Somatic SV reference (PCAWG)

Download consensus PCAWG calls from:

```
https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_sv/final_consensus_sv_bedpe_passonly.icgc.public.tgz
https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_sv/final_consensus_sv_bedpe_passonly.tcga.public.tgz
```

The PCAWG consensus call set for SVs (v1.6, n=309,246) is used as a reference for somatic SVs. The PCAWG SVs are in hg19 genome coordinates, thus is necessary to perform a liftover to GRCh38 using [crossmap](https://crossmap.readthedocs.io/en/latest/) (v0.3.9). PCAWG SVs at sample level (n=2,748) were merged into a non-redundant cohort callset with [SURVIVOR (v1.0.7)](https://github.com/fritzsedlazeck/SURVIVOR) (merge subcommand), leading to a total of 283,980 non-redundant somatic SVs.

### Germline SV reference (GNOMAD)

The gnomAD database (v2.1, n=299,211) GRCh38 liftover.
```
https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd166.GRCh38.variant_call.vcf.gz
``` 

### Genomic regions 
The following databases, BED format, should be provided,

1. [Cosmic v90](https://cosmic-blog.sanger.ac.uk/cosmic-release-v90/)
2. [gencode v33](https://www.gencodegenes.org/human/release_33.html)
3. [conserved genomic regions (100-way PhastCons)](https://gist.github.com/darencard/21860562a6edbc9fa12180f9df00381b) 
4. [centromeric genomics regions](http://genome.ucsc.edu/cgi-bin/hgTables) 

are used to annotate the location/impact of SVs. The directory ***bedfiles*** provide examples for each of the aforementioned BED files. Note that files should be uncompressed.


## STEP 1: Building a Custom Panel-of-Normal (PON)

Call germline SVs using [SvABA (v1.1.0)](https://github.com/walaj/svaba), [Manta (v1.6.0)](https://github.com/Illumina/manta), and [Delly (v0.8.3)](https://github.com/dellytools/delly) for building a custom panel of normal for each caller, follow the best practices for calling germline variants provided by each tool. Finally, make a consensus for each caller by integrating several samples with SURVIVOR (merge subcommand).
 

## STEP 2: Build random forest matrices

The matrices for training and predictions are built using the script ***svtools.pl***, the script should be called for each tumor and matched sample:

```
# manta example, similar calls should be done for DElly and SVaba

perl ssvht/svtools.pl -a Sample.tumor_only.pass.vcf \ 
						  -b Manta.PoN.vcf -c gnomad.2.1.GRCh38.vcf \ 
						  -d all_pcawg.hg38.non_redundant.vcf -s Sample.matched.pass.vcf \ 
						  -p <prefix> -x ssvht/bedfiles/bedfiles_list.txt 
```
In brief, the -a option provides the VCF calls that PASS filter for the tumor-only sample, -b provides the Panel-of-Normal file (VCF), -c and -d provide the reference SV databases, -s the somatic calls in case of training (or empty file), -p is the prefix and -x the path to bedfiles needed to annotate the location/impact of SVs.

The output is a matrix at the sample level with all the features needed by the model (see [BioRXiv](https://www.biorxiv.org/content/10.1101/2022.07.06.499003v1)), including:

```
Tumor ID CHROM POS TYPE SVLEN PE_SR PE SR PON PON_SUPP PON_TYPE PON_IDS PON_BC1 PON_BC2 GNOMAD GNOMAD_AC GNOMAD_TYPE GNOMAD_IDS GNOMAD_BC1 GNOMAD_BC2 PCAWG PCAWG_SUP PCAWG_TYPE PCAWG_IDS PCAWG_BC1 PCAWG_BC2 SOMATIC SOMATIC_TYPE SOMATIC_IDS COSMIC_GENE CENTROMER EXON CONSERVATION_BC1 CONSERVATION_BC2 CNV_TCN1 CNV_TCN2 CNV_CF1 CNV_CF2 RAF RFS
SAMPLE_T INV00000004 chr1 10433 INV 248935200 3 0 3 1 5 INV INV00000001 5 11 0 0 0 0 0 0 0 0 0 0 3 2 0 0 0 0 0 0 0 0 2 2 0 1 0.122222222222222 90
SAMPLE_T DUP00000008 chr1 789482 DUP 223225120 43 43 0 1 45 DUP DUP00000008 6 12 0 0 0 0 5 0 0 0 0 0 1 6 0 0 0 0 0 0 6 11 2 2 1 1 0.686567164179104 67
SAMPLE_T DEL00000009 chr1 789502 DEL 223222907 17 17 0 1 28 DEL DEL00000006 6 12 0 0 0 0 5 2 0 0 0 0 1 6 0 0 0 0 0 0 6 11 2 2 1 1 0.390243902439024 41
SAMPLE_T DEL00000018 chr1 934105 DEL 740 44 44 0 1 40 DEL DEL00000014 2 2 1 4 DUP nssv15963806 9 9 1 4 DEL SV0196854 9 9 0 0 0 0 0 0 6 6 2 2 1 1 0.978723404255319 47
SAMPLE_T DEL00000022 chr1 1028337 DEL 737 3 3 0 1 22 DEL DEL00000019 5 5 1 2134 DEL nssv15852976 4 6 0 0 0 0 2 2 0 0 0 0 0 0 6 6 2 2 1 1 0.0833333333333333 36
SAMPLE_T INV00000026 chr1 1067280 INV 1462 7 0 7 1 9 INV INV00000022 3 3 0 0 0 0 2 2 0 0 0 0 2 2 0 0 0 0 0 0 8 8 2 2 1 1 0.448979591836735 98
SAMPLE_T DEL00000027 chr1 1288706 DEL 1162 6 6 0 1 6 DEL DEL00000019 7 6 1 6 DUP nssv15960599 19 16 0 0 0 0 5 5 0 0 0 0 0 0 5 5 2 2 1 1 1 4
SAMPLE_T DEL00000028 chr1 1350077 DEL 1338 7 7 0 1 17 DEL DEL00000025 4 4 1 3456 DEL nssv15848677 5 5 0 0 0 0 2 2 0 0 0 0 0 0 5 5 2 2 1 1 0.538461538461538 13
SAMPLE_T DEL00000031 chr1 1454714 DEL 30726 3 3 0 0 0 0 0 0 1 0 0 0 0 9 12 0 0 0 0 4 3 0 0 0 0 0 1 5 5 2 2 1 1 0.0461538461538462 65
```

The same script can be used to build a matrix containing the set of predicted somatic SVs from matched genomic data (-s option). 

## STEP 3: Train and Predict Somatic SVs with the random forest models.

The directory ***rcode*** provide the RScripts for training, evaluation and classification of the random forest models for each caller. The input are the matrices generated on STEP 2. 


## STEP 4: Build a consensus with the predicted somatic calls 
	
After classifying the tumor-only SVs into somatic and germline for each of the SV callers, a consensus per sample can be built using the SURVIVOR code (merge subcommand).


# Final Note

Some of the scripts migth contain hard-coded paths which need to be replacing to match your data.


