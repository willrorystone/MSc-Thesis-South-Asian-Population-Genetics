# Pipeline for admixture analysis of South Asian populations
This is the pipeline for analysing the population structure and admixture patterns of South Asian populations from my MSc Thesis. 

## Software required
I carry out all analyses in MacOS Terminal and Rstudio. Here is a list of software programs I use in chronological order with links for downloads and documentation:
* **VCFTOOLS** (https://vcftools.github.io/index.html)
* **PLINK** (https://www.cog-genomics.org/plink2)
* **ADMIXTURE** (http://software.genetics.ucla.edu/admixture/download.html)
* **BCFTOOLS** (http://www.htslib.org/download/)
* **TreeMix** (https://bitbucket.org/nygcresearch/treemix/wiki/Home)
* **Dsuite** (https://github.com/millanek/Dsuite)

## Data Acquisition
I use SNP data in VCF format from the 1000Genomes project. I downloaded the Chromosome 22 file containing SNP calls for every individual in the 1000Genomes project data set, to be trimmed to the populations I would be analysing.
I downloaded the gzipped **ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz** file from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ .

## Processing VCF file for my research purposes
The raw VCF file contains SNP data for every population in the 1000Genomes project, I am only interested in South Asian populations and a Han Chinese outgroup population.
In order to trim the VCF file accordingly, I use VCFTOOLS.

First, I generate a text file containing the ID numbers for every individual I want to retain for my analyses. I use the spreadsheet provided by 1000Genomes (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx) and keep individuals from the five South Asian populations in the data set (BEB, GIH, ITU, PJL, STU) and Southern Han Chinese (CHS).

I then recode a new VCF file containing only my selected sample IDs:

```
vcftools --gzvcf file.vcf.gz --keep listnames.txt --recode --out new_file_vcf
```

I then use site filtering options in VCFTOOLS to filter genotypes with coverage below 50% and SNPs with minor allele count below 3:

```
vcftools --vcf all_pops_22_vcf.recode.vcf --max-missing 0.5 --mac 3 --recode --recode-INFO-all --out all_pops_22_filtered
```

The VCF file is ready for analysis.


## Principle Components Analysis

The first analytical technique I use to explore the dataset for population structure is Principle Components Analysis (PCA). I conduct PCA in RStudio using the package SNPRelate. The R code for this analysis is in the R file PCA_script.R.


## ADMIXTURE Analysis

I use ADMIXTURE to infer ancestry components in the populations.

Before running ADMIXTURE software, I convert the VCF file into PLINK format for linkage disequilibrium (LD) pruning:

```
vcftools --vcf all_imputed_filtered.recode.vcf --plink --out all_pops_22_filtered
```

The following code conducts LD pruning on the PLINK snp file using the ```--indep-pairwise``` option.

```
/Users/willstone/Downloads/plink_mac_20200428/plink --file all_pops_22_filtered --indep-pairwise 50 5 0.5
```

The three numbers following ```--indep-pairwise``` stipulate that,  within an SNP window size of 50, one of any pair of SNPS with LD above 0.5 must be removed and then the window shifted 5 SNPs and the process repeated. The list of SNPs to be kept is contained within a file called plink.prune.in, and a new set PLINK files that contain these SNPs in .bed, .bim and .fam format are produced by:

```
/Users/willstone/Downloads/plink_mac_20200428/plink --file indian_22_filtered --extract plink.prune.in --make-bed --out all_pops_22_pruned
```

Now that the .bed, .bim and .fam files with LD pruned SNPs for my set of sample IDs are ready, I can run the ADMIXTURE software:

```
for K in 2 3 4 5 6

do ./admixture --cv all_pops_22_pruned.bed $K

done
```

This loop runs ADMIXTURE using user-defined values for K (representing the number of putative ancestral populations) from 2 to 6. The ```--cv``` option will output the cross validation error for each run into standard output.
