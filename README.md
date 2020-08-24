# Pipeline for admixture analysis of South Asian populations
This is the pipeline I used for analysing the population structure and admixture patterns of South Asian populations for my MSc Thesis. 

## Data Acquisition
I used SNP data in VCF format from the 1000Genomes project. I downloaded the Chromosome 22 file containing SNP calls for every individual in the 1000Genomes project data set, to be trimmed to the populations I would be analysing.
I downloaded the gzipped **ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz** file from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ .

## Processing VCF file for my research purposes
The raw VCF file contains SNP data for every population in the 1000Genomes project, I am only interested in South Asian populations and a Han Chinese outgroup population.
In order to trim the VCF file accordingly, I used VCFTOOLS ( installation details at https://vcftools.github.io/index.html).

First, I generated a text file containing the ID numbers for every individual I wanted to retain for my analyses. I used the spreadsheet provided by 1000Genomes (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx) and kept individuals from the five South Asian populations in the data set (BEB, GIH, ITU, PJL, STU) and Southern Han Chinese (CHS).
I then recoded a new VCF file containing only my selected sample IDs:

```
vcftools --gzvcf file.vcf.gz --keep listnames.txt --recode --out new_file_vcf
```

I then used site filtering options in VCFTOOLS to filter genotypes with coverage below 50% and SNPs with minor allele count below 3:

```
vcftools --vcf all_pops_22_vcf.recode.vcf --max-missing 0.5 --mac 3 --recode --recode-INFO-all --out all_pops_22_filtered
```

The VCF file is ready for analysis.

## Principle Components Analysis

The first analytical technique I used to explore the dataset for population structure was Principle Components Analysis (PCA). I carried out PCA in RStudio using the package SNPRelate. The R script I used for this analysis 
