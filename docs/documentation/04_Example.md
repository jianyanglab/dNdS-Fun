---
layout: page
title: Example 
description: ~
---
This tutorial provides an illustrative analysis of somatic mutation data using simulation data available here [link](https://yanglab.westlake.edu.cn/data/dNdSFun/example.txt) with dNdS-Fun. Before running the analysis, please ensure that dNdS-Fun is installed. For installation instructions, refer to this [link](https://jianyanglab.github.io/dNdS-Fun/documentation/02_installation.html).


## Input Data
`dNdS-Fun` supports two types of input data:
- Somatic mutation data with functional impact scores (annotated in the sixth column).
- Somatic mutation data without functional impact scores. This type of data requires an additional CADD score database, available here [link](https://jianyanglab.github.io/dNdS-Fun/documentation/03_data.html). Place the database in the path specified by our software or provide the correct path for the software to locate it. We support both GRCh37 and GRCh38 human genome; please use the corresponding model to ensure compatibility.

For instructions on annotating the functional impact score in the input data, please refer to the `Preparing Reference Data` section of this tutorial.

The example data required for this tutorial can be downloaded from the following [link](https://yanglab.westlake.edu.cn/data/dNdS-Fun/examples.tar.gz). Detailed information about the input data is provided below.


### 1. Input data format without functional impact scores
```r
# Load the example cds-exon data
input_dat = readRDS("../example_cds-exon.rds")

head(vcf_cds)

Sample1	chr14	20844257	G	T
Sample1	chr19	21558661	A	G
Sample2	chr1	192335033	C	A
Sample1	chr1	241875099	T	C
Sample2	chr16	48265863	G	A
```
A tab-delimited text file containing all somatic mutations sequenced from one or more individuals. Each row corresponds to a specific variant, and the five columns are: "IndividualID," "Chromosome," "Position," "Ref," and "Alt."

A table-delimited text file contained all somatic mutation sequenced from multiple individuals. Each row corresponds to a specific variant, and the 5 columns are separate "IndividualID, Chromosome, Position, Ref, Alt".

Note: Our software automatically adds CADD scores, but you can also provide mutation files with CADD scores. See below for the detailed format. If your input data format without functional impact scores, you will need to provide the CADD scores database, which can be downloaded from here [(GRCh37)](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz) or here [(GRCh38)](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz). You can use the parameter [`score_database`](https://jianyanglab.github.io/dNdS-Fun/documentation/03_data.html) to provide your CADD scores database path.


### 2. Split the CADD scores into different chromosomes

To speed up your analysis and reduce the resource footprint of your tasks, we offer multi-threaded processing. You can split the downloaded CADD scores into 22 chromosomes using the following code, allowing you to analyze them simultaneously with multiple threads.

```r

chr=22 #range[1-22]
wkDir="/your_path/"
mkdir=${wkDir}/chr/
wholefil="whole_genome_SNVs.tsv.gz" # downloaded from CADD scores, put into the current directory

rm ${wkDir}/chr/whole_genome_SNVs.tsv.gz.chr${chr}.gz.rankRawScore.GRCh37
echo -ne "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n" >${wkDir}/chr/whole_genome_SNVs.tsv.gz.chr${chr}.gz.rankRawScore.GRCh37
zcat ${wkDir}/${wholefil} | awk -v var=$chr '{if($1==var)print $0}' | sort | uniq >>${wkDir}/chr/whole_genome_SNVs.tsv.gz.chr${chr}.gz.rankRawScore.GRCh37
bgzip -c ${wkDir}/chr/whole_genome_SNVs.tsv.gz.chr${chr}.gz.rankRawScore.GRCh37 > ${wkDir}/chr/whole_genome_SNVs.tsv.gz.chr${chr}.gz.rankRawScore.GRCh37.gz
tabix -fp vcf ${wkDir}/chr/whole_genome_SNVs.tsv.gz.chr${chr}.gz.rankRawScore.GRCh37.gz
```
The path `${wkDir}/chr/` can be provided to dNdS-Fun as the parameter `score_database`. The maximum memory required for this step may reach up to 35GB.


### 3. Input data format with functional impact scores
```r
# Load the example cds-exon data
input_dat = readRDS("../example_cds-exon.rds")

head(impScore_cds)

Sample1	chr14	20844257	G	T	0.957112
Sample1	chr19	21558661	A	G	0.544228
Sample2	chr1	192335033	C	A	0.312798
Sample1	chr1	241875099	T	C	0.221462
Sample2	chr16	48265863	G	A	0.973034
```
The sixth column is the raw CADD score, annotated from the above the CADD database.


## Basic usage
```r
library(dNdSFun)

# Documentations
help(dNdSFun)
``` 
### 1. Global selection estimation of dNdS-Fun
This section provides an basic usage of dNdS-Fun.
- mutsFile:  
- reg: 
- dichotomylog: 
- globaldnds_outFile: 
- genelevel_outFile: 
- score:
- score_database: 
- iscv: NULL (The default).
- model: 1 or other numberic variables.

For further details about the parameters, please refer to this [link](https://jianyanglab.github.io/dNdS-Fun/documentation/03_data.html).
```r
#  
module load R/4.0.5
library(dNdSFun)

export R_LIBS_USER="/your_path/R_LIB_4.0.5"
Rscript src/cal_CADD_dndsWGS.NEG.R ${mutsFile} data/GRCh37/${reg}.GRCh37.rda data/GRCh37/${dichotomylog} ${reg} ${genelevel_outFile} ${iscv} ${model} ${score} ${score_database}
#or
Rscript src/cal_CADD_dndsWGS.NEG.R ${mutsFile} data/GRCh37/${reg}.GRCh37.rda data/GRCh38/${dichotomylog} ${reg} ${genelevel_outFile} ${iscv} ${model} ${score} ${score_database}
```
If you provide the somatic mutations with functional impact scores, you can ignore the parameter `score` and `score_database`.


### 2. Global selection estimation with dNdS-Fun
This section provides a global (genome-wide) selection estimation of somatic mutations using dNdS-Fun.
```r
# Set the gene selection method in dNdS-Fun function 
library(dNdSFun)
readRDS("../outFile_cds-exon.rds")

head(globaldnds_outFile)

region wsel ci_low ci_high pval selection
cds 1.29 1.23 1.37 0.001 positive
```
The results are stored in CADD_dndsWGSout.
- Global selection estimate: CADD_dndsWGSout$globaldnds_outFile
- wsel: ω ratio
- ci_low: lower 95% confidence of the ω ratio
- ci_high: lower 95% confidence of the ω ratio
- pval: p-value from the Chi-square test to determine whether the ω ratio differs from 1 (indicating no selection).
- selection: predicted selection direction

### 3. Local selection estimation with dNdS-Fun
The dNdS-Fun allows users to estimate selection at the local (gene or element-specific) scales. For five types of genomic regions, dNdS-Fun performs separate analyses for each gene. Here, we demonstrate one type of results: coding sequence (cds).

```r
# Set the gene selection method in dNdS-Fun function 
library(dNdSFun)
readRDS("../outFile_cds-exon.rds")

head(genelevel_outFile)

gene_name	n_neutral    n_select    wsel_cv   psel_cv qsel_cv
ENSG00000115524.11::SF3B1   0   15  148.10    9.06e-12    1.71e-07
ENSG00000160201.7::U2AF1    0   7   194.33    9.73e-09    0.00018
ENSG00000171862.5::PTEN     0   6   201.18    4.42e-08    0.00083
ENSG00000163554.7::SPTA1    2   11  18.27     9.67e-07    0.01834
```
The results are stored in CADD_dndsWGSout.
- Local selection estimate: CADD_dndsWGSout$genelevel_outFile.
- gene_name: ENSEMBL::GeneID.
- n_neutral: Number of somatic mutations in the less functional group
- n_select: Number of somatic mutations in the more functional group
- wsel: ω ratio for a specific gene
- pval: Statistical significance for a specific gene
- qval: Adjusted significance for a specific gene, using the Bonferroni multiple correction test

The dndsout data may contain genes without any mutations, showing a p-value of "0," or genes with only a few somatic mutations. We recommend selecting genes with more than five somatic mutations (n_neutral + n_select > 5) and an adjusted significance (qval < 0.05).

