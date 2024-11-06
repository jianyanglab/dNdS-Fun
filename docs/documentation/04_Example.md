---
layout: page
title: Example (dNdS-Fun)
description: ~
---
This tutorial offers an illustrative analysis of the somatic mutation data from simulation data [link](https://yanglab.westlake.edu.cn/data/dNdS-Fun/examples.tar.gz) using dNdS-Fun. Prior to running the analysis, it is important to ensure that the dNdS-Fun has been installed. For installation instructions, please refer to the following [link](https://jianyanglab.github.io/dNdS-Fun/documentation/02_installation.html).


## Input Data
`dNdS-Fun` are available for two types of input data:
- Somatic mutation data with functional impact scores (Annotated at the sixth column).
- Somatic mutation data without functional impact scores. (This type data are required an additional CADD score database [link](https://jianyanglab.github.io/dNdS-Fun/documentation/03_data.html). Place it in the path specified by our software or provide the software with a correct path to find it. We support the GRCh37 and GRCh38 both, please use the corresponding model to match.)

For how to annotate the functional impact score into the input data, please read the section of `Preparing Reference Data` in this tutorial. 

The example data required for running this tutorial can be downloaded from the following [link](https://yanglab.westlake.edu.cn/data/dNdS-Fun/examples.tar.gz). 
Detailed information regarding the input data is provided as follows.

### 1. The format of input data without functional impact scores
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
A table-delimited text file contained all somatic mutation sequenced from one or multiple individuals. Each row corresponds to a specific variant, and the 5 columns are separate "IndividualID, Chromosome, Position, Ref, Alt".

Note: Our software provides automatic addition of CADD scores, or you can also provide mutation files with CADD scores. See below for a detailed format. If your mutation withou functional impact score, you need to provide the CADD scores database, collected from [link](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz)(GRCh37) or [link](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz)(GRCh38). You can use the parameter `score_database` [link](https://jianyanglab.github.io/dNdS-Fun/documentation/01_About.html) to provide your CADD scres database.


### 2. split the CADD scores into different chromosomes

To speed up your analysis and reduce the resource footprint of your tasks, we offer multi-threaded tasks. You can divide the downloaded CADD scores into 22 chromosomes according to the following code so that you can analyze them at the same time with multiple threads.

```r

chr=22 #range[1-22]
wkDir="/your_path/"
mkdir=${wkDir}/chr/
wholefil="whole_genome_SNVs.tsv.gz" # downloaded from CADD scores, put into the current directory

rm ${wkDir}/chr/whole_genome_SNVs.tsv.gz.chr${chr}.gz.rankRawScore.GRCh38
echo -ne "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n" >${wkDir}/chr/whole_genome_SNVs.tsv.gz.chr${chr}.gz.rankRawScore.GRCh38
zcat ${wkDir}/${wholefil} | awk -v var=$chr '{if($1==var)print $0}' | sort | uniq >>${wkDir}/chr/whole_genome_SNVs.tsv.gz.chr${chr}.gz.rankRawScore.GRCh38
bgzip -c ${wkDir}/chr/whole_genome_SNVs.tsv.gz.chr${chr}.gz.rankRawScore.GRCh38 > ${wkDir}/chr/whole_genome_SNVs.tsv.gz.chr${chr}.gz.rankRawScore.GRCh38.gz
tabix -fp vcf ${wkDir}/chr/whole_genome_SNVs.tsv.gz.chr${chr}.gz.rankRawScore.GRCh38.gz
```
The path `${wkDir}/chr/` could be provided to dNdS-Fun as parameter `score_database`. The largest memory of this step should reach to 35Gb.


### 3. The format of input data with functional impact scores
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
A table-delimited text file contained all somatic mutation sequenced from one or multiple individuals. Each row corresponds to a specific variant, and the five columns are separate "IndividualID, Chromosome, Position, Ref, Alt, Functional_impact_scores (raw score). In the dNdS-Fun with CADD model, the raw score should be the CADD scores download from the [link](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz)(GRCh37) or [link](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz)(GRCh38).



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

For further details about the parameters, please refer to this [link](https://jianyanglab.github.io/dNdS-Fun/documentation/01_About.html).
```r
#  
module load R/4.0.5
library(dNdSFun)

export R_LIBS_USER="/your_path/R_LIB_4.0.5"
Rscript src/cal_CADD_dndsWGS.NEG.R ${mutsFile} data/GRCh37/${reg}.GRCh37.rda data/GRCh37/${dichotomylog} ${reg} ${genelevel_outFile} ${iscv} ${model} ${score} ${score_database}
#or
Rscript src/cal_CADD_dndsWGS.NEG.R ${mutsFile} data/GRCh38/${reg}.GRCh38.rda data/GRCh38/${dichotomylog} ${reg} ${genelevel_outFile} ${iscv} ${model} ${score} ${score_database}
```
If you provide the somatic mutations with functional impact scores, you can ignore the parameter `score` and `score_database`.


### 2. Global selection estimation of dNdS-Fun
This section provides an global selection estimation of somatic mutation of dNdS-Fun.
```r
# Set the gene selection method in dNdS-Fun function 
library(dNdSFun)
readRDS("../outFile_cds-exon.rds")

head(globaldnds_outFile)

region wsel    ci_low  ci_high AIC wsel_qua    ci_low_qua  ci_high_qua iscv    pval
cds-exon 1.299900469 1.232365546 1.371136377 2202.259634 1.299900469 1.21040382  1.396014456 nocv
```
The results are stored in CADD_dndsWGSout.
- The estimated global selection: CADD_dndsWGSout$globaldnds_outFile
- wsel: dN/dS ratio
- ci_low: lower 95% confidence
- ci_high: lower 95% confidence
- wsel_qua: dN/dS ratio (function: quaPoission. It is used when the number of mutations is low, to avoid overdispersion)
- ci_low_qua: lower 95% confidence (function: quaPoission. It is used when the number of mutations is low, to avoid overdispersion.)
- ci_high_qua: lower 95% confidence (function: quaPoission. It is used when the number of mutations is low, to avoid overdispersion.)
- iscv: Declare whether covariates are used in the model.
- pval: The significance of Chiq-square test to estimate whether the wsel is different from 1 (neutral selection).


### 3. Gene/element-level selection estimation of dNdS-Fun
dNdS-Fun provides users to estimate the gene/element-level selection. As for five types of genomic regions, dNdS-Fun would do the separate analysis for each gene. Here, we demonstrate the two types of results, including cds-exon and promoter regions.

```r
# Set the gene selection method in dNdS-Fun function 
library(dNdSFun)
readRDS("../outFile_cds-exon.rds")

head(genelevel_outFile)

gene_name   n_neutral   n_select    wsel_cv psel_cv qsel_cv
ENSG00000141510.11::TP53    0   136 6060.85135038872    0   0
ENSG00000133703.7::KRAS 2   260 3458.86768343209    0   0
ENSG00000138413.9::IDH1 0   32  975.940758569441    0   0
ENSG00000168036.12::CTNNB1  0   21  319.417061732786    0   0
ENSG00000121879.3::PIK3CA   0   23  304.785430640954    0   0
ENSG00000115524.11::SF3B1   0   15  148.106979984788    9.06019703705851e-12    1.71699794049296e-07
ENSG00000160201.7::U2AF1    0   7   194.339763938046    9.73021996308887e-09    0.000184387668300534
ENSG00000171862.5::PTEN 0   6   201.183923868826    4.4216568628741e-08 0.000837859758946013
ENSG00000163554.7::SPTA1    2   11  18.2768279520214    9.6787494618944e-07 0.0183392944803975
```
The results are stored in CADD_dndsWGSout.
- The estimated genelevel selection: CADD_dndsWGSout$genelevel_outFile.
- gene_name: ENSEMBL::GeneID.
- n_neutral: The number of somatic mutations in the non-functional group.
- n_select: The number of somatic mutations in the functional group.
- wsel: The dN/dS ratio of a specific gene.
- pval: The significance of a specific gene.
- qval: The adjust significance of a specific gene using Benjamini-Hochberg multiple correction test.

The dndsout data could contained the genes withou any mutations with a pval with "0" or contained genes with a few somatic mutations. We suggest that you should select the genes with somatic mutations more than 5 (n_neutral+n_select>5) and the adjust significance (qval<0.05).

Note: we do not provide a visulization function in dNdS-Fun.

