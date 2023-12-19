---
layout: page
title: Example (hPSC)
description: ~
---

This tutorial provides an illustrative analysis of the hPSC dataset from [Chu et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1033-x) using MeDuSA. 


In this tutorial, we will use a dataset obtained from the hPSC cell line to estimate cell-state abundance along the hPSC differentiation trajectory in bulk RNA-seq data using MeDuSA. We will then validate the performance of MeDuSA by comparing the estimated cell-state abundance with that measured from scRNA-seq data. 

Prior to running the analysis, it is important to ensure that the MeDuSA package has been installed. For installation instructions, please refer to the following [link](https://leonsong1995.github.io/MeDuSA/documentation/02_installation.html).


## Input Data
`MeDuSA` requires two types of input data:
- Bulk RNA-seq data. 
- Single cell RNA-seq (scRNA-seq) data, which should be provided in the form of a Seurat object that includes the annotated cell-state trajectory and cell types. 

For how to prepare the cell-state trajectory data, please read the section of `Preparing Reference Data` in this tutorial. 

The input data required for running this tutorial can be downloaded from the following [link](https://yanglab.westlake.edu.cn/data/MeDuSA_data/hPSC.tar.gz). 
Detailed information regarding the input data is provided as follows.

### 1. Bulk RNA-seq Data
```r
# Load the example bulk RNA-seq data
bulk = readRDS("../hPSC_bulk.rds")
```
The bulk RNA-seq data is represented in a matrix format, where each row corresponds to a specific gene and each column corresponds to a particular sample.

### 2. Reference scRNA-seq Data
```r
# Load the example scRNA-seq data
sce = readRDS("./hPSC_sce.rds")
class(sce)
[1] "Seurat"
attr(,"package")
[1] "SeuratObject"

sce@assays$RNA@counts[1:3,1:3]
      H1_Exp1.001 H1_Exp1.002 H1_Exp1.003
MKL2  0.006680284 0.074591629 0.001734291
CD109 0.004262021 0.001206358 0.096426585
ABTB1 0.000000000 0.012892380 0.000000000       

sce$cell_trajectory[1:3]
H1_Exp1.001 H1_Exp1.002 H1_Exp1.003 
  0.8623402   0.7571288   0.6784661
		 
sce$cell_type[1:3]
H1_Exp1.001 H1_Exp1.002 H1_Exp1.003 
     "hPSC"      "hPSC"      "hPSC"
```
For compatibility with MeDuSA, the reference scRNA-seq data must be in the Seurat object format. Specifically, the reference data should be stored in `sce@assays$RNA@counts`, the cell-state trajectory in `sce$cell_trajectory`, and the cell-type in `sce$cell_type`. For more information about Seurat, please refer to the following [resource](https://satijalab.org/seurat/).


## Cell-State Deconvolution Analysis
```r
library(MeDuSA)

# Documentations
help(MeDuSA)
``` 
### 1. Basic Usage of MeDuSA
This section provides an introduction to the basic usage of MeDuSA.
- bulk: A matrix of bulk RNA-seq data. 
- sce: A Seurat object of scRNA-seq data.  
- select.ct: A character variable indicating the focal cell type.
- markerGene: A character vector containing the marker genes across the cell-state trajectory.If not provided, MeDuSA will utilize the `MeDuSA_marker` function to select marker genes for the analysis.
- resolution: A numeric variable used to specify the number of cell-state bins along the cell trajectory.
- smooth: A boolean variable to determine whether to smooth the estimated cell-state abundance.
- span: A numeric variable to control the degree of smoothing.
- fractional: A boolean variable to determine whether to normalize the estimated cell-state abundance to the fractional abundance (0-1).
- ncpu: The number of CPU cores to be used. 

For further details about the parameters, please refer to this [link](https://leonsong1995.github.io/MeDuSA/documentation/01_About.html).
```r
MeDuSA_obj = MeDuSA(bulk,sce,
                  select.ct = 'embry',markerGene = NULL,span = 0.35,
		  resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)		 
```
The results are stored in `MeDuSA_obj@Estimation`.
- The estimated cell-state abundance: `MeDuSA_obj@Estimation$cell_state_abundance`
- The median state (pseudo-time) of cell-state bins: `MeDuSA_obj@Estimation$TimeBin`
- The used marker genes: `MeDuSA_obj@Estimation$markerGene`

### 2. P-Values of the Random Effects Component
After completing the deconvolution analysis using MeDuSA, users can utilize the MeDuSA_VarExplain function to obtain the explained variance of the bulk data by the reference scRNA-seq data, as well as the corresponding p-values.
```R
MeDuSA_obj = MeDuSA_VarExplain(MeDuSA_obj)
```
The results is stored in `MeDuSA_obj@VarianceExplain`. 


## Preparing Reference Data
It is important to note that in real-world applications, users should annotate the cell-state trajectory based on their own data and research interests. There are many methods to infer the cell trajectory in scRNA-seq data, such as: 

- [Slingshot](https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html)
- [CytoTRACE](https://cytotrace.stanford.edu/)
- [Monocle3](https://cole-trapnell-lab.github.io/monocle3/)
- [scVelo](https://github.com/theislab/scvelo)

In this tutorial, we use the CytoTRACE to infer the differentiation trajectory of hPSCs.

### 1. Download Raw scRNA-seq Data
We will download the raw data from the GEO database. 
```bash
#the scRNA-seq data
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75748/suppl/GSE75748_sc_cell_type_ec.csv.gz
```
### 2.Processing Raw scRNA-seq Data
```R
library(Seurat)
library(WaveCrest)
library(data.table)
library(ggplot2)

#1. Load the count data
sce = fread('/Users/songliyang/Documents/MeDuSA_new/revision/embry/GSE75748_sc_cell_type_ec.csv.gz')
sce = as.data.frame(sce)
rownames(sce) = sce[,1];sce = sce[,-1]

#2. Infer the cell-state trajectory
cyto = CytoTRACE(sce,enableFast = F,ncores = 4)
pseudotime = cyto$CytoTRACE

#3. Build the reference scRNA-seq data
sce = CreateSeuratObject(sce)
sce$cell_type = 'embry'
sce$cell_trajectory = pseudotime
sce$sample = as.vector(Idents(sce))

# We suggest users to normalize cell-size in the reference data before running deconvolution analysis, although MeDuSA is generally robust to varying data scales
sce@assays$RNA@counts = sweep(as.matrix(sce@assays$RNA@counts),2,colSums(sce@assays$RNA@counts),'/')*1e+3
```

## Validation of the MeDuSA Method
This hPSC dataset includes both bulk RNA-seq data and scRNA-seq data from the same sample. It is expected that the cell-state abundance would strongly correlate between the two types of data, despite potential variations in the sequenced specimens. To validate the MeDuSA method, we will compare the estimated cell-state abundance from the bulk data to that measured from the scRNA-seq data.


```r
# Load the data
bulk = readRDS("../hPSC_bulk.rds")
sce = readRDS("./hPSC_sce.rds")

# Estimate cell-state abundance from the scRNA-seq data 
pseudotime = sce$cell_trajectory
sampleID = sce$sample
pseudotime = pseudotime[names(sampleID)]
abundance_expect = sapply(unique(sampleID),function(id){
	pseudotime_temp = sort(pseudotime[which(sampleID == id)])
	pseudotime_temp = pseudotime_temp
	count_temp = hist(pseudotime_temp,breaks = seq(0,1,(1/50)))$counts
	abundance_temp  =  count_temp/sum(count_temp)
	return(abundance_temp)
})
rownames(abundance_expect) = paste0('bin',seq(1,nrow(abundance_expect)))

# Run MeDuSA
MeDuSA_obj = MeDuSA(bulk,sce,
                  select.ct = 'embry',markerGene = NULL,span = 0.35,
		  resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)	

markerGene = MeDuSA_obj@Estimation$markerGene
abundance_estimate = MeDuSA_obj@Estimation$cell_state_abundance
TimeBin = MeDuSA_obj@Estimation$TimeBin

# Visualize
commonId = intersect(colnames(abundance_expect),colnames(abundance_estimate))
abundance_expect = abundance_expect[,commonId]
abundance_estimate = abundance_estimate[,commonId]
dat = data.frame('MeDuSA' = c(abundance_estimate),'CytoTRACE' = c(abundance_expect))
p1 = ggplot(dat,aes(x=CytoTRACE,y=MeDuSA))+
  geom_point(col='#feb24c')+
  geom_smooth(method = 'lm',col='black',se=F)
print(p1)
```
