---
layout: page
title: Example (dNdSFun)
description: ~
---
This tutorial offers an illustrative analysis of the human monocytes data from [Oetjen et al., 2018](https://insight.jci.org/articles/view/124928) using MeDuSA. Prior to running the analysis, it is important to ensure that the MeDuSA package has been installed. For installation instructions, please refer to the following [link](https://jianyanglab.github.io/dNdSFun/documentation/02_installation.html).


## Input Data
`MeDuSA` requires two types of input data:
- Bulk RNA-seq data. 
- Single cell RNA-seq (scRNA-seq) data, which should be provided in the form of a Seurat object that includes the annotated cell-state trajectory and cell types. 

For how to prepare the cell-state trajectory data, please read the section of `Preparing Reference Data` in this tutorial. 

The input data required for running this tutorial can be downloaded from the following [link](https://yanglab.westlake.edu.cn/data/MeDuSA_data/Monocytes.tar.gz). 
Detailed information regarding the input data is provided as follows.

### 1. Bulk RNA-seq Data
```r
# Load the example bulk RNA-seq data
bulk = readRDS("../Monocytes_bulk.rds")
```
The bulk RNA-seq data is represented in a matrix format, where each row corresponds to a specific gene and each column corresponds to a particular sample.

### 2. Reference scRNA-seq Data
```r
# Load the example scRNA-seq data
sce = readRDS("./Monocytes_sce.rds")
class(sce)
[1] "Seurat"
attr(,"package")
[1] "SeuratObject"

sce@assays$RNA@counts[1:3,1:3]
              A/A.rds_AAACCTGCAGCGAACA-1 A/A.rds_AAACCTGGTCGACTGC-1 A/A.rds_AAACCTGGTCGCTTCT-1
FO538757.2                             .                          .                          .
AP006222.2                             .                          .                          .
RP4-669L17.10                          .                          .                          .

sce$cell_trajectory[1:3]
A/A.rds_AAACCTGCAGCGAACA-1 A/A.rds_AAACCTGGTCGACTGC-1 A/A.rds_AAACCTGGTCGCTTCT-1 
                 0.9220488                  0.5167408                  0.4567616 
		 
sce$cell_type[1:3]
A/A.rds_AAACCTGCAGCGAACA-1 A/A.rds_AAACCTGGTCGACTGC-1 A/A.rds_AAACCTGGTCGCTTCT-1 
                     "mon"                      "mon"                      "mon" 
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

For further details about the parameters, please refer to this [link](https://jianyanglab.github.io/dNdSFun/documentation/01_About.html).
```r
MeDuSA_obj = MeDuSA(bulk,sce,
                  select.ct = 'mon',markerGene = NULL,span = 0.35,
		  resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)		 
```
The results are stored in MeDuSA_obj@Estimation.
- The estimated cell-state abundance: MeDuSA_obj@Estimation$cell_state_abundance
- The median state (pseudo-time) of cell-state bins: MeDuSA_obj@Estimation$TimeBin
- The used marker genes: MeDuSA_obj@Estimation$markerGene


### 2. How to Select Marker Genes
MeDuSA provides users with the option to input their own cultivated marker genes. Additionally, MeDuSA offers two methods for selecting representative marker genes along the cell-state trajectory.

- wilcox test: MeDuSA divides the cells in the trajectory into a specified number of bins and uses the Wilcoxon rank sum test (wilcox) to identify marker genes for each bin. The Wilcoxon test is applied to compare gene expression levels between the cells in the current bin and all other cells in the trajectory. Significant genes are then identified as marker genes for that particular bin. MeDuSA performs this test for each gene at every bin along the cell-state trajectory, enabling the identification of marker genes that are specific to each stage of the trajectory. The Wilcoxon test is implemented in the Seurat::FindMarkers function.

- gam-wald test: MeDuSA employs the generalized additive model (gam) to associate genes along the cell-state trajectory and considers only genes with an FDR-adjusted p-value less than 0.01. Significant genes are ranked based on their association strength, facilitating the identification of relevant genes that are associated with the cell-state trajectory. To avoid overrepresentation of certain cell-states, MeDuSA divides the trajectory into intervals and assigns each gene to the interval where it has the highest mean expression. MeDuSA then selects a set of top informative genes as signature genes for each interval.

Users can specify the `method` in the `MeDuSA` function as either `wilcox` or `gam` to use these two methods. Alternatively, users can select the marker genes using the `MeDuSA_marker` function before running the deconvolution analysis. 

Note that if your reference scRNA-seq data was collected from multiple batches, we recommend using the 'pseudo-bulk' method, adjusting for the batches as covariates, to identify marker genes. In this scenario, the pseudo-bulk method can enhance the robustness of the marker genes.

```r
# Set the gene selection method in MeDuSA function 
library(MeDuSA)
# The wilcox test
MeDuSA_obj = MeDuSA(bulk,sce,
		    select.ct = 'mon',markerGene = NULL,span = 0.35,method = "wilcox",
		    resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)	
# The gam-wald test
MeDuSA_obj = MeDuSA(bulk,sce,
		    select.ct = 'mon',markerGene = NULL,span = 0.35,method = "gam",
		    resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)
		    
# Select marker genes using MeDuSA_marker before running deconvolution analysis
marker = MeDuSA_marker(sce[,which(sce$cell_type=='mon')],bulk,
                               geneNumber = 200,nbins = 10,
                               family ="gaussian",k = 10,ncpu = 4,method = "wilcox")
# Documentations			       
help(MeDuSA_marker)			       
```
### 3. How to Include Other Cell Types as Covariates
MeDuSA offers two ways to incorporate other cell types as covariates to address confounding factors. 

- To save memory during the deconvolution analysis, we recommend that users construct the covariates matrix before running the analysis.

```r
# Load the data
sce_otherCT = readRDS("../Monocytes_OtherCell.rds")
cov_otherCT = Seurat::AverageExpression(object = sce_otherCT,group.by = 'cell_type',assays ='RNA',slot='counts')$RNA
remove(sce_otherCT)

# To input the covariates matrix into MeDuSA, users can specify the parameter of fixCov. 
MeDuSA_obj = MeDuSA(bulk,sce,fixCov = cov_otherCT,
		    select.ct = 'mon',markerGene = NULL,span = 0.35,method = "wilcox",
		    resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)	
```

- Alternatively, users can merge the data of the focal cell-type and other cell-types into a single Seurat object.

```r
# Load the data
sce_otherCT = readRDS("../Monocytes_OtherCell.rds")
sce_big = merge(sce,sce_otherCT)
cell_type = c(sce$cell_type,sce_otherCT$cell_type)
sce_big$cell_type = cell_type[colnames(sce_big)]
remove(sce_otherCT); remove(sce)

# MeDuSA automatically constructs the covariates matrix based on the cell-type labels stored in sce_big$cell_type.
MeDuSA_obj = MeDuSA(bulk,sce = sce_big,
		    select.ct = 'mon',markerGene = NULL,span = 0.35,method = "wilcox",
		    resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)	
```

### 4. How to Use the Mode of Conditional Autoregressive (CAR)
In many LMM applications, random effects are assumed to be independent and identically distributed. However, the abundances of cells at adjacent states are likely to be correlated. MeDuSA incorporates the CAR to accommodate such correlations. Users can set the parameter of `CAR` as TRUE to use the CAR mode. Further, MeDuSA provides users with the ability to set the range of possible correlation strengths through the `phi` parameter.  Based on this range, MeDuSA then automatically searches for the optimal correlation strength that maximizes the likelihood function. 
```R
# The default phi is [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
MeDuSA_obj = MeDuSA(bulk,sce, CAR = TRUE,
                   select.ct = 'mon',markerGene = NULL,span = 0.35,
		   resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)	

# Change the parameter space of phi
phi = c(0,0.1,0.99)
MeDuSA_obj = MeDuSA(bulk,sce, CAR = TRUE, phi = phi,
                   select.ct = 'mon',markerGene = NULL,span = 0.35,
		   resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)	
``` 
It's worth noting that when using the CAR mode in MeDuSA, the computational speed may slow down due to the computational burden of inverting the covariance matrix between cells. This can become particularly significant when working with large reference datasets.

### 5. How to Normalize the Data
It's important to normalize the reference and bulk data to the same scale before running deconvolution analysis using MeDuSA. MeDuSA does not perform normalization on the input data because the data may be in various scales, such as raw counts, CPM, TPM, FPKM, or log-transformed data. While MeDuSA is generally robust to different scales, heterogeneity in data scale between the bulk and reference data may negatively impact performance. Therefore, users must carefully check and normalize their data appropriately before running MeDuSA to ensure accurate and reliable results. For instance, in this tutorial, we normalized the data to CPM scale.
```r
bulk = sweep(bulk,2,colSums(bulk),'/')*1e+4
sce@assays$RNA@counts = sweep(sce@assays$RNA@counts,2,colSums(sce@assays$RNA@counts),'/')*1e+4
```

### 6. How to Obtain the P-Value of the Random Effects Component
Once the deconvolution analysis is complete, users can use the `MeDuSA_VarExplain` function to obtain the explained variance of the bulk data by the reference scRNA-seq data, as well as the corresponding p-values.
```R
MeDuSA_obj = MeDuSA_VarExplain(MeDuSA_obj)
```
The results is stored in `MeDuSA_obj@VarianceExplain`. 


## Preparing Reference Data
In this section, we will walk through the steps involved in preparing the reference scRNA-seq data used in this tutorial. 
### 1. Downloading Raw scRNA-seq Data
We download the raw data from the GEO database and rename them based on their respective sample names
```bash
#!/bin/bash
#1. Download the data from the GEO database
mkdir JCI
cd JCI
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120221/suppl/GSE120221_RAW.tar
tar -xvf GSE120221_RAW.tar

#2. Rename the data based on sample id 
ls *mtx.gz | while read file
do
  id=$(echo $file | cut -d "_" -f 3 | cut -d "." -f 1)
  mkdir `pwd`/${id}
  mv *barcodes_${id}.tsv.gz `pwd`/${id}/barcodes.tsv.gz 
  mv *genes_${id}.tsv.gz `pwd`/${id}/features.tsv.gz 
  mv *matrix_${id}.mtx.gz `pwd`/${id}/matrix.mtx.gz 
done
```
### 2. Processing Raw scRNA-seq Data
We use [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) to detect possible doublets in each scRNA-seq data. After filtering out the identified doublets, the data will be merged into a single Seurat object.

```R
library(Seurat)
library(DoubletFinder)
library(dplyr)

#1. Read the data and do quality control
setwd("../JCI")
file = list.files()
for(id in file){
	print(id)
# Read the data
	path = paste0("../JCI/",id)
	data = Read10X(data.dir = path)
	data = CreateSeuratObject(counts = data,min.cells = 3, min.features = 200)
	data[["percent.mt"]] = PercentageFeatureSet(data, pattern = "^MT-")
	
# Standard process
	data = NormalizeData(data) %>%
	       FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
	       ScaleData() %>%
	       RunPCA() %>%
	       RunUMAP(dims = 1:15) %>%
	       FindNeighbors(reduction = "pca", dims = 1:15) %>%
	       FindClusters() %>%
	
# Define the doublet rate (follow the possible doublet rate provided by 10x)
	ncell = ncol(data)
	if(ncell<500){dbrate = 0.4/100}
	if(ncell>=500 && ncell<1000){dbrate = 0.4/100}
	if(ncell>=1000 && ncell<2000){dbrate = 0.8/100}
	if(ncell>=2000 && ncell<3000){dbrate = 1.6/100}
	if(ncell>=3000 && ncell<4000){dbrate = 2.3/100}
	if(ncell>=4000 && ncell<5000){dbrate = 3.1/100}
	if(ncell>=5000 && ncell<6000){dbrate = 3.9/100}
	if(ncell>=6000 && ncell<7000){dbrate = 4.6/100}
	if(ncell>=7000 && ncell<8000){dbrate = 5.4/100}
	if(ncell>=8000 && ncell<9000){dbrate = 6.1/100}
	if(ncell>=9000 && ncell<10000){dbrate = 6.9/100}
	if(ncell>=10000){dbrate = 7.6/100}

# Remove doublet
	homotypic.prop = DoubletFinder::modelHomotypic(Idents(data))  
	nExp_poi = round(dbrate*ncell) 
	nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))
	data = DoubletFinder::doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	
# Save the data
	out = paste(path,paste0(id,'.rds'),sep='/')
	saveRDS(data,out)		
}

#2. Merge the data
setwd("../JCI")
file = list.files()
data = sapply(file,function(id){
  print(id)
  path = paste(id,paste0(id,'.rds'),sep='/')	
  dat = readRDS(path)
  db = c(table(dat@meta.data[,grep('DF',colnames(dat@meta.data))])[1])
  dat.qc = dat[,dat@meta.data[,grep('DF',colnames(dat@meta.data))]=='Singlet']
  dat.qc = dat.qc[,dat.qc$percent.mt<20]
  dat.qc = RenameCells(dat.qc,new.names=paste(path,colnames(dat.qc),sep='_'))
  dat.qc$sample = id
  dat.qc 
})
data.merge = Reduce(function(x, y) merge(x, y), data)
saveRDS(data.merge,'../Human_BoneMarrow_JCI_Insight.rds')
```
The merged scRNA-seq data can be obtained either by following the above pipeline or by downloading it directly from the [link] (https://yanglab.westlake.edu.cn/data/MeDuSA_data/Monocytes.tar.gz)

### 3. Cell Type Annotation
We will perform cell clustering and assign cell-types based on expression pattern of marker genes. To account for potential confounding factors during single-cell RNA sequencing, the black gene list, provided by [Xue et al.](https://www.nature.com/articles/s41586-022-05400-x) will be utilized. 
```R
library(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggplot2)

# Load the data
BlackGene = read.csv('../Gene_BlackList.csv',fill=T)
BM = readRDS('../Human_BoneMarrow_JCI_Insight.rds')

# Rename the sample C1 as C
BM$sample[which(BM$sample=='C1')]='C'

# Construct the confounding score
mt_gene = intersect(unique(BlackGene[,'Mitochondria']),rownames(BM))
hsp_gene = intersect(unique(BlackGene[,'Heat.shock.protein']),rownames(BM))
rib_gene = intersect(unique(BlackGene[,'Ribosome']),rownames(BM))
disso_gene = intersect(unique(BlackGene[,'Dissociation']),rownames(BM))
noisy_gene = c(mt_gene,hsp_gene,rib_gene,disso_gene)
BM = AddModuleScore(BM,features=list(mt_gene),name = 'mt_score',nbin=10)
BM = AddModuleScore(BM,features=list(hsp_gene),name = 'hsp_score',nbin=10)
BM = AddModuleScore(BM,features=list(rib_gene),name = 'rib_score',nbin=10)
BM = AddModuleScore(BM,features=list(disso_gene),name = 'disso_score',nbin=10)

# Regress out the confounding score
out = c('mt_score1','hsp_score1','rib_score1','disso_score1')
BM = BM[!rownames(BM) %in% noisy_gene,]
BM = NormalizeData(BM) %>% 
     FindVariableFeatures() %>% 
     ScaleData(vars.to.regress = out) %>% 
     RunPCA(verbose = FALSE) %>%
     RunUMAP(reduction = "pca", dims = 1:50) %>% 
     FindNeighbors(reduction = "pca", dims = 1:50) %>%
     FindClusters(resolution=2)

mk_gene = c('HBD','AVP','CD3E','CD3G','CD3D','CD79B','CD79A','NKG7','GATA1','SPI1','CD4','CD8A','FOXP3','CCR7','CD74','TMSB10')

# Check the expression level of marker genes
figure.dim = DimPlot(BM,pt.size=0.05,label=T,label.size = 6,repel = TRUE)+ 
               theme(legend.position='none') +
               ggtitle(NULL) +
               xlab('UMAP-1')+ylab('UMAP-2')

p = FeaturePlot(BM,features = mk_gene,combine = FALSE)
for(i in 1:length(p)) {p[[i]] = p[[i]] + NoLegend() + NoAxes()}
figure.feature  = cowplot::plot_grid(plotlist = p,ncol = 4)

# Annotate cell types
HSPC = c('21')
HSPCtoMon = c('22','23','26','30','1','16')
ct = as.vector(Idents(BM))
ct[ct %in% HSPC] = 'HSPC'
ct[ct %in% HSPCtoMon] = 'HSPCtoMon'
BM$ct = ct

# Splict the data based on cell types and save the data. 
Mon = BM[,BM$ct %in% c('HSPC','HSPCtoMon')]
OtherCell = BM[,BM$ct %in% setdiff(cell_type_annotated,c('HSPC','HSPCtoMon'))]
saveRDS(Mon,'../Mon.rds')
saveRDS(OtherCell,'../OtherCell.rds')
```
For user convenience, we have provided the processed scRNA-seq data in the link above.

### 4. Infering the Cell-State Trajectory
In this analysis, we use [Slingshot](https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html) to infer the development trajectory of monocytes. 

```R
library(Seurat)
library(dplyr)
library(slingshot)
library(ggplot2)

# Load the data
Mon = readRDS('../Mon.rds')

# Regress out the confounding factors and use the SCT transformation. 
confounding = c("mt_score1","hsp_score1","rib_score1","disso_score1","nCount_RNA","sample")
Mon = NormalizeData(Mon) %>%
      SCTransform(vars.to.regress = confounding,return.only.var.genes = T) %>%
      RunPCA(verbose = FALSE) %>%
      RunUMAP(Mon, reduction = "pca", dims = 1:15) %>%
      FindNeighbors(Mon, reduction = "pca", dims = 1:15)  %>%
      FindClusters(Mon,resolution = 2)
	
# Remove clusters that are not connect to the main trajectory.    
sce = Mon[,!Idents(Mon) %in% c('18','12','27','2','25')]
FeaturePlot(Mon,features=c('CD34','FGL2'))
remove(Mon)

# Annotate cell trajectory using the slingshot
umap = as.data.frame(Embeddings(sce,'umap'))
umap$ct = as.vector(Idents(Mon))
sds = getLineages(umap[,1:2], umap$ct, start.clus = '21')
sds = getCurves(sds)
path = as.data.frame(sds@metadata$curves$Lineage1$s)
pseudo_time = sds@metadata$curves$Lineage1$lambda
pseudo_time = pseudo_time/max(pseudo_time)

# Add the cell-type and cell-state trajectory into the meta data. 
sce$cell_type = 'mon'
sce$cell_trajectory = pseudo_time[colnames(sce)]

# Visualize the cell-state trajectory
colors = c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2")
colors = colors[seq(length(colors),1,-1)]
p1 = ggplot(umap,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(col=pseudo_time),size=0.5)+
  xlab('PATH-1')+
  ylab('PATH-2')+
  theme(legend.position = 'right',
        legend.justification = "left",
        panel.border = element_blank(), axis.line = element_line(size = 0.8))+
  annotate('text',x=-2,y=-9,label='GMPs',size=5)+
  annotate('text',x=3,y=-3,label='Non-classical monocytes',size=5)+
  scale_color_gradientn(colours = colors,name='Pseudotime',labels = scales::number_format(accuracy = 0.1))
p1

# Save the reference scRNA-seq data
saveRDS(sce,'../Monocytes_sce.rds')
```

## Validation of the MeDuSA Method
This dataset includes both bulk RNA-seq data and scRNA-seq data from the same sample. It is expected that the cell-state abundance would strongly correlate between the two types of data, despite potential variations in the sequenced specimens. To validate the MeDuSA method, we will compare the estimated cell-state abundance from the bulk data to that measured from the scRNA-seq data.

```r
bulk = readRDS("../Monocytes_bulk.rds")
sce = readRDS("./Monocytes_sce.rds")

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
abundance_expect = abundance_expect[, colnames(abundance_expect) %in% colnames(bulk)]

# Compare the cell-state abundance estimated from bulk data to that from scRNA-seq data
abundance_estimate = MeDuSA_obj@Estimation$cell_state_abundance
state = MeDuSA_obj@Estimation$TimeBin
rownames(abundance_estimate) = state
rownames(abundance_expect) = state

commonId = intersect(colnames(abundance_expect),colnames(abundance_estimate))
abundance_expect = abundance_expect[,commonId]
abundance_estimate = abundance_estimate[,commonId]
dat = data.frame('MeDuSA' = c(abundance_estimate),'Slingshot' = c(abundance_expect))
p2 = ggplot(dat,aes(x=Slingshot,y=MeDuSA))+
  geom_point(col='#feb24c')+
  geom_smooth(method = 'lm',col='black',se=F)+
  xlab('scRNA-seq (slingshot)')
  ylab('MeDuSA')
print(p2)
```

