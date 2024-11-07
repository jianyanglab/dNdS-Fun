# dNdS-Fun: Genome-wide dN/dS framework for somatic detection in cancer genomes

**dNdS-Fun** is a genome-wide `selection` detection algorithm, integrating established genome-wide functional impact scores into the conventional dN/dS framework for cancer genome study. **dNdS-Fun** significantly improves our understanding of selection in noncoding regions of cancer genomes, which account for more 98.5% of the genome than coding sequences. Similar to synonymous sites in dN/dS, **dNdS-Fun** define variants at the bottom 50% functional impact scores in whole genome, most of which are assumed selectively neutral, as nonfunctional class of sites to control background mutation rates; as nonsynonymous sites, other variants at the top 50% functional impact scores are defined as functional class of sites. Then, selection can be quantified as the ratio between the probability of a mutation occurring at either class of sites. To correct context-dependent effects of mutations, we also fit all 192 trinucleotide mutational types (all possible combinations for one base upstream and downstream from the mutant base in either transcribed or non-transcribed strand) in the model as previously described.

## Installation
### Work with R packages
```R

# Install the dNdS-Fun (R version > 4.0.5)

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jianyanglab/dNdS-Fun", build_vignettes=F)
```
### Work with code
  ####   Requirements
```R
# Install the dNdS-Fun (R version > 4.0.5)

1.parallel
2.data.tabl
3.MASS
4.doParallel
5.foreach
6.GenomicRanges (install.packages("BiocManager") -> BiocManager::install("GenomicRanges"))
7.Biostrings (BiocManager::install("Biostrings"))
8.tabix
```

## How to Use 
See [tutorial.](https://jianyanglab.github.io/dNdS-Fun/)


## Contact
If you have any questions for dNdS-Fun, please feel free to leave messages on the github issues or contact us <zhengmengyue@westlake.edu.cn> or <houjunren@westlake.edu.cn>.   


## Citation
