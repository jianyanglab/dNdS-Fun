---
layout: full
homepage: true
disable_anchors: true

---
## dNdSFun
![iDEA\_pipeline](dNdS_Overview.png)

### Introductions
dNdSFun is a genome-wide `selection` detection algorithm, integrating established genome-wide functional impact scores into the conventional dN/dS framework for cancer genome study. dNdSFun significantly improves our understanding of selection in noncoding regions of cancer genomes, which account for more 98.5% of the genome than coding sequences. Similar to synonymous sites in dN/dS, dNdSFun define variants at the bottom 50% functional impact scores in whole genome, most of which are assumed selectively neutral, as nonfunctional class of sites to control background mutation rates; as nonsynonymous sites, other variants at the top 50% functional impact scores are defined as functional class of sites. Then, selection can be quantified as the ratio between the probability of a mutation occurring at either class of sites. To correct context-dependent effects of mutations, we also fit all 192 trinucleotide mutational types (all possible combinations for one base upstream and downstream from the mutant base in either transcribed or non-transcribed strand) in the model as previously described.



### Installations
* [here](https://jianyanglab.github.io/dNdSFun/documentation/02_installation.html)

### Example analysis
* [Example](https://jianyanglab.github.io/dNdSFun/documentation/04_Example.html)

