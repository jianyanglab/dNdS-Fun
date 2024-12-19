---
layout: page
title: Installation
description: ~
---

**dNdS-Fun** is implemented as an R package, which can be installed from GitHub by:

## Installation
### Install as an R Package
**Requirements:**R version > 4.0.5   
To install dNdS-Fun from GitHub using devtools, run the following commands in R:
```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jianyanglab/dNdS-Fun", build_vignettes = FALSE)
```
### Install from Source Code
**Requirements:**    
1. parallel (version >= 4.0.5)
2. data.table (version >= 1.16.4)
3. MASS (version >= 7.3.56)
4. doParallel (version >= 1.0.17)
5. foreach (version >= 1.5.2)
6. GenomicRanges (install via Bioconductor, version >= 1.42.0)
7. Biostrings (install via Bioconductor, version >= 2.58.0)
8. **tabix** (command-line tool, version >= 0.2.6)

### Installation Steps:
1. Install necessary R packages:   
```R 
    # Install Bioconductor manager if not already installed     
    if (!requireNamespace("BiocManager", quietly = TRUE))    
         install.packages("BiocManager")   

    # Install required CRAN packages    
    install.packages(c("parallel", "data.table", "MASS", "doParallel", "foreach"))   

    # Install required Bioconductor packages     
    BiocManager::install(c("GenomicRanges", "Biostrings"))
```
2. Install **tabix**:
- On Linux: Use your package manager, e.g., sudo apt-get install tabix.
- On macOS: Use Homebrew, e.g., brew install htslib.
- On Windows: Precompiled binaries are available; ensure tabix is added to your system PATH.
3. Clone the dNdS-Fun repository and install:   
bash  
```R
git clone https://github.com/jianyanglab/dNdS-Fun.git
```
Then, from within R:   
```R
setwd("path_to_dNdS-Fun")   
install.packages(".", repos = NULL, type = "source")   
```