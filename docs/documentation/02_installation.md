---
layout: page
title: Installation
description: ~
---

`dNdSFun` is implemented as an R package, which can be installed from GitHub by:

### Dependencies 
* R version >= 4.0.5
* R packages: foreach, parallel, MASS, doParallel

#### 1. Installing Dependent Packages
```r
install.packages("devtools")
```

#### 2. Installing `dNdSFun`
```r
devtools::install_github("jianyanglab/dNdSFun", build_vignettes=F)
```

#### 3. Loading Package
```r
library(dNdSFun)
```

