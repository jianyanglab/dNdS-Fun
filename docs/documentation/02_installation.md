---
layout: page
title: Installation
description: ~
---

`dNdS-Fun` is implemented as an R package, which can be installed from GitHub by:

### Dependencies 
* R version >= 4.0.5
* R packages: foreach, parallel, MASS, doParallel

#### 1. Installing Dependent Packages
```r
install.packages("devtools")
```

#### 2. Installing `dNdS-Fun`
```r
devtools::install_github("jianyanglab/dNdS-Fun", build_vignettes=F)
```

#### 3. Loading Package
```r
library(dNdSFun)
```

