devtools::install_github("jianyanglab/dNdS-Fun")
library(dNdSFun)

dNdSFun(mutsFile="data-test/CADD.cds.GRCh37.rda",refDb_element="data-test/CADD.cds.GRCh37.rda", reg="cds", geneDB="GRCh37", globaldnds_outFile="dNdS_CADD.global.out", genelevel_outFile="dNdS_CADD.element.gp.out", thread_num = 22)
