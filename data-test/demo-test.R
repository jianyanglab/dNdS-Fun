devtools::install_github("jianyanglab/dNdS-Fun")
library(dNdSFun)

dNdSFun(mutsFile="Example.GRCh37.cds.score.txt",refDb_element="CADD.cds.GRCh37.rda", reg="cds", GenoVersion="GRCh37", globaldnds_outFile="dNdS_CADD.global.out", genelevel_outFile="dNdS_CADD.element.gp.out", thread_num = 22)
