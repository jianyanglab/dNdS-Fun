---
layout: page
title: Data
description: Link to download the RefElement dataset and Score dataset
---
The dNdS-Fun sub-model 
-------------------
### 1. RefElement
  * [GRCh37: whole-gene](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh37/CADD.whole-gene.GRCh37.rda)  (78.0MB) 
  * [GRCh37: cds-exon](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh37/CADD.cds-exon.GRCh37.rda)  (36.5MB)
  * [GRCh37: ss](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh37/CADD.ss.GRCh37.rda)  (18.4MB)
  * [GRCh37: 5utr](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh37/CADD.5utr.GRCh37.rda)  (15.6MB)
  * [GRCh37: 3utr](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh37/CADD.3utr.GRCh37.rda)  (27.4MB)
  * [GRCh37: prom](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh37/CADD.prom.GRCh37.rda)  (30.0MB)
  * [GRCh37: intron-gencode](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh37/CADD.intron-gencode.GRCh37.rda)  (1.2GB)
  * [GRCh37: intergenic-gencode](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh37/CADD.intergenic-gencode.GRCh37.rda)  (813MB)
  * [GRCh37: DichotomyLog](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh37/Dichotomy.GRCh37.log)  (259B)
  <br>
  * [GRCh38: whole-gene](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh38/CADD.whole-gene.GRCh38.rda)  (78.0MB) 
  * [GRCh38: cds-exon](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh38/CADD.cds.GRCh38.rda)  (36.5MB)
  * [GRCh38: ss](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh38/CADD.ss.GRCh38.rda)  (18.4MB)
  * [GRCh38: 5utr](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh38/CADD.5utr.GRCh38.rda)  (15.6MB)
  * [GRCh38: 3utr](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh38/CADD.3utr.GRCh38.rda)  (27.4MB)
  * [GRCh38: prom](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh38/CADD.prom.GRCh38.rda)  (30.0MB)
  * [GRCh38: intron-gencode](https://jianyanglab.github.io/dNdS-Fun/documentation/03_data.html)  (ongoing)
  * [GRCh38: intergenic-gencode](https://jianyanglab.github.io/dNdS-Fun/documentation/03_data.html)  (ongoing)
  * [GRCh38: DichotomyLog](https://yanglab.westlake.edu.cn/data/dNdSFun/RefElement/GRCh38/Dichotomy.GRCh38.log)  (259B)
  
### 2. Score
  * [CADD scores (GRCh37)](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz)  (78GB)
  * [CADD scores (GRCh38)](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz)  (ongoing)






Parameters of `dNdS-Fun`
-------------------
- `mutsFile`: The input table of somatic mutations is a table-delimited text file. Each row corresponds to a somatic mutation, and six columns are separate "IndividualID, Chromosome, Position, Ref, Alt, Functional_impact_scores (raw score)". (Attention: The IndividualID is needed be provided to filtered out the duplicate mutations from the same indivudual. If you have not exact IndividualID for the recurrent mutations, you can use different pseudo-IndividualID to avoid the filteration of recurrent mutations from different individuals.)
- `refDb_element`: The trinucleotide models hava been constructed, including coding region (cds-exon), proximatory regions of splicing site (ss), 5'UTR (5utr), 3'UTR (3utr), promoter (prom), intron (intron-gencode) and intergenic (intergenic-gencode) regions.
- `reg`: A character variable to specify the corresponding constructed models, including coding region (cds-exon), proximatory regions of splicing site (ss), 5'UTR (5utr), 3'UTR (3utr), promoter (prom), intron (intron-gencode) and intergenic (intergenic-gencode) regions. Please use the values in parentheses for assignment.
- `globaldnds_outFile`: The output of global selection estimation.
(revised) - `genelevel_selcv_outFile`:
- `genelevel_outFile`: The output of gene/element-level selection estimation.
- `iscv`: The covariates parameter. The default is NULL. This parameter is designed to use covariates to reduce the uncertainty in the variation of the mutation rate across genes. However, current no covariates is currently suitable for the dNdS-Fun model. We only kept the NULL parameter until it's available.
- `score`: A bool variable for informing software whether the input mutsFile  needs to add a column named "score". This field means functional impact score for somatic mutation. Default is false (We suppose your mutsFile includes this column). If this field is missing, please set it as "true", then this parameter should be combined with the parameter "score database" (showing below).
- `score_database`: This database is used to obtain the score value for mutsFile when "socre" column is missing. If your "score" parameter is set as true, please download the [data](https://jianyanglab.github.io/dNdS-Fun/documentation/03_data.html) for input.
- `model`: A numeric variable to determine whether the gene/element-level analysis could be conducted. If the variable equal to 1, it only estimate the global selection, or it would estimate the global and gene/element-level selection.