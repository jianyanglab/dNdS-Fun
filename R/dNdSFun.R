
dNdSFun <- function(mutsFile,refDb_element, reg, positive, negative, positiveThreshold, negativeThreshold, gene_group, globaldnds_outFile,
                 genelevel_selloc_outFile, genelevel_selcv_outFile, iscv, score = "ture", model = 3, negmu = 1){
    rm(list=ls())
    source("CADD_dndsWGS.NEG.R")
    # .libPaths("/storage/yangjianLab/zhengmengyue/SOFTWARE.bak/R_LIB_4.0.5")
    library(parallel)
    library(data.table)
    library(MASS)
    library(doParallel)
    positiveThreshold <- as.numeric(positiveThreshold)
    negativeThreshold <- as.numeric(negativeThreshold)
    negbeta <- as.numeric(negmu)
    model <- as.character(model)
    elements_list=NULL
    cv=NULL
    max_muts_per_element_per_sample=3
    max_element_muts_per_sample=100000000
    maf_data <- fread(mutsFile,header = FALSE)
    if(model=="1"){
        outp=1
    }else{
        outp=3
    }

    # Group by the second column
    grouped_data <- split(maf_data, maf_data[, 2])
    sorted_keys <- names(grouped_data)[order(as.numeric(gsub("chr", "", names(grouped_data))))]
    data_sorted <- grouped_data[sorted_keys]
    group_array <- lapply(data_sorted, identity)

    # Load reference data of noncoding elements
    load(refDb_element) # gr_elements, RefElement
    data_classified <- split(RefElement, unlist(lapply(RefElement, "[[", "chr")))
    sorted_keys <- names(data_classified)[order(as.numeric(gsub("chr", "", names(data_classified))))]
    data_sorted <- data_classified[sorted_keys]
    RefElement_array <- lapply(data_sorted, identity)
    thread_data <- Map(function(x, y) list(x, y), group_array, RefElement_array)
    score <- toupper(score)

    # Define trinucleotide contexts
    nt <- c("A","C","G","T")
    trinucs <- paste(rep(nt,each=16,times=1),rep(nt,each=4,times=4),rep(nt,each=1,times=16), sep="")
    trinucidxs <- setNames(1:64, trinucs)
    trinucMuts <- NULL
    for (j in 1:length(trinucs)) {
    trinucMuts <- c(trinucMuts, paste(trinucs[j], paste(substr(trinucs[j],1,1), setdiff(nt,substr(trinucs[j],2,2)), substr(trinucs[j],3,3), sep=""), sep=">"))
    }
    trinucMutsidx <- setNames(1:192, trinucMuts)

    #Thread function
    ncpu <- length(thread_data)
    ncpu = min(ncpu,parallel::detectCores())
    cl = parallel::makeCluster(ncpu)
    parallel::clusterExport(cl=cl, varlist=c("gr_elements", "nt", "score", "trinucMutsidx", "elements_list", "cv", "max_muts_per_element_per_sample",
                                            "max_element_muts_per_sample", "positiveThreshold", "negativeThreshold"), envir=environment())
    registerDoParallel(cl)
    `%dopar2%` = foreach::`%dopar%`

    result = foreach::foreach(data = thread_data) %dopar2% {
        group <- data[[1]]
        RefElement <- data[[2]]
        if (score == "FALSE"){
            # STEPT2 Add score to input
            chr_name <- group[1,2]
            # Open log file
            if (!file.exists("tmp")) {
                dir.create("tmp")
            }
            log_file_path <- paste0(chr_name,"log.txt")
            file_path <- file.path("tmp", log_file_path)
            log_file <- file(file_path, open = "wt")
            sink(log_file, type = "message")
            maf <- data.frame()
            for (data_row in 1:nrow(group)) {
                row_data <- group[data_row,]
                first_four_columns <- row_data[,2:5]
                input_data_comparison <- paste(first_four_columns, collapse = "")
                chr_num <- gsub("\\D", "", row_data[,2])
                postion <- row_data[,3]
                commend <- paste0("tabix -h /storage/yangjianLab/sunxiwei/data/annotations/CADD/chr/whole_genome_SNVs.tsv.gz.",chr_name,".gz.rankRawScore.gz ",chr_num,":",as.integer(postion),"-",as.integer(postion),"| grep -v \"^#\"")
                database <- system(commend, intern = TRUE)
                count <- 0
                for (database_row in database) {
                    # row_database <- database[data_row]
                    split_database_row <- strsplit(database_row, "\t")[[1]]
                    combined_database <- paste(split_database_row[1:4], collapse = "")
                    combined_database <- paste0("chr",combined_database)
                    # print(combined_database)

                    if (combined_database == input_data_comparison){
                        input_result = c(row_data,split_database_row[5])
                        colnames(maf) <- colnames(input_result)
                        maf <- rbind(maf,input_result)
                        count <- count + 1

                    }
                }
                if (count == 0){
                    message("Failed to match",row_data," in database!")
                }
            }
            sink(type = "message")
            close(log_file)
            maf <- na.omit(maf)
        }else{
            maf <- na.omit(group)
        }

        # Step 1: Variables required
        message("Step 1: Loading the variables required ......")
        # [Input] Gene list (The user can input a gene list as a character vector)
        if (is.null(elements_list)) {
            elements_list = sapply(RefElement, function(x) x$gene_name) # All genes [default]
        } else { # Using only genes in the input gene list
            allg = sapply(RefElement,function(x) x$gene_name)
            nonex = elements_list[!(elements_list %in% allg)]
            if (length(nonex)>0) { stop(sprintf("The following input gene names are not in the RefCDS database: %s", paste(nonex,collapse=", "))) }
            RefElement = RefElement[allg %in% elements_list] # Only input genes
            gr_genes = gr_genes[gr_genes$names %in% elements_list] # Only input genes
        }
        maf <- maf[,1:6] #  Requiring the first 6 columns of input data
        maf[,c(1,2,3,4,5)] <- lapply(maf[,c(1,2,3,4,5)], as.character) # Convert factors to character
        maf[[3]] <- as.numeric(maf[[3]]) # Convert position as numeric
        maf[[6]] <- as.numeric(maf[[6]]) # Convert impact scores as numeric
        maf <- maf[maf[,4]!=maf[,5],] # Remove mutations with identical reference and mutant base
        colnames(maf) <- c("sampleID","chr","pos","ref","alt","impScore")
        # Remove rows with NA values
        idxna <- which(is.na(maf),arr.ind=TRUE)
        if (nrow(idxna)>0) {
            maf <- maf[-unique(idxna[,1]),]
            warning(sprintf("%0.0f rows in the input table contained NA entries and have been removed. Please investigate.",length(unique(idxna[,1]))))
        }

        # Expanding the reference sequences [for faster access]
        for (j in 1:length(RefElement)) {
            RefElement[[j]]$seq_element <- base::strsplit(as.character(RefElement[[j]]$seq_element), split="")[[1]]
            RefElement[[j]]$seq_element1up <- base::strsplit(as.character(RefElement[[j]]$seq_element1up), split="")[[1]]
            RefElement[[j]]$seq_element1down <- base::strsplit(as.character(RefElement[[j]]$seq_element1down), split="")[[1]]
        }

        idx <- setNames(1:length(RefElement), sapply(RefElement,function(x) x$gene_name))
        gr_elements_idx <- idx[gr_elements$names]
        # Start and end position of each mutation
        maf$end <- maf$start <- maf$pos
        # Mapping mutations to genes
        gr_muts <- GenomicRanges::GRanges(maf$chr, IRanges::IRanges(maf$start,maf$end))
        ol <- as.data.frame(GenomicRanges::findOverlaps(gr_muts, gr_elements, type="any", select="all"))
        maf <- maf[ol[,1],] # Duplicating subs if they hit more than one gene
        maf$elementidx <- gr_elements_idx[ol[,2]]
        maf$element <- sapply(RefElement,function(x) x$gene_name)[maf$elementidx]

        #maf <- unique(maf)
        # Excluding samples with more than the defined maximum number of mutations per sample (optional)
        nsampl <- sort(table(maf$sampleID))
        exclsamples <- NULL
        if (any(nsampl>max_element_muts_per_sample)) {
            message(sprintf('    Note: %0.0f samples excluded for more than the defined maximum number of mutations per sample (see the max_element_muts_per_sample argument). %0.0f samples left after filtering.',sum(nsampl>max_element_muts_per_sample),sum(nsampl<=max_element_muts_per_sample)))
            exclsamples <- names(nsampl[nsampl>max_element_muts_per_sample])
            maf <- maf[!(maf$sampleID %in% names(nsampl[nsampl>max_element_muts_per_sample])),]
        }
        exclmuts <- NULL
        maf$strand <- sapply(RefElement,function(x) x$strand)[maf$elementidx]
        snv <- (maf$ref %in% nt & maf$alt %in% nt)
        if (!any(snv)) { stop("Zero coding substitutions found in this dataset. Common causes for this error are inputting only indels or using chromosome names different to those in the reference database (e.g. chr1 vs 1)") }
        maf <- maf[snv,]
        maf$impidx <- rep(NA,dim(maf)[1])
        maf[maf$impScore<as.numeric(negativeThreshold),"impidx"] <- 1 # mutations group: neutral set
        ##mid set
        maf[maf$impScore>=as.numeric(positiveThreshold),"impidx"] <- 2 # mutations group: selection target set
        idxna <- which(is.na(maf),arr.ind=TRUE)
        if (nrow(idxna)>0) {
            maf <- maf[-unique(idxna[,1]),]
            warning(sprintf("%0.0f mutations contained NA impact score values and have been removed. Please investigate.",length(unique(idxna[,1]))))
        }
        maf$ref_cod <- maf$ref
        maf$mut_cod <- maf$alt
        compnt <- setNames(rev(nt), nt)
        isminus <- (maf$strand==-1 | maf$strand=="-")
        maf$ref_cod[isminus] <- compnt[maf$ref[isminus]]
        maf$mut_cod[isminus] <- compnt[maf$alt[isminus]]
        for (j in 1:length(RefElement)) {
            RefElement[[j]]$N = array(0, dim=c(192,2)) # Initialising the N matrices
        }

        # Subfunction: obtaining the positions of a noncoding mutation given the intervals of the elements
        chr2element <- function(pos,element_int,strand) {
            if (strand==1 || strand=="+") {
            return(which(unlist(apply(element_int, 1, function(x) x[1]:x[2])) %in% pos))
            } else if (strand==-1 || strand=="-") {
            return(which(rev(unlist(apply(element_int, 1, function(x) x[1]:x[2]))) %in% pos))
            }
        }

        ref3_cod <- mut3_cod <- wrong_ref <- impact <- codonsub <- array(NA, nrow(maf))

        for (j in 1:nrow(maf)) {
            elementidx <- maf$elementidx[j]
            pos <- maf$pos[j]
            impidx <- maf$impidx[j]
            pos_idx <- chr2element(pos, RefElement[[elementidx]]$intervals_element, RefElement[[elementidx]]$strand)
            element <- RefElement[[elementidx]]$seq_element[pos_idx]
            ref3_cod[j] <- sprintf("%s%s%s", RefElement[[elementidx]]$seq_element1up[pos_idx], RefElement[[elementidx]]$seq_element[pos_idx], RefElement[[elementidx]]$seq_element1down[pos_idx])
            mut3_cod[j] <- sprintf("%s%s%s", RefElement[[elementidx]]$seq_element1up[pos_idx], maf$mut_cod[j], RefElement[[elementidx]]$seq_element1down[pos_idx])

            if (maf$ref_cod[j] != as.character(element)) { # Incorrect base annotation in the input mutation file (the mutation will be excluded with a warning)
            wrong_ref[j] <- 1
            } else if (!is.na(impidx)) { # Correct base annotation in the input mutation file
            triMut <- trinucMutsidx[ paste(ref3_cod[j], mut3_cod[j], sep=">") ]
            RefElement[[elementidx]]$N[triMut,impidx] <- RefElement[[elementidx]]$N[triMut,impidx] + 1 # Adding the mutation to the N matrices
            # RefElement_addN <- append(RefElement_addN, list(RefElement[[elementidx]]))
            }

            # if (round(j/1e4)==(j/1e4)) { message(sprintf('    %0.3g%% ...', round(j/nrow(maf),2)*100)) }
        }

        wrong_refbase <- NULL
        if (any(!is.na(wrong_ref))) {
            if (mean(!is.na(wrong_ref)) < 0.1) { # If fewer than 10% of mutations have a wrong reference base, we warn the user
            warning(sprintf('%0.0f (%0.2g%%) mutations have a wrong reference base (see the affected mutations in dndsout$wrongmuts).', sum(!is.na(wrong_ref)), 100*mean(!is.na(wrong_ref))))
            } else { # If more than 10% of mutations have a wrong reference base, we stop the execution (likely wrong assembly or a serious problem with the data)
            stop(sprintf('%0.0f (%0.2g%%) mutations have a wrong reference base. Please confirm that you are not running data from a different assembly or species.', sum(!is.na(wrong_ref)), 100*mean(!is.na(wrong_ref))))
            }
            wrong_refbase <- maf[!is.na(wrong_ref), 1:6]
            maf <- maf[is.na(wrong_ref),]
        }
        return(list(maf = maf, RefElement1 = RefElement, exclsamples = exclsamples))

    }

    # Closing thread cluster
    stopCluster(cl)
    maf_result <- lapply(result, function(x) x$maf)
    maf <- do.call(rbind, maf_result)
    RefElement1 <- unlist(lapply(result, function(x) x$RefElement1), recursive = FALSE)
    # RefElement <- unlist(RefElement_result, recursive = FALSE)
    exclsamples <- unlist(lapply(result, function(x) x$exclsamples))

    if (score == "FALSE"){
        output_logfile <- file("log.txt", "w")
        file_list <- list()
        for (log_num in 1:22){
            logfile_name <- paste0("chr",log_num,"log.txt")
            file_list <- c(file_list,logfile_name)
        }

        for (file_name in file_list) {
            file_path <- file.path("tmp", file_name)
            file <- file(file_path, "r")
            cat(readLines(file), file = output_logfile, sep = "\n")
            file.remove(file_path)
            close(file)
        }
        close(output_logfile)
    }

    CADD_dndsWGSout <- dnds2wgs.noncoding(maf, RefElement1, exclsamples, negbeta, trinucMuts, outp)

    if(model=="1"){
        globaldnds <- CADD_dndsWGSout$globaldnds
        ####output of globaldnds
        info <- c(positive,negative,positiveThreshold,negativeThreshold,gene_group,negmu,iscv)
        globaldnds_res <- as.data.frame(matrix(c(unlist(globaldnds),info),nrow=1))
        colnames(globaldnds_res) <- c("mle","ci_low","ci_high","AIC","deviance",
                                    "overdis_chisq","overdis_ratio","overdis_p",
                                    "mle_qua","ci_low_qua","ci_high_qua",
                                    "positive","negative","positiveThreshold","negativeThreshold",
                                    "gene_group","negbeta","iscv")
        write.table(globaldnds_res,globaldnds_outFile,sep="\t",row.names=F,quote=F)

    }else{
        globaldnds <- CADD_dndsWGSout$globaldnds
        sel_cv <- CADD_dndsWGSout$sel_cv
        sel_loc <- CADD_dndsWGSout$sel_loc
        annotmuts <- CADD_dndsWGSout$annotmuts
        genemuts <- CADD_dndsWGSout$genemuts
        mle_submodel <- CADD_dndsWGSout$mle_submodel
        nbreg <- CADD_dndsWGSout$nbreg
        nbregind <- CADD_dndsWGSout$nbregind
        possmodel <- CADD_dndsWGSout$poissmodel

        ####output of globaldnds
        info <- c(positive,negative,positiveThreshold,negativeThreshold,gene_group,negmu,iscv)
        globaldnds_res <- as.data.frame(matrix(c(unlist(globaldnds),info),nrow=1))
        colnames(globaldnds_res) <- c("mle","ci_low","ci_high","AIC","deviance",
                                    "overdis_chisq","overdis_ratio","overdis_p",
                                    "mle_qua","ci_low_qua","ci_high_qua",
                                    "positive","negative","positiveThreshold","negativeThreshold",
                                    "gene_group","negbeta","iscv")
        write.table(globaldnds_res,globaldnds_outFile,sep="\t",row.names=F,quote=F)

        ####output of sel_loc
        selloc_res <- as.data.frame(sel_loc)
        write.table(selloc_res,genelevel_selloc_outFile,sep="\t",row.names=F,quote=F)

        ####output of sel_cv
        selcv_res <- as.data.frame(sel_cv)
        write.table(selcv_res,genelevel_selcv_outFile,sep="\t",row.names=F,quote=F)
    }

}


