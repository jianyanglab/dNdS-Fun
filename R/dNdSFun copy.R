dNdSFun <- function(mutsFile, refDb_element, reg, geneDB, globaldnds_outFile,
                  genelevel_selcv_outFile, iscv = NULL, score = "false", 
                  score_database = NULL, model = 2, thread_num = 22) 
{
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        BiocManager::install("GenomicRanges")
    }
    
    if (!requireNamespace("Biostrings", quietly = TRUE)) {
        BiocManager::install("Biostrings")
    }

    library(parallel)
    library(data.table)
    library(MASS)
    library(doParallel)
    library(foreach)
    library(iterators)
    if (geneDB == "GRCh38"){
        Dichotomy_path <- system.file("Dichotomy.GRCh38.log", package = "dNdSFun")
    }else{
        Dichotomy_path <- system.file("Dichotomy.GRCh37.log", package = "dNdSFun")
    }
    options_file <- read.table(Dichotomy_path, header = TRUE, stringsAsFactors = FALSE)
    positive <- options_file[options_file$Region == reg, 2]
    negative <- options_file[options_file$Region == reg, 2]
    positiveThreshold <- options_file[options_file$Region == reg, 3]
    negativeThreshold <- options_file[options_file$Region == reg, 3]
    if (is.null(positive) || is.null(negativeThreshold)){
        error_message <- "Please enter the correct reg option."
        stop(error_message)
    }
    model <- as.character(model)
    negbeta <- 1
    iscv = NULL
    elements_list=NULL
    cv=NULL
    max_muts_per_element_per_sample=3
    max_element_muts_per_sample=100000000
    maf_data <- fread(mutsFile,header = FALSE)
    maf_data <- na.omit(maf_data)
    maf_col <- ncol(maf_data)
    if (maf_col < 6 && score == "false") {
        error_message <- "Please check your input file, it should have at least 6 columns or --score should not be 'false'."
        stop(error_message)
    }
    if(model=="1"){
        outp=1
    } else if (model=="2") {
        outp=2
    } else {
        outp=3
    }

    score <- toupper(score)
    if (score == "TRUE" & is.null(score_database))
    {
      error_message <- "Please enter score_database file."
      stop(error_message)
    }

    # Group by the second column
    grouped_data <- split(maf_data, maf_data[, 2])
    chrs <- names(grouped_data)
    filt_chrs <- chrs[grepl("^(chr)?([1-9]|1[0-9]|2[0-2])$", chrs)]

    tmp_folder <- paste0(reg, "_", nrow(maf_data), "_tmp")
    if (!file.exists(tmp_folder)) {
        dir.create(tmp_folder)
    }

    # Load reference data of noncoding elements
    load(refDb_element) # gr_elements, RefElement
    data_classified <- split(RefElement, unlist(lapply(RefElement, "[[", "chr")))
    rm(RefElement)
    RefNames <- names(data_classified)
    sorted_keys <- RefNames[order(as.numeric(gsub("chr", "", RefNames)))]
    data_sorted <- data_classified[sorted_keys]
    rm(data_classified)
    RefElement_array <- lapply(data_sorted, identity)
    rm(data_sorted)
    for (idx in filt_chrs)
    {
        idx_num <- as.numeric(gsub("chr", "", idx))
        assign(paste0("RefElement_", idx), RefElement_array[[idx_num]])
        save(list = paste0("RefElement_", idx), file = paste0(tmp_folder, "/RefElement_", idx, ".rda"))
    }
    save(gr_elements, file = paste0(tmp_folder, "/gr_elements.rda"))
    rm(RefElement_array)
    rm(gr_elements)
    if (score == "TRUE") 
    {   
        sorted_keys <- names(grouped_data)[order(as.numeric(gsub("chr", "", names(grouped_data))))]
        data_sorted <- grouped_data[sorted_keys]
        group_array <- lapply(data_sorted, identity)

        # generate data chunk
        lengths_array <- sapply(group_array, nrow)
        min_length_index <- which.min(lengths_array)
        min <- nrow(group_array[[min_length_index]])
        sub_blocks <- sapply(lengths_array, function(length) ceiling(length / min))
        sub_block_keys <- unlist(sapply(1:length(sub_blocks), function(i) paste0(i, ".", 1:sub_blocks[i])))
        
        rm(sorted_keys, data_sorted, group_array)

        chrs <- length(keys)
        ncpu = min(chrs, thread_num, parallel::detectCores())

        cl = parallel::makeCluster(ncpu)
        parallel::clusterExport(cl=cl, varlist=c("tmp_folder", "mutsFile"), envir=environment())
        registerDoParallel(cl)
        `%dopar2%` = foreach::`%dopar%`
        result = foreach::foreach(data = sub_block_keys) %dopar2% {
            library(data.table)

            split_data <- strsplit(data, "\\.")
            chr_num <- as.numeric(split_data[[1]][1])
            chr_name <- paste0("chr", chr_num)
            chunk <- as.numeric(split_data[[1]][2])

            db_file <- paste0(tmp_folder, "/database_", data, ".txt")
            split_file <- paste0(tmp_folder, "/data_", data, ".txt")

            maf_data <- fread(mutsFile,header = FALSE)
            grouped_data <- split(maf_data, maf_data[, 2])
            sorted_keys <- names(grouped_data)[order(as.numeric(gsub("chr", "", names(grouped_data))))]
            data_sorted <- grouped_data[sorted_keys]
            group_array <- lapply(data_sorted, identity)

            split_group <- lapply(group_array, function(data) {
            split_data <- split(data, 
                      rep(1:ceiling(nrow(data) / min), 
                      each = min, 
                      length.out = nrow(data)))
            })
            group <- copy(split_group[[chr_num]][[chunk]])

            rm(maf_data, grouped_data, sorted_keys, 
                data_sorted, group_array, split_group)
            
            group_file <- paste0(tmp_folder, "/chunk_", data, ".txt")
            write.table(group, file = group_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

            miss_file_path <- "/storage/yangjianLab/westlakechat/yanglabpipe_online/scripts/dNdSFun/log/missing.txt"
	    if (geneDB == "GRCh38"){
		    tbi_file <- paste0(score_database, "/whole_genome_SNVs.tsv.gz.", chr_name, ".gz.rankRawScore.GRCh38.gz")
	    }else{
                tbi_file <- paste0(score_database, "/whole_genome_SNVs.tsv.gz.", chr_name, ".gz.rankRawScore.gz")
	    }
            command <- paste0("cat ", group_file, " | ",
                            "awk -v var=", chr_name, " '{if($2==var)print $2\":\"$3\"-\"$3}' | ",
                            "sed 's/^chr//g'",
                            " > ", split_file)
            result <- system(command, intern = TRUE)

            command <- paste0("cat ", split_file, " | ",
                              "xargs tabix ", tbi_file, " | ",
                              "awk '{print $3,$4,$5}'",
                              " > ", db_file)
            result <- system(command, intern = TRUE)
            database <- fread(db_file, header = FALSE)

            maf <- data.frame()
            for (data_idx in 1:nrow(group)) 
            {   
                find <- FALSE
                data_row <- group[data_idx,]
                data_gene1 <- data_row[,4]
                data_gene2 <- data_row[,5]
                position <- data_row[, 3]

                if (!is.na(as.integer(position)))
                {
                    start_idx <- data_idx * 3 - 2
                    end_idx <- data_idx * 3
                    for (db_idx in start_idx:end_idx)
                    {
                        db_row <- database[db_idx, ]
                        db_gene1 <- db_row[, 1]
                        db_gene2 <- db_row[, 2]
                        if (data_gene1 == db_gene1 && data_gene2 == db_gene2)
                        {
                            input_result = c(data_row, db_row[, 3])
                            colnames(maf) <- colnames(input_result)
                            maf <- rbind(maf, input_result)
                            find <- TRUE
                            break
                        }
                    }
                }
                if (!find)
                {
                    row_string <- paste(data_row, collapse = " ")
                    message <- paste0(row_string, " not found.")
                    writeLines(message, miss_file_path)
                }
            }
            maf <- na.omit(maf)
            outfile <- paste0(tmp_folder, "/", data, ".txt")
            write.table(maf, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
            rm(group, database, maf)
        }
        stopCluster(cl)

        numbers <- sapply(strsplit(sub_block_keys, "\\."), function(x) as.numeric(x[1]))
        output <- data.frame()
        for (i in unique(numbers)) {
            files <- sub_block_keys[numbers == i]
            data <- do.call(rbind, lapply(files, function(file) read.table(paste0(tmp_folder, "/", file, ".txt"), header = FALSE)))
            output <- rbind(output, data)
            outfile <- paste0(tmp_folder, "/chr", i, "_all.txt")
            write.table(output, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
            output <- data.frame()
        }
    }
    rm(maf_data, grouped_data)

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
    chr_nums <- length(filt_chrs)
    ncpu = min(chr_nums, thread_num, parallel::detectCores())
    message(sprintf("\nAccording to the environment, %d threads will be created for calculation.\n", ncpu))
    cl = parallel::makeCluster(ncpu)
    parallel::clusterExport(cl=cl, varlist=c("nt", 
                                            "trinucMutsidx", "elements_list", 
                                            "cv", "max_muts_per_element_per_sample",
                                            "max_element_muts_per_sample", 
                                            "positiveThreshold", 
                                            "negativeThreshold",
                                            "mutsFile"), envir=environment())

    registerDoParallel(cl)
    `%dopar2%` = foreach::`%dopar%`
    result = foreach::foreach(data = filt_chrs) %dopar2% {
        library(data.table)
        library(GenomicRanges)
        library(GenomeInfoDb)

        if (score == "TRUE") {
            scoreFile <- paste0(tmp_folder, "/chr", data, "_all.txt")
            maf <- fread(scoreFile, header = FALSE)
        } else  {
            maf_data <- fread(mutsFile,header = FALSE)
            grouped_data <- split(maf_data, maf_data[, 2])
            rm(maf_data)
            maf <- copy(na.omit(grouped_data[[data]]))
            rm(grouped_data)
        }

        refDb_element <- paste0(tmp_folder, "/RefElement_", data, ".rda")
        ref_name <- load(refDb_element)
        load(paste0(tmp_folder, "/gr_elements.rda"))
        RefElement <- get(ref_name)
        rm(ref_name)
        
        elements_list=NULL
        
        # Step 1: Variables required
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
        maf <- maf[maf$V4 != maf$V5, ] # Remove mutations with identical reference and mutant base
        colnames(maf) <- c("sampleID","chr","pos","ref","alt","impScore")
        
        # Remove rows with NA values
        idxna <- which(is.na(maf),arr.ind=TRUE)
        if (nrow(idxna)>0) {
            maf <- maf[-unique(idxna[,1]),]
            warning(sprintf("%0.0f rows in the input table contained NA entries and have been removed. Please investigate.",length(unique(idxna[,1]))))
        }

        maf$chr <- ifelse(
            grepl("^chr", maf$chr),  # 检查是否以"chr"开头
            maf$chr,                 # 如果以"chr"开头，则保持不变
            paste0("chr", maf$chr)   # 否则，在原字符串前添加"chr"
        )

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
        if (nrow(ol) == 0) {
            stop("The overlap data frame is empty.")
        }

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
        maf$impidx <- as.integer(maf$impidx) 
        maf[maf$impScore<as.numeric(negativeThreshold),"impidx"] <- as.integer(1) # mutations group: neutral set
        ##mid set
        maf[maf$impScore>=as.numeric(positiveThreshold),"impidx"] <- as.integer(2) # mutations group: selection target set
        idxna <- which(is.na(maf),arr.ind=TRUE)
        # stop(maf$elementidx)
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

    unlink(tmp_folder, recursive = TRUE)

    maf_result <- lapply(result, function(x) x$maf)
    maf <- do.call(rbind, maf_result)
    RefElement1 <- unlist(lapply(result, function(x) x$RefElement1), recursive = FALSE)
    # RefElement <- unlist(RefElement_result, recursive = FALSE)
    exclsamples <- unlist(lapply(result, function(x) x$exclsamples))

    CADD_dndsWGSout <- dnds2wgs.noncoding(maf, RefElement1, exclsamples, negbeta, trinucMuts, outp)

    if(model=="1"){
        globaldnds <- CADD_dndsWGSout$globaldnds
        ####output of globaldnds
        info <- c(positive,negative,positiveThreshold,negativeThreshold,negbeta)
        globaldnds_res <- as.data.frame(matrix(c(unlist(globaldnds),info),nrow=1))
        colnames(globaldnds_res) <- c("mle","ci_low","ci_high","AIC","deviance",
                                    "overdis_chisq","overdis_ratio","overdis_p",
                                    "mle_qua","ci_low_qua","ci_high_qua",
                                    "positive","negative","positiveThreshold","negativeThreshold",
                                    "negbeta")
        write.table(globaldnds_res,globaldnds_outFile,sep="\t",row.names=F,quote=F)

    } else if (model=="2"){
        globaldnds <- CADD_dndsWGSout$globaldnds
        sel_cv <- CADD_dndsWGSout$sel_cv
        annotmuts <- CADD_dndsWGSout$annotmuts
        genemuts <- CADD_dndsWGSout$genemuts
        mle_submodel <- CADD_dndsWGSout$mle_submodel
        nbreg <- CADD_dndsWGSout$nbreg
        nbregind <- CADD_dndsWGSout$nbregind
        possmodel <- CADD_dndsWGSout$poissmodel

        ####output of globaldnds
        info <- c(positive,negative,positiveThreshold,negativeThreshold,negbeta)
        globaldnds_res <- as.data.frame(matrix(c(unlist(globaldnds),info),nrow=1))
        colnames(globaldnds_res) <- c("mle","ci_low","ci_high","AIC","deviance",
                                    "overdis_chisq","overdis_ratio","overdis_p",
                                    "mle_qua","ci_low_qua","ci_high_qua",
                                    "positive","negative","positiveThreshold","negativeThreshold",
                                    "negbeta")
        write.table(globaldnds_res,globaldnds_outFile,sep="\t",row.names=F,quote=F)

        ####output of sel_cv
        selcv_res <- as.data.frame(sel_cv)
        write.table(selcv_res,genelevel_selcv_outFile,sep="\t",row.names=F,quote=F)


    }    
    else{
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
        info <- c(positive,negative,positiveThreshold,negativeThreshold,negbeta)
        globaldnds_res <- as.data.frame(matrix(c(unlist(globaldnds),info),nrow=1))
        colnames(globaldnds_res) <- c("mle","ci_low","ci_high","AIC","deviance",
                                    "overdis_chisq","overdis_ratio","overdis_p",
                                    "mle_qua","ci_low_qua","ci_high_qua",
                                    "positive","negative","positiveThreshold",
                                    "negativeThreshold", "negbeta")
        write.table(globaldnds_res,globaldnds_outFile,sep="\t",row.names=F,quote=F)

        ####output of sel_loc
        selloc_res <- as.data.frame(sel_loc)
        write.table(selloc_res,"../OUT/test/dNdS_CADD.element.nb.out",sep="\t",row.names=F,quote=F)

        ####output of sel_cv
        selcv_res <- as.data.frame(sel_cv)
        write.table(selcv_res,genelevel_selcv_outFile,sep="\t",row.names=F,quote=F)
    }
    message(sprintf('Analysis is complete!'))
}
