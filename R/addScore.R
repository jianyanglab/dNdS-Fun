addScore <- function(mutsFile, score_database = NULL, outFile)
{
    library(data.table)
    library(doParallel)

    if (is.null(score_database))
    {
      error_message <- "Please enter score_database file."
      stop(error_message)
    }

    tmp_folder <- paste0("score_tmp")
    if (!file.exists(tmp_folder)) {
        dir.create(tmp_folder)
    }

    # 读取，分类文件
    maf_data <- fread(mutsFile,header = FALSE)
    grouped_data <- split(maf_data, maf_data[, 2]) # chr
    chrs <- names(grouped_data)
    filt_chrs <- chrs[grepl("^(chr)?([1-9]|1[0-9]|2[0-2])$", chrs)]
    sorted_keys <- filt_chrs[order(as.numeric(filt_chrs))]
    data_sorted <- grouped_data[sorted_keys]
    group_array <- lapply(data_sorted, identity)

    # 数据以最小值分块，方便输入线程
    lengths_array <- sapply(group_array, nrow)
    min_length_index <- which.min(lengths_array)
    min <- nrow(group_array[[min_length_index]])
    sub_blocks <- sapply(lengths_array, function(length) ceiling(length / min))
    sub_block_names <- names(sub_blocks)
    sub_block_keys <- unlist(sapply(sub_block_names, 
                                function(i) paste0(gsub("chr", "", i), ".", 1:sub_blocks[i])))
    rm(data_sorted, group_array)

    # 开线程，拆分文件添加score
    chr_nums <- length(filt_chrs)
    ncpu = min(chr_nums, parallel::detectCores())
    cl = parallel::makeCluster(ncpu)
    parallel::clusterExport(cl=cl, varlist=c("tmp_folder", "mutsFile", "sorted_keys"), envir=environment())
    registerDoParallel(cl)
    `%dopar2%` = foreach::`%dopar%`
    result = foreach::foreach(data = sub_block_keys) %dopar2% 
    {
        library(data.table)

        split_data <- strsplit(data, "\\.")
        chr_num <- split_data[[1]][1]
        chr_name <- paste0("chr", chr_num)
        chunk <- split_data[[1]][2]

        db_file <- paste0(tmp_folder, "/database_", data, ".txt")
        split_file <- paste0(tmp_folder, "/data_", data, ".txt")

        maf_data <- fread(mutsFile,header = FALSE)
        grouped_data <- split(maf_data, maf_data[, 2])
        data_sorted <- grouped_data[sorted_keys]
        group_array <- lapply(data_sorted, identity)

        split_group <- lapply(group_array, function(data) 
        {
            split_data <- split(data, 
                  rep(1:ceiling(nrow(data) / min), 
                  each = min, 
                  length.out = nrow(data)))
        })
        if (grepl("chr", names(split_group)))
        {
            group <- copy(split_group[[chr_name]][[chunk]])
        } else 
        {
            group <- copy(split_group[[chr_num]][[chunk]])
        }
        rm(maf_data, grouped_data, data_sorted, group_array, split_group)
        
        group_file <- paste0(tmp_folder, "/chunk_", data, ".txt")
        write.table(group, file = group_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

        miss_file_path <- "missing.log"
        tbi_file <- paste0(score_database, "/whole_genome_SNVs.tsv.gz.", chr_name, ".gz.rankRawScore.gz")
        # 接受 chr1 vs 1
        command <- paste0("cat ", group_file, " | ",
                  "awk -v var=", chr_name, " '{if($2==var || \"chr\"$2==var)print $2\":\"$3\"-\"$3}' | ",
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
                file_conn <- file(miss_file_path, open = "a")
                writeLines(message, file_conn)
                close(file_conn)
            }
        }
        maf <- na.omit(maf)
        outfile <- paste0(tmp_folder, "/", data, ".txt")
        write.table(maf, file = outfile, sep = "\t", row.names = FALSE, col.names = FALSE)
        rm(group, database, maf)
    }
    stopCluster(cl)

    # 合并各染色体文件
    output <- data.frame()
    for (file in sub_block_keys) {
        file_path <- paste0(tmp_folder, "/", file, ".txt")
        # 检查文件是否存在并且非空
        if (file.exists(file_path) && !file.info(file_path)$size == 0) {
            data <- read.table(file_path, header = FALSE, colClasses = c("character", "character", "numeric", "character", "character", "numeric"))
            output <- rbind(output, data)
        }
    }
    outfile <- paste0(outFile, ".txt")
    write.table(output, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    rm(maf_data, grouped_data)
    unlink(tmp_folder, recursive = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
mutsFile <- args[1]
if (length(args) > 1)
{
  outFile <- args[2]
} else 
{
  outFile <- "output"
}
db <- "/storage/yangjianLab/sunxiwei/data/annotations/CADD/chr"
addScore(mutsFile, db, outFile)
