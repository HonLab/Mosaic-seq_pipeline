#!/usr/bin/env Rscript
library(org.Hs.eg.db)

args <- commandArgs(TRUE)


df <- read.table(args[1], stringsAsFactors = F)

entry_ids <- df$V1
transcript_info <- select(org.Hs.eg.db, keys = entry_ids, 
                          columns = c("SYMBOL", "GENENAME"),
                          keytype = "REFSEQ")

output_dataframe <- cbind(transcript_info, df)
output_dataframe_name <- paste(sub(".txt", "", args[1]), 
                               '_wAnnotation.txt', sep = '')

write.table(output_dataframe[, -1], output_dataframe_name,
            quote = F, row.names = F, col.names = F, sep = '\t')


