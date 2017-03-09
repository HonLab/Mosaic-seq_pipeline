#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)


args <- commandArgs(trailingOnly = TRUE)

df <- read.table(args[1])
df <- subset(df, V12 != 0)

for (i in 3:11) {
    df[[i]] <- (df[[i]] - df[, 12]) / df[, 12]
}


colnames(df) <- c('transcript', 'length', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%')


df <- df[, -12]
dfm <- melt(df, id.vars = c(1: 2))

pdf(sub('txt', 'pdf', args[1]))
ggplot(dfm, aes(variable, value)) + geom_boxplot() + labs(x = '% of total reads', y = 'Percent relative error')
dev.off

sink(sub('txt', 'log', args[1]), append = FALSE, split = FALSE)
print(colMeans(abs(df[, c(3:11)])))






