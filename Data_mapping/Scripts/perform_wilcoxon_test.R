#!/usr/bin/env Rscript

PerformWilcoxTest <- function(x, m, n) {
    a <- x[1 : m]
    b <- x[(m + 1) : (m + n)]
    result <- wilcox.test(a, b)
    return(result$p.value)
}

CalculateCV <- function(x) {
    return(sd(x) / mean(x))
}

args <- commandArgs(TRUE)
# output_file_name <- args[3]

expr_a <- read.table(args[1], row.names = 1, header = T, comment.char = "")
num_cells_a <- ncol(expr_a)
expr_b <- read.table(args[2], row.names = 1, header = T, comment.char = "")
num_cells_b <- ncol(expr_b)

combined_expr <- cbind(expr_a, expr_b)


pvalues <- apply(combined_expr, 1, PerformWilcoxTest, num_cells_a, num_cells_b)
adjusted_pvalues = p.adjust(pvalues, method = 'BH')

median_a <- apply(expr_a, 1, median)
median_b <- apply(expr_b, 1, median)
cv_a <-  apply(expr_a, 1, CalculateCV)
cv_b <-  apply(expr_b, 1, CalculateCV)

output_data <- data.frame('transcript' = row.names(combined_expr),
                          'pvalue' = pvalues, 'BH' = adjusted_pvalues,
                          'mean_a' = rowMeans(expr_a),
                          'mean_b' = rowMeans(expr_b),
                          'median_a' = median_a,
                          'median_b' = median_b,
                          'cv_a' = cv_a,
                          'cv_b' = cv_b)

write.table(output_data, , quote = F, sep = '\t', row.names = F)
