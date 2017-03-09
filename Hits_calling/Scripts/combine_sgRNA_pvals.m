function [pval_combined] = combine_sgRNA_pvals(pval_array)
  [num_genes, num_sgRNAs] = size(pval_array);
  p_log = -2*log(pval_array);
  p_sum = sum(p_log, 2);
  p_chi2 = chi2cdf(p_sum, 2*num_sgRNAs, 'upper');
  %p_fdr = mafdr(p_chi2, 'BHFDR', 'true');
  pval_combined = p_chi2;
end
