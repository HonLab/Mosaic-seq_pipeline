function [] = save_stats(pval,merged_nonsgRNA_expr, merged_sgRNA_expr, log2ratio, region, experiment, sgRNA)

  description = sprintf('%s.%s.%s', region, experiment, sgRNA);

  % write p-values for each sgRNA
  out_file = sprintf('pvals/pvals.%s.txt', description);
  dlmwrite(out_file, pval', '\t')
  out_file = sprintf('cpm/nonsgRNA_cpm.%s.txt', description);
  dlmwrite(out_file, merged_nonsgRNA_expr', '\t')
  out_file = sprintf('cpm/sgRNA_cpm.%s.txt', description);
  dlmwrite(out_file, merged_sgRNA_expr', '\t')
  out_file = sprintf('cpm/log2ratio.%s.txt', description);
  dlmwrite(out_file, log2ratio', '\t')
end
