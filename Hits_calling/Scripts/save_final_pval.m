function [] = save_final_pval(pval_final, region)
  out_file = sprintf('pvals/pvals.%s.combined_pval.txt', region);
  dlmwrite(out_file, pval_final, '\n')
end
