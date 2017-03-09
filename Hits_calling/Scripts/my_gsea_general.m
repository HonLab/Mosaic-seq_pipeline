clear;

setup_and_load_constants;

experiments = {'Batch_A', 'Batch_B', 'Batch_C'};

% read expression data
[data_norm, cell_barcodes] = load_all_expr(experiments);

for region_ind = 1:length(regions)
  region = regions{region_ind};

  pval_array = [];
  for experiment_ind = 1:length(experiments)
    experiment = experiments{experiment_ind};
    [num_genes, num_cells] = size(data_norm{experiment_ind});

    % get list of sgRNAs
    all_sgRNAs = get_sgRNAs(experiment, region);

    for sgRNA_ind = 1:length(all_sgRNAs)
      sgRNA = all_sgRNAs{sgRNA_ind};
      fprintf('region = %s, experiment = %s, sgRNA = %s\n', region, experiment, sgRNA);

      % load cell barcodes for the given sgRNA
      sgRNA_cell_barcodes = load_sgRNA_cell_barcodes(experiment, region, sgRNA);
      if (length(sgRNA_cell_barcodes) == 0)
        continue;
      end

      % perform differential expression over all genes
      [pval, merged_nonsgRNA_expr, merged_sgRNA_expr, log2ratio] =...
          perform_DE(data_norm{experiment_ind},...
                     cell_barcodes{experiment_ind},...
                     sgRNA_cell_barcodes);

      % save p-values
      next_index = size(pval_array,2) + 1;
      pval_array(:,next_index) = pval;

      % save output
      save_stats(pval,merged_nonsgRNA_expr, merged_sgRNA_expr, log2ratio,...
                 region, experiment, sgRNA);
      
      % Manhattan plot on sgRNA alone
      plot_Manhattan(pval, sprintf('%s.%s.%s', region, experiment, sgRNA));
    end
  end

  % combine pvals, save, and plot
  pval_sgRNAs_combined = combine_sgRNA_pvals(pval_array);
  save_final_pval(pval_sgRNAs_combined, region);
  plot_Manhattan(pval_sgRNAs_combined, region);
end
