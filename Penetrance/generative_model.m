function [] = generative_model(cell_barcodes, data_norm, ...
                               sgRNA_cell_barcodes, gene_ind, ...
                               num_rand, description)

  [num_genes, num_cells] = size(data_norm);

  all_fraction_penetrance = 0.05 : 0.05 : 0.95;
  all_fraction_repression = 0.05 : 0.05 : 0.95;

  % expression
  norm_expr = data_norm(gene_ind, :);

  % get expression of sgRNA cells
  sgRNA_cell_index = get_cell_barcode_index(cell_barcodes, sgRNA_cell_barcodes);
  sgRNA_cell_expr  = norm_expr ( sgRNA_cell_index ) ;
  sgRNA_pdf        = get_pdf(sgRNA_cell_expr);

  % Get pdf from all data. This is the null.
  %null_pdf = get_pdf(norm_expr_sorted);
  null_pdf = get_pdf(norm_expr);

  % evaluate null at sgRNA cell expression
  null_prob = pdf(null_pdf, sgRNA_cell_expr);

  all_loglik = zeros(num_rand, length(all_fraction_penetrance), length(all_fraction_repression));

  % iterate through penetrance
  for penetrance_ind = 1:length(all_fraction_penetrance)
    penetrance = all_fraction_penetrance(penetrance_ind);

    % iterate through change of expression
    for repression_ind = 1:length(all_fraction_repression)
      repression = all_fraction_repression(repression_ind);
      fprintf('penetrance = %.2f, repression = %.2f\n', penetrance, repression);

      %sub_sim_pdf = {};
      sub_log_lik = [];
      parfor rand_ind = 1:num_rand
        %fprintf('rand_ind = %d\n', rand_ind);

        % create a new expression distribution using rand_sub:
        % 1. choose "penetrance" percent of cells
        % 2. reduce expression using "repression"
        simulated = reduce_expression(norm_expr, penetrance, repression);

        % create a pdf from simulated data. Evaluate at sgRNA expressed cell expression
        sim_pdf  = get_pdf(simulated);
        sim_prob = pdf(sim_pdf, sgRNA_cell_expr);

        % compare pdf to null using LRT
        log_lik = sum(log10(sim_prob ./ null_prob));
      
        sub_log_lik(rand_ind) = log_lik;
        %sub_sim_pdf{rand_ind} = sim_pdf;
      end
    
      %print_pdf_dist(sub_sim_pdf, null_pdf, sgRNA_pdf, description, ...
      %               penetrance, repression);
      all_loglik(:, penetrance_ind, repression_ind) = sub_log_lik;
    end
  end

  summarize_loglik(all_fraction_penetrance,...
                   all_fraction_repression,...
                   all_loglik,...
                   description);
end
