function [pval, median_nonsgRNA_expr, median_sgRNA_expr, log2ratio] ...
    = perform_DE(data_norm, cell_barcodes, sgRNA_cell_barcodes)

  index_sgRNA = get_cell_barcode_index(cell_barcodes, sgRNA_cell_barcodes);
  [num_genes, num_cells] = size(data_norm);

  % get merged datasets
  index_all_cells = 1:num_cells;
  index_nonsgRNA  = setdiff(index_all_cells, index_sgRNA);

  parfor gene_ind = 1:num_genes-1   % less one to remove sgRNA gene

    current_median = median(data_norm(gene_ind,:));
    low_cell_ind  = find(data_norm(gene_ind,:) <= current_median);
    overlap_count = length(intersect(cell_barcodes(low_cell_ind), sgRNA_cell_barcodes));
    current_pval  = hygecdf(overlap_count - 1,...
                            num_cells,...
                            length(sgRNA_cell_barcodes),...
                            length(low_cell_ind),...
                            'upper');

    pval(gene_ind)           = current_pval;

    median_sgRNA_expr(gene_ind)    = median( data_norm(gene_ind, index_sgRNA));
    median_nonsgRNA_expr(gene_ind) = median( data_norm(gene_ind, index_nonsgRNA));
    log2ratio(gene_ind)            = log2(median_sgRNA_expr(gene_ind) / median_nonsgRNA_expr(gene_ind));
  end
end
