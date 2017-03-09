function [all_expr, all_cell_barcodes] = load_all_expr(experiments)
  all_expr = [];
  all_cell_barcodes = [];
  num_cells_total = 0;
  for i = 1:length(experiments)
    [all_expr{i}, all_cell_barcodes{i}] = load_expr(experiments{i});
    num_cells_total = num_cells_total + length(all_cell_barcodes{i});
  end
  fprintf('Final dataset: %d cells\n\n', num_cells_total);
end
