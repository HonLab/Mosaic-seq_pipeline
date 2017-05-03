function [simulated] = reduce_expression(rand_sub, ...
                                         penetrance,...
                                         repression)

  % choose cells to change expression
  num_cells = length(rand_sub);
  cells_chosen = find(rand(num_cells, 1) <= penetrance);
  cells_remain = setdiff(1:num_cells, cells_chosen);

  % determine amount of expression to change
  cells_chosen_expr = rand_sub(cells_chosen) * (1-repression);
  cells_remain_expr = rand_sub(cells_remain);
  simulated = [cells_remain_expr cells_chosen_expr];
end
