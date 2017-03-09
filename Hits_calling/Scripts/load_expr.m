function [all_data, all_cell_barcodes] = load_expr(experiment)
  prefix = '../step1.data/Normalized_data';

  % get list of all replicates
  search_string = sprintf('%s/combined.%s.rep*',...
                          prefix,...
                          experiment);
  rep_files = dir(search_string);

  all_rep = [];
  for i = 1:length(rep_files)
    ind = findstr(rep_files(i).name, '.rep');
    all_rep{i} = rep_files(i).name(ind+1 : ind+1+3);
  end

  % load all the data and barcodes
  all_data = [];
  all_cell_barcodes = [];
  for i = 1:length(all_rep)
    full_experiment = sprintf('%s.%s', experiment, all_rep{i});
    fprintf('loading experiment: %s\n', full_experiment);

    % data
    rep_data = load(sprintf('%s/combined.%s/normalized.all.matrix.txt.no_header',...
                            prefix,...
                            full_experiment));
    all_data = [all_data rep_data];
    fprintf('  loaded %d transcriptomes\n', size(rep_data, 2));
  
    % barcodes
    rep_cell_barcodes = load_cell_barcodes(full_experiment);
    all_cell_barcodes = [all_cell_barcodes ; rep_cell_barcodes];
    fprintf('  loaded %d cell barcodes\n', length(rep_cell_barcodes));
  
    if (size(all_data,2) ~= length(all_cell_barcodes))
      error ('Error. \nThe cell barcode number does not match the matrix size')
    end
  end
  
  fprintf('Done loading %d transcriptomes\n', size(all_data, 2));
end
