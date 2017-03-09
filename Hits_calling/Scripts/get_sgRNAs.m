function [unique_sgRNAs] = get_sgRNAs(experiment, region)

  prefix = '../step1.data/Normalized_data';

  % get all replicates
  dir_search_string = sprintf('%s/combined.%s.*',...
                              prefix,...
                              experiment);
  dir_files = dir(dir_search_string);

  all_sgRNAs = [];
  for rep_ind = 1:length(dir_files)
    ind = findstr(dir_files(rep_ind).name, '.rep');
    rep = dir_files(rep_ind).name(ind+1 : ind+1+3);
    
    % within each replicate, find sgRNAs
    search_string = sprintf('%s/combined.%s.%s/%s.*.cell_barcodes',...
                            prefix,...
                            experiment,...
                            rep,...
                            region);
    sgRNA_files = dir(search_string);

    for i = 1:length(sgRNA_files)
      ind = findstr(sgRNA_files(i).name, '.sgRNA');
      sgRNA = sgRNA_files(i).name(ind+1 : ind+1+5);

      all_sgRNAs{ length(all_sgRNAs) + 1 } = sgRNA;
    end
  end
  
  unique_sgRNAs = unique(all_sgRNAs);
end
