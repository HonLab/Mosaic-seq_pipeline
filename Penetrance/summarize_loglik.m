function [] = summarize_loglik(all_fraction_penetrance, ...
                               all_fraction_repression,...
                               all_loglik,...
                               description)

  % get ave
  for penetrance_ind = 1:length(all_fraction_penetrance)
    penetrance = all_fraction_penetrance(penetrance_ind);

    % iterate through change of expression
    for repression_ind = 1:length(all_fraction_repression)
      repression = all_fraction_repression(repression_ind);

      median_loglik(penetrance_ind, repression_ind)...
          = median(all_loglik(:, penetrance_ind, repression_ind));
    end
  end
  
  median_loglik

  % make heatmap
  figure
  max_y = max( [max(max(median_loglik)) 5] );
  imagesc(median_loglik, [0 max_y]);
  %imagesc(median_loglik, [-max_y max_y]);
  set_colormap_red;
  set(gca, 'XTick',      1:length(all_fraction_repression));
  set(gca, 'XTickLabel', all_fraction_repression * 100);
  set(gca, 'YTick',      1:length(all_fraction_penetrance));
  set(gca, 'YTickLabel', all_fraction_penetrance * 100);
  xlabel('%repression');
  ylabel('%penetrance');
  title(description);
  colorbar();
  
  % plot max point
  [max_val, max_ind] = max(median_loglik(:));
  [ind_penetrance, ind_repression] = ind2sub(size(median_loglik), max_ind);
  hold on;
  plot(ind_repression, ind_penetrance, 'ko');
  text(ind_repression + 0.5, ind_penetrance + 0.5, ...
       sprintf('%.1f', median_loglik(ind_penetrance, ind_repression)));
  hold off;
  
  print('-dpdf', sprintf('plots/loglik.%s.pdf',...
                         description));

  % save stats
  out_file = sprintf('stats/median_loglik.%s.txt', description);
  dlmwrite(out_file, median_loglik, '\t')

  % save all iterations
  [num_iter, num_penetrance, num_repression] = size(all_loglik);
  for iter = 1:num_iter
    out_file = sprintf('all_stats/loglik.%s.iter=%d.txt', description, iter);
    out = reshape(all_loglik(iter,:,:), num_penetrance, num_repression);
    dlmwrite(out_file, out, '\t');
  end
end
