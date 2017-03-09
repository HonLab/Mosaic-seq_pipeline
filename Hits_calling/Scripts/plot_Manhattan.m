function [] = plot_Manhattan(pvals, description)

  SMYD3_ind = [5051 5052];
  CD52_ind  = 828;
  PIM1_ind  = [16284 16285];
  HBG2_ind  = [26940 26941];
  HBE1_ind  = 26942;
  FADS1_ind = 27995;
  UFD1L_ind = [47727 47728];
  PRKAR2B_ind = 19849;

  scatter([1:length(pvals)],-log10(pvals), 10,...
          'o',...
          'MarkerFaceColor','k',...
          'MarkerEdgeColor','k',...
          'MarkerFaceAlpha', 0.1,...
          'MarkerEdgeAlpha', 0.1);

  ymax = max([10 max(-log10(pvals)) + 1]);
  xlim([-1000 length(pvals)+1000]);
  ylim([0 ymax]);
  hold on;

  scatter(PIM1_ind, -log10(pvals(PIM1_ind)), 30,...
          'MarkerFaceColor',[1 0 0],...
          'MarkerEdgeAlpha', 0);

  scatter(SMYD3_ind,-log10(pvals(SMYD3_ind)), 30,...
          'MarkerFaceColor',[0 1 0],...
          'MarkerEdgeAlpha', 0);
  
  scatter(HBG2_ind, -log10(pvals(HBG2_ind)), 30,...
          'MarkerFaceColor',[0 0 1],...
          'MarkerEdgeAlpha', 0);
  
  scatter(HBE1_ind, -log10(pvals(HBE1_ind)), 30,...
          'MarkerFaceColor',[0.5 0 1],...
          'MarkerEdgeAlpha', 0);
  
  scatter(CD52_ind, -log10(pvals(CD52_ind)), 30,...
          'MarkerFaceColor',[1 1 0],...
          'MarkerEdgeAlpha', 0);
  
  scatter(FADS1_ind,-log10(pvals(FADS1_ind)), 30,...
          'MarkerFaceColor',[1 0.5 1],...
          'MarkerEdgeAlpha', 0);
  
  scatter(UFD1L_ind,-log10(pvals(UFD1L_ind)), 30,...
          'MarkerFaceColor',[1 0.5 0],...
          'MarkerEdgeAlpha', 0);
  
  scatter(PRKAR2B_ind,-log10(pvals(PRKAR2B_ind)), 30,...
          'MarkerFaceColor',[0 0.5 0.8],...
          'MarkerEdgeAlpha', 0);
  
  legend('Others','PIM1','SMYD3','HBG2','HBE1', 'CD52','FADS1','UFD1L','PRKAR2B',...
         'Location', 'eastoutside');
  hold off;

  out_file = sprintf('plots/Manhattan_plot.%s.p_val.combined.tiff', description);
  %print (gcf, out_file, '-dpdf');
  print (gcf, out_file, '-r600', '-dtiff');
end