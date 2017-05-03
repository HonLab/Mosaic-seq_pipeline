function [prob_density] = get_pdf(x)
  % probability density
  prob_density = fitdist(x', 'kernel', 'Bandwidth', 0.2);

%  % evaluate probability density at lim
%  pdf_array  = pdf(pd, lim);
%  pdf_array  = pdf_array / sum(pdf_array);
end
