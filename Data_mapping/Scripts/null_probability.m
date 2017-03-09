function [p] = null_probability(x, B)
  B_0 = sum(x > 0);
  p = (1/B)^sum(x) * factorial(B) / factorial(B - B_0);
end
