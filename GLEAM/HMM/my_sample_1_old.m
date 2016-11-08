function x = my_sample_1(p)
%my_sample_1    Sample from categorical distribution.

  cdf = cumsum(p(:));
  x = sum(cdf < rand*cdf(end)) + 1;

