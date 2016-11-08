%PURPOSE: update the gamma
%INPUT:
%gamma_0: 1*J the prior of gamma based on poisson process
%mu_0: prior variance
%R: I*J recombination counts
%OUTPUT:
%gamma: 1*J recombination probability

function [gamma]=gamma_update(gamma_0, mu_0, R)

[I,J]=size(R);

alpha_0=gamma_0.^2.*(1-gamma_0)./mu_0-gamma_0;
beta_0=alpha_0.*(1./gamma_0-1);

R_colsum=sum(R,1);

gamma=randbeta(alpha_0+R_colsum, beta_0+2*I-R_colsum);




    