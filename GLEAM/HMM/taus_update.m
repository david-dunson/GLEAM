% PURPOSE
% updata tau_a_old tau_b_old
% INPUT:
% pa: 1*J current population A's proportion of variant alleles
% pb: 1*J current population B's proportion of variant alleles
% sigma: SD for normal proposal
% OUTPUT:
% tau_a, tau_b
% tau_a_old=300; tau_b_old=300; sigma=10;
% tau_a=300; tau_b=300; sigma=10;
function [tau_a, tau_b]=taus_update(pa, pb,pa_0, pb_0, tau_a_old, tau_b_old,sigma, range)

%range=[50,1000];
tau_a_prop = tau_a_old + TruncatedGaussian(sigma, range-tau_a_old,1); 


%tau_a_prop=tau_a_old+randn(1)*sigma;




logR=sum(log(betapdf(pa,tau_a_prop*pa_0,tau_a_prop*(1-pa_0))))-...
     sum(log(betapdf(pa,tau_a_old*pa_0, tau_a_old*(1-pa_0))));

u=rand(1);

if log(u)<min(0,logR)
    tau_a=tau_a_prop;
else
    tau_a=tau_a_old;
end

tau_b_prop = tau_b_old + TruncatedGaussian(sigma, range-tau_b_old,1); 


%tau_b_prop=tau_b_old+randn(1)*sigma;

logR=sum(log(betapdf(pb,tau_b_prop*pb_0,tau_b_prop*(1-pb_0))))-...
     sum(log(betapdf(pb,tau_b_old*pb_0,tau_b_old*(1-pb_0))));

u=rand(1);

if log(u)<min(0,logR)
    tau_b=tau_b_prop;
else
    tau_b=tau_b_old;
end




