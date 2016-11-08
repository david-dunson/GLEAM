%PURPOSE: Initiate R matrix for recombination counts \in {0,1,2}
%INPUT: 
%gamma: 1*J recombination probability
%I: number of subject
%rand_seed: random seed
%OUTPUT: 
%R: I*J state matrix

function [R] = R_initiate(gamma, I, rand_seed)

RandState=rand_seed;
rand('state',RandState); % set arbitrary seed for uniform draws
randn('state',RandState); % set arbitrary seed for normal draws

J=length(gamma);

R=zeros(I,J);

for j=1:J
    R(:,j) = ...
        sample([(1-gamma(j))^2, 2*gamma(j)*(1-gamma(j)), gamma(j)^2]',I)-1;
end

