%PURPOSE: simulate the latent states
%INPUT: 
%Q: 3*3*J transition matrix
%Q_0: initial probability mass probability
%P: 3*3*J observation mass matrix
%rand_seed: random seed
%OUTPUT:
%S: 1*J states

function [X,S]=state_simu(Q, Q_0,P,rand_seed)


RandState=rand_seed;
rand('state',RandState); % set arbitrary seed for uniform draws
randn('state',RandState); % set arbitrary seed for normal draws


[tmp,tmp,J]=size(Q);

S=zeros(1,J);
X=zeros(1,J);

S_prev=randsample(3,1,true,Q_0);
for j=1:J
S_prev=randsample(3,1,true,Q(S_prev,:,j));
X(j)=randsample(3,1,true,P(S_prev,:,j))-1;
S(j)=S_prev-1;
end
