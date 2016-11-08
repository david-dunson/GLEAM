%PURPOSE: caculate Q_0 matrix
%INPUT: 
%rho: I*1 population proportion
%OUTPUT: 
%Q_0 I*3 matrix of initial states mass matrix 

function [Q_0] = Q_0_get(rho)

I=length(rho);

Q_0=zeros(I,3);

Q_0(:,1)=(1-rho).^2;
Q_0(:,2)=2.*rho.*(1-rho);
Q_0(:,3)=rho.^2;

