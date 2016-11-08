%PURPOSE: caculate Q_0 matrix
%INPUT: 
%rho: I*1 population proportion
%OUTPUT: 
%Q_0 I*3 matrix of initial states mass matrix 

function [Q_0] = Q_0_get_old(rho)

I=length(rho);

Q_0=[];

for i=1:I
    Q_0_i=[(1-rho(i))^2, 2*rho(i)*(1-rho(i)),rho(i)^2];
    Q_0=[Q_0;Q_0_i];
end
