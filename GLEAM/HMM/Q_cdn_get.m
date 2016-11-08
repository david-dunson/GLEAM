%PURPOSE: caculate Q matrix conditioned on recombination
%INPUT: 
%rho: I*1 population proportion
%OUTPUT: 
%Q: 3*3*I*3 matrix state transition matrix

function [Q_cdn] = Q_cdn_get(rho)

I=length(rho);

Q_cdn=zeros(3,3,I,3);

for i=1:I
one_rho=1-rho(i);

Q_cdn(:,:,i,1)=diag([1,1,1]);
Q_cdn(:,:,i,2)=[one_rho,rho(i),0;one_rho/2,1/2,rho(i)/2;0,one_rho,rho(i)];
Q_cdn(:,:,i,3)=repmat([one_rho^2 2*rho(i)*one_rho rho(i)^2],3,1);

end

