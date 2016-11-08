%PURPOSE: caculate Q matrix
%INPUT: 
%gamma: J*1 probability of recombination for ONE chromosome
%rho: I*1 population proportion
%OUTPUT: 
%Q: 3*3*I*J matrix state transition matrix


function [Q] = Q_get(gamma, rho)

I=length(rho);
J=length(gamma);

Q=zeros(3,3,I,J);

for i=1:I
    for j=1:J
        ab=gamma(j)*(1-rho(i));
        aa=1-ab;
        ba=gamma(j)*rho(i);
        bb=1-ba;
        Q(1,1,i,j)= bb^2;
        Q(1,2,i,j)=2*ba*bb;
        Q(1,3,i,j)=ba^2;
        Q(2,1,i,j)=ab*bb;
        Q(2,2,i,j)=aa*bb+ba*ab;
        Q(2,3,i,j)=ba*aa;
        Q(3,1,i,j)=ab^2;
        Q(3,2,i,j)=2*ab*aa;
        Q(3,3,i,j)=aa^2;
    end
end
