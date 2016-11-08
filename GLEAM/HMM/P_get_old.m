%PURPOSE:
%calcualte P matrix
%INPUT:
%pop_a_var:J by 1 population A variant allele frequncies
%pop_b_var: J by 1 population B variant allele frequncies
%OUTPUT: 
%P 3*3*J matrix observation probability mass matrix
function [P] = P_get_old(pop_a_var,pop_b_var)

J=length(pop_a_var);

P=zeros(3,3,J);
j=1;
for j=1:J
    pa=pop_a_var(j);
    pb=pop_b_var(j);
    P(1,1,j)=(1-pb)^2;
    P(1,2,j)=2*pb*(1-pb);
    P(1,3,j)=pb^2;
    P(2,1,j)=(1-pa)*(1-pb);
    P(2,2,j)=pa*(1-pb)+pb*(1-pa);
    P(2,3,j)=pa*pb;
    P(3,1,j)=(1-pa)^2;
    P(3,2,j)=2*pa*(1-pa);
    P(3,3,j)=pa^2;
end