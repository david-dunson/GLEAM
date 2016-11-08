%PURPOSE: update the R
%INPUT
%S: I*J states
%S_0: I*23 initial state
%Q_cdn: 3*3*I*3 conditional transition matrix for subject i
%gamma: 1*J recombination probability
%chrom: J*1 chromosome number
%OUTPUT
%R: I*J recombination counts

function [R]=R_update(S, S_0, Q_cdn,gamma,chrom)

[I, J]=size(S);

R=zeros(I,J);
R_i=zeros(1,J);

for i=1:I
    for k=1:23
         idx=(chrom==k); %for chromosome k;
         R_ik=R_ik_update(S(i,idx), S_0(i,k), Q_cdn(:,:,i,:),gamma(idx));
         R_i(1,idx)=R_ik;
    end
    R(i,:)=R_i;
end
