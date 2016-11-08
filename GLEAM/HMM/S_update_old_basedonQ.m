%PURPOSE: update the S matrix
%INPUT
%SNPs: I*J SNPs without missing
%P 3*3*J matrix observation probability mass matrix
%Q: 3*3*I*J matrix state transition matrix
%Q_0: I*3 matrix of initial states mass matrix
%chrom: J*1 chromosome number
%OUTPUT:
%S: I*J state matrix
%logL: I*J log likelihood
%Q_B_marg: 1*J*I posterior marginal densities

%function [S, logL, Q_B_marg]=S_update(SNPs, P, Q, Q_0,chrom)

function [S, logL]=S_update(SNPs, P, Q, Q_0,chrom)

[tmp1,tmp2, I, J]=size(Q);

S=zeros(I,J);
logL=zeros(I,J);
%Q_B_marg=zeros(3, J, I);
S_i=zeros(1,J);
logL_i=zeros(1,J);

for i=1:I
    %S_i=[];
    %logL_i=[];
    %Q_B_marg_i=[];
    for k=1:23
        idx=(chrom==k); %for chromosome k;
        %[S_ik, logL_ik, Q_B_marg_ik]=FFBS(P(:,:,idx),Q(:,:,i,idx),Q_0(i,:)',SNPs(i,idx));
        [S_ik, logL_ik]=FFBS2(P(:,:,idx),Q(:,:,i,idx),Q_0(i,:)',SNPs(i,idx));
        
        S_i(1,idx)=S_ik;
        logL_i(1,idx)=logL_ik;
        %S_i=[S_i S_ik];
        %logL_i=[logL_i logL_ik];
        %Q_B_marg_i=[Q_B_marg_i Q_B_marg_ik];
    end
    S(i,:)=S_i;
    logL(i,:)=logL_i;
    %Q_B_marg(:,:,i)=Q_B_marg_i;
end

