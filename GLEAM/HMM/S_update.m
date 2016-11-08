%PURPOSE: update the S matrix
%INPUT
%SNPs: I*J SNPs without missing
%P: 3*3*J matrix observation probability mass matrix
%Q_cdn: 3*3*I*3 conditional state transition matrix
%Q_0: I*3  initial states mass matrix
%R: I*J recombination counts
%chrom: J*1 chromosome number
%OUTPUT:
%S: I*J state matrix
%S_0: I*23;
%logL: I*J log likelihood

function [S, S_0, L]=S_update(SNPs2, P, Q_cdn, Q_0, R, chrom)

[I, J]=size(R);

S=zeros(I,J);
S_0=zeros(I, 23);
L=zeros(I,J);

S_i=zeros(1,J);
L_i=zeros(1,J);

for i=1:I
    
    for k=1:23
    
        idx=(chrom==k); %for chromosome k;
        [S_ik,S_ik_0, L_ik]=FFBS2(P(:,:,idx),Q_cdn(:,:,i,:),Q_0(i,:)', R(i,idx),SNPs2(i,idx));
        
        S_i(1,idx)=S_ik;
        S_0(i,k)=S_ik_0;
        L_i(1,idx)=L_ik;
    
    end
    
    S(i,:)=S_i;
    L(i,:)=L_i;
end

