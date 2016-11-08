%PURPOSE: impute the missing SNPs
%INPUT: 
%SNPs: I*J matrix of SNPs with missing
%P 3*3*J matrix observation probability mass matrix
%S: I*J state matrix 
%SNPs_missing_idx: I*J missing index matrix 
%OUTPUT:
%SNPS: I*J matrix of SNPs without missing

function SNPs_missing_impute(obj,S,P,SNPs_missing_idx)
    
[I,J]=size(S);

j=1;
for j=1:J
    idx=SNPs_missing_idx(:,j);
    S_missing=S(idx,j);
    P_j=P(:,:,j);
    P_j_missing=zeros(length(S_missing),1);
    for k=1:length(S_missing)
        %P_j_missing(k)=sample(P_j(S_missing(k)+1,:),1);
        P_j_missing(k)=my_sample_1(P_j(S_missing(k)+1,:));
    end
    obj.SNPs(idx,j)=P_j_missing-1;
end

end 
 
