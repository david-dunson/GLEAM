%PURPOSE: update the R_ik for subject i and choromsome k
%INPUT
%S_ik: 1*J states
%S_0_ik: initial state
%Q_cdn_i: 3*3*1*3 conditional transition matrix for subject i
%gamma: 1_J recombination probabilities
%OUTPUT
%R_ik: 1*J recombination counts

%S_ik=S(i,idx); S_0_ik=S_0(i,k); Q_cdn_i=Q_cdn(:,:,i,:);gamma=gamma(idx);

function [R_ik]=R_ik_update(S_ik, S_0_ik, Q_cdn_i,gamma)

J=length(S_ik);
S_0ik=[S_0_ik,S_ik];
R_ik=zeros(1,J);

for j=1:J
   
    R_ik(j)=my_sample_1([Q_cdn_i(S_0ik(j)+1,S_0ik(j+1)+1,1,1)*(1-gamma(j))^2,...
                         Q_cdn_i(S_0ik(j)+1,S_0ik(j+1)+1,1,2)*2*gamma(j)*(1-gamma(j)),...
                         Q_cdn_i(S_0ik(j)+1,S_0ik(j+1)+1,1,3)*gamma(j)^2])-1;
   
end

