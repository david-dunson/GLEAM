%PURPOSE: update the rho
%INPUT:
%S: I*J states
%S_0: I*23 initial states
%R: I*J recombination counts
%rho_0: I*1 prior proportion of population A
%nu_0: prior variance of rho_0
%chrom: J*1 chromosome number
%OUTPUT:
%rho: I*1 proportion of population A

function [rho]=rho_update(S, S_0, R, rho_0, nu_0,chrom)

[I,J]=size(R);

alpha_0=rho_0.^2.*(1-rho_0)./nu_0-rho_0;
beta_0=alpha_0.*(1./rho_0-1);

S_pre=zeros(I,J);
rho=zeros(I,1);

for k=1:23
    idx=(chrom==k); %for chromosome k;
    tmp=S(:,idx);
    S_pre(:,idx)=[S_0(:,k) tmp(:,1:(end-1)) ];
end

for i = 1:I
    idx=(R(i,:)~=0); %R_ij =/= 0
    idx_pre_0=(S_pre(i,idx)==0);
    idx_pre_1=(S_pre(i,idx)==1);
    idx_pre_2=(S_pre(i,idx)==2);
    idx_S_0=(S(i,idx)==0);
    idx_S_1=(S(i,idx)==1);
    idx_S_2=(S(i,idx)==2);
    idx_R_1=(R(i,idx)==1);
    idx_R_2=(R(i,idx)==2);
    idx_rho_a=(idx_pre_0&idx_S_1&idx_R_1)| ...
              (idx_pre_1&idx_S_2&idx_R_1)| ...
              (idx_pre_2&idx_S_2&idx_R_1)| ...
              (idx_S_1&idx_R_2);
    idx_rho_b=(idx_pre_0&idx_S_0&idx_R_1)| ...
              (idx_pre_1&idx_S_0&idx_R_1)| ...
              (idx_pre_2&idx_S_1&idx_R_1)| ...
              (idx_S_1&idx_R_2);
    
    rho(i)=randbeta(alpha_0(i) + sum(idx_rho_a)+2*sum(idx_S_2&idx_R_2),...   % for rho_i
                    beta_0(i)  + sum(idx_rho_b)+2*sum(idx_S_0&idx_R_2));   % for 1-rho_i
end



