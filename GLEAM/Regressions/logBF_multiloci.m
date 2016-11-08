function [max_idx,max_logBF]=logBF_multiloci(logBF,Y,S,response, cutoff)

nonnull_idx=logBF>cutoff;

I=size(S,1);
K=sum(nonnull_idx);


E=zeros(I,1);
multiloci_logBF=[];

S_tmp=S(:,nonnull_idx);

for j=1:min([K,4])
comb=nchoosek(1:K,j);
L=size(comb,1);
logBF_tmp=zeros(L,1);
for i=1:L
[x,fval]=fminbnd(@(x)-logBF_QNM_gprior(Y,S_tmp(:,comb(i,:)),E,x,response),0.0001,0.2);
logBF_tmp(i)=-fval;    
end
multiloci_logBF=[multiloci_logBF;[ones(L,1)*j,(1:L)',logBF_tmp]];
end

size(multiloci_logBF);
[C,I]=max(multiloci_logBF(:,3));
tmp=multiloci_logBF(I,:);
comb_tmp=nchoosek(1:K,tmp(1));

tmp2=find(nonnull_idx);

max_logBF=tmp(3);
max_idx=tmp2(comb_tmp(tmp(2),:));