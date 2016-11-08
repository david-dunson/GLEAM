function logBLR=logBLR_ADMIXTURE(local_prop, global_prop, r)

% Purpose:
% calcuate the log baeysian likelihood ration test statistics of ADMIXTURE.

% input:
% local_prop: I*3, risk propulation proportion of 0,1,2 at locus j
% global_prop: I*3, risk propulation proportion of 0,1,2 across the genome
% r: 1:k assumed risk model 

% ouput: 
% logBLR: 1*1 scaler of log10 likelihood ratio.

I=size(local_prop,1);
k=size(r,2);

lr_m=zeros(1,k);

for l=1:k
r_m=ones(I,3);
r_m(:,2:3)=repmat([r(l),r(l)^2],I,1);

lr_m(l)=sum(log10(sum(local_prop.*r_m,2)./...
                  sum(global_prop.*r_m,2)));

end

logBLR=mean(lr_m);