%null models:
logBLR=zeros(100,J);
logBF=zeros(100,J);

for iter=1:100
%impute the proportion sampels:
idx=randsample(1296,J,true);

for j=1:J
   sample_prop(:,j)=randsample(Afr_prop(:,idx(j)),I,true);
end

sample_sub_prop=mean(sample_prop,2); %genome-wise average for each sub
sample_locus_prop=mean(sample_prop,1); % population averate for each locus
%%%%%%%%%%%%%%%%%%%%%%%%%%
%case1: no E;
E=zeros(I,1);
Y=binornd(1,0.5,I,1); %null model
casegroup=(Y==1);

%Admixture:
logBLR(iter,:)=logBLR_ADMIXTURE_GW(sample_prop, sample_sub_prop,casegroup);
end;
file_name=strcat(save_dir,'nullmodel_noE_adm.mat');
save(file_name,'logBLR','logBF');


for iter=1:100
%impute the proportion sampels:
idx=randsample(1296,J,true);

for j=1:J
   sample_prop(:,j)=randsample(Afr_prop(:,idx(j)),I,true);
end

sample_sub_prop=mean(sample_prop,2); %genome-wise average for each sub
sample_locus_prop=mean(sample_prop,1); % population averate for each locus

E=normrnd(0,1,I,1);
eta=0.5*E;
Y=binornd(1,exp(eta)./(1+exp(eta)));
casegroup=(Y==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%case2: alpha=0.5 
%Admixture:
logBLR(iter,:)=logBLR_ADMIXTURE_GW(sample_prop, sample_sub_prop,casegroup);
% %BF_QNM
% logBF(iter,:)=logBF_QNM_gprior_GW_binary(sample_prop,E,casegroup);
end
file_name=strcat(save_dir,'nullmodel_E0_5_adm.mat');
save(file_name,'logBLR','logBF');

for iter=1:100
%impute the proportion sampels:
idx=randsample(1296,J,true);

for j=1:J
   sample_prop(:,j)=randsample(Afr_prop(:,idx(j)),I,true);
end

sample_sub_prop=mean(sample_prop,2); %genome-wise average for each sub
sample_locus_prop=mean(sample_prop,1); % population averate for each locus

E=normrnd(0,1,I,1);
eta=1*E;
Y=binornd(1,exp(eta)./(1+exp(eta)));
casegroup=(Y==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%case3: alpha=1 
%Admixture:
logBLR(iter,:)=logBLR_ADMIXTURE_GW(sample_prop, sample_sub_prop,casegroup);
% %BF_QNM
% logBF(iter,:)=logBF_QNM_gprior_GW_binary(sample_prop,E,casegroup);
end
file_name=strcat(save_dir,'nullmodel_E1_adm.mat');
save(file_name,'logBLR','logBF');

for iter=1:100
%impute the proportion sampels:
idx=randsample(1296,J,true);

for j=1:J
   sample_prop(:,j)=randsample(Afr_prop(:,idx(j)),I,true);
end

sample_sub_prop=mean(sample_prop,2); %genome-wise average for each sub
sample_locus_prop=mean(sample_prop,1); % population averate for each locus

E=normrnd(0,1,I,1);
eta=2*E;
Y=binornd(1,exp(eta)./(1+exp(eta)));
casegroup=(Y==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%cas4: alpha=2
%Admixture:
logBLR(iter,:)=logBLR_ADMIXTURE_GW(sample_prop, sample_sub_prop,casegroup);
% %BF_QNM
% logBF(iter,:)=logBF_QNM_gprior_GW_binary(sample_prop,E,casegroup);
end
file_name=strcat(save_dir,'nullmodel_E2_adm.mat');
save(file_name,'logBLR','logBF');