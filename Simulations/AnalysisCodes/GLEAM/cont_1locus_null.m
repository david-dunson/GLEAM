% Purpose:
% run single locus for continuous  outcomes in Cluster
 
idnum_str = getenv('SGE_TASK_ID');
idnum = str2num(idnum_str);
 
%outcome=getenv('JOB_NAME');
addpath(genpath('/home/chg/bz27/MATLAB'))
data_dir='/home/chg/bz27/Admixture/Datasets/Sdata/Real/';
save_dir=strcat('/home/chg/bz27/Admixture/Results/');
 
RandState=4000;
c=0;
rand('state',RandState+idnum*100+c*100000); % set arbitrary seed for uniform draws
randn('state',RandState+idnum*100+c*100000);
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=1000;
J=10000;
 
file_name=strcat(data_dir,'Afr_prop.mat');
load(file_name);
 
locus_prop=mean(Afr_prop,1);
 
sub_prop=mean(Afr_prop,2);
 
top20_idx=locus_prop>quantile(locus_prop,0.8);
 
sample_prop=zeros(I,J);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%impute the proportion sampels:
idx=randsample(find(~top20_idx),J,true);

for j=1:J
   sample_prop(:,j)=randsample(Afr_prop(:,idx(j)),I,true);
end

sample_sub_prop=mean(sample_prop,2); %genome-wise average for each sub
sample_locus_prop=mean(sample_prop,1); % population averate for each locus

%E=normrnd(0,1,I,1);
E=zeros(I,1);
eta=c*E;
Y=normrnd(eta,1,I,1); % GAM use continuous.

%Admixture:
logBLR=logBLR_ADMIXTURE_GW(sample_prop, sample_sub_prop,Y>quantile(Y,0.8));
%BF_QNM
logBF=logBF_QNM_gprior_GW(sample_prop,E,Y,'continuous');
file_name=strcat(save_dir,'null_E0_',idnum_str,'_cont.mat');
save(file_name,'Y', 'logBLR', 'logBF');
 
