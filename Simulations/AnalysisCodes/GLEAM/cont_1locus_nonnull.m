% run single locus for continuous outcomes in Cluster
 
idnum_str = getenv('SGE_TASK_ID');
idnum = str2num(idnum_str);
 
%outcome=getenv('JOB_NAME');
addpath(genpath('/home/chg/bz27/MATLAB'))
data_dir='/home/chg/bz27/Admixture/Datasets/Sdata/Real/';
save_dir=strcat('/home/chg/bz27/Admixture/Results/');
 
RandState=4000;
c=1;
rand('state',RandState+idnum*100+c*100000); % set arbitrary seed for uniform draws
randn('state',RandState+idnum*100+c*100000);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=1000;
J=100;
 
file_name=strcat(data_dir,'Afr_prop.mat');
load(file_name);
 
locus_prop=mean(Afr_prop,1);
 
sub_prop=mean(Afr_prop,2);
 
top20_idx=locus_prop>quantile(locus_prop,0.8);
 
sample_prop=zeros(I,J);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%single loci nonull models: continuous NO E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
response='continuous';
logBLR=zeros(1,5*J);
logBF=zeros(1,5*J);


%impute the nonnull models:
idx=randsample(find(top20_idx),J,true);

for j=1:J
   sample_prop(:,j)=randsample(Afr_prop(:,idx(j)),I,true);
end

sample_locus_prop=mean(sample_prop,1); % population averate for each locus

S=zeros(I,J);
for j=1:J
S(:,j)=mnrnd(1,[(1-sample_prop(:,j)).^2,2.*sample_prop(:,j).*(1-sample_prop(:,j)),sample_prop(:,j).^2])*[0;1;2];
end;
S=S-repmat(mean(S),I,1);

%E=zeros(I,J);
E=normrnd(0,1,I,J);
%%%%%%%%%%%%%%%%%%%%%%
%1
eta=0.2.*repmat(sample_locus_prop,I,1).*S+c.*E;
Y=normrnd(eta,1,I,J);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Admixture:
logBLR(1,1:J)=logBLR_ADMIXTURE_GW_nonull(sample_prop,unifrnd(0.81,0.84,I,1),Y>repmat(quantile(Y,0.8),I,1));
%BF_QNM
logBF(1,1:J)=logBF_QNM_gprior_GW_nonull(S,E,Y,response);

%%%%%%%%%%%%%%%%%%%%%%
%2
eta=0.25.*repmat(sample_locus_prop,I,1).*S+c.*E;
Y=normrnd(eta,1,I,J);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Admixture:
logBLR(1,1+J:2*J)=logBLR_ADMIXTURE_GW_nonull(sample_prop,unifrnd(0.81,0.84,I,1),Y>repmat(quantile(Y,0.8),I,1));
%BF_QNM
logBF(1,1+J:2*J)=logBF_QNM_gprior_GW_nonull(S,E,Y,response);

%%%%%%%%%%%%%%%%%%%%%%
%3
eta=0.3.*repmat(sample_locus_prop,I,1).*S+c.*E;
Y=normrnd(eta,1,I,J);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Admixture:
logBLR(1,1+2*J:3*J)=logBLR_ADMIXTURE_GW_nonull(sample_prop,unifrnd(0.81,0.84,I,1),Y>repmat(quantile(Y,0.8),I,1));
%BF_QNM
logBF(1,1+2*J:3*J)=logBF_QNM_gprior_GW_nonull(S,E,Y,response);

%%%%%%%%%%%%%%%%%%%%%%
%4
eta=0.35.*repmat(sample_locus_prop,I,1).*S+c.*E;
Y=normrnd(eta,1,I,J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Admixture:
logBLR(1,1+3*J:4*J)=logBLR_ADMIXTURE_GW_nonull(sample_prop,unifrnd(0.81,0.84,I,1),Y>repmat(quantile(Y,0.8),I,1));
%BF_QNM
logBF(1,1+3*J:4*J)=logBF_QNM_gprior_GW_nonull(S,E,Y,response);

%%%%%%%%%%%%%%%%%%%%%%
%5
eta=0.4.*repmat(sample_locus_prop,I,1).*S+c.*E;
Y=normrnd(eta,1,I,J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Admixture:
logBLR(1,1+4*J:5*J)=logBLR_ADMIXTURE_GW_nonull(sample_prop,unifrnd(0.81,0.84,I,1),Y>repmat(quantile(Y,0.8),I,1));
%BF_QNM
logBF(1,1+4*J:5*J)=logBF_QNM_gprior_GW_nonull(S,E,Y,response);


file_name=strcat(save_dir,'nonull_E1_',idnum_str,'_cont.mat');
save(file_name,'logBLR', 'logBF');


