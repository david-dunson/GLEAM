% Purpose:
% run multiloci simualtion  in Cluster
 
idnum_str = getenv('SGE_TASK_ID');
idnum = str2num(idnum_str);
 
%outcome=getenv('JOB_NAME');
addpath(genpath('/home/chg/bz27/MATLAB'))
data_dir='/home/chg/bz27/Admixture/Datasets/Sdata/Real/';
save_dir=strcat('/home/chg/bz27/Admixture/Results/');

RandState=4000;
rand('state',RandState+idnum); % set arbitrary seed for uniform draws
randn('state',RandState+idnum);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_name=strcat(data_dir,'Afr_prop.mat');
load(file_name);
 
locus_prop=mean(Afr_prop,1);
 
sub_prop=mean(Afr_prop,2);
 
top20_idx=locus_prop>quantile(locus_prop,0.8);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Multiple loci: nonnull at 289 and 39  binary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[I,J]=size(Afr_prop);

E=zeros(I,1);

file_name=strcat(data_dir, 'SR250.mat');
load(file_name);

S=S-repmat(mean(S),I,1);

eta=0.7.*S(:,[39,289])*locus_prop([39,289])';

Y=binornd(1,exp(eta)./(1+exp(eta)));
Y=(Y==1);
cutoff=2;

%Admixture:
logBLR=logBLR_ADMIXTURE_GW(S, sub_prop,Y);
%BF_QNM
logBF=logBF_QNM_gprior_GW2(S,E,Y,'binary');

max_idx=0;
max_logBF=0;

if(sum(logBF>2)>0)
[max_idx,max_logBF]=logBF_multiloci(logBF,Y,S,'binary',cutoff);
end;
file_name=strcat(save_dir,'multiloci_E0_',idnum_str,'_binary2.mat');
save(file_name,'Y', 'logBLR', 'logBF','max_idx','max_logBF');
 
