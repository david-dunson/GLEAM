%Simualte the null model
save_dir='C:\MyProjects\Admixture\Toshare\Simulations\SimulatedData\';
data_dir='C:\MyProjects\Admixture\Toshare\Simulations\';

%tic;
%idnum=1;
for idnum = 1:1 %for 1 simulation dataset only. In the paper, we simulated 100 dataset.


RandState=4000;
rand('state',RandState+idnum*1000); % set arbitrary seed for uniform draws
randn('state',RandState+idnum*1000);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=1000;
J=1000; %For 1000 makers only. In the paper, we set J=10000, because the type I error is quit small

file_name=strcat(data_dir,'Afr_prop.mat'); % Afr_prop for each locus: 1001 subjects, 1296 AIM
load(file_name);

locus_prop=mean(Afr_prop,1);

sub_prop=mean(Afr_prop,2);

top20_idx=locus_prop>quantile(locus_prop,0.8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_all = zeros(I,J);
sample_prop=zeros(I,J);

idx=randsample(find(~top20_idx),J,true); % sample with replacement for J loci

for j=1:J
   sample_prop(:,j)=randsample(Afr_prop(:,idx(j)),I,true); % sample with replacement for I loci
end

for j=1:J
S=mnrnd(1,[(1-sample_prop(:,j)).^2,2.*sample_prop(:,j).*(1-sample_prop(:,j)),sample_prop(:,j).^2])*[0;1;2];
S_all(:,j)=S-mean(S); %sample the ancestry number 0 1 or 2;
end
% Mean zero S

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 0;
E0=normrnd(0,1,I,1);
eta=c*E0; % without Enviomental effects


Y=binornd(1,exp(eta)./(1+exp(eta)));
Y_binary_E0=(Y==1);

Y_cont_E0=normrnd(eta,1,I,1);

c = 1;
E1=normrnd(0,1,I,1);
eta=c*E1; % with Enviomental effects


Y=binornd(1,exp(eta)./(1+exp(eta)));
Y_binary_E1=(Y==1);

Y_cont_E1=normrnd(eta,1,I,1);

file_name=strcat(save_dir,'simu_null_',num2str(idnum),'.csv');
csvwrite(file_name,[Y_binary_E0,Y_binary_E1,Y_cont_E0,Y_cont_E1,E0,E1,S_all]);
end
%toc;