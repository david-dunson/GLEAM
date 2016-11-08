save_dir='C:\MyProjects\Admixture\Toshare\Simulations\SimulatedData\';
data_dir='C:\MyProjects\Admixture\Toshare\Simulations\';

for idnum = 1:1; %for 1 simulation dataset only. In the paper, we simulated 100 dataset.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%single loci nonull models: binary, no E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RandState=4000;
c=0;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y_binary_nocov=zeros(I,5*J);

%impute the nonnull models:
idx=randsample(find(top20_idx),J,true);

for j=1:J
   sample_prop(:,j)=randsample(Afr_prop(:,idx(j)),I,true);
end

sample_locus_prop=mean(sample_prop,1); % population average for each locus

S=zeros(I,J);
for j=1:J
S(:,j)=mnrnd(1,[(1-sample_prop(:,j)).^2,2.*sample_prop(:,j).*(1-sample_prop(:,j)),sample_prop(:,j).^2])*[0;1;2];
end;
S=S-repmat(mean(S),I,1);


E=normrnd(0,1,I,J);
%%%%%%%%%%%%%%%%%%%%%%
%1
eta=0.4.*repmat(sample_locus_prop,I,1).*S;
Y=binornd(1,exp(eta)./(1+exp(eta)));
Y_binary_nocov(:,1:100)=(Y==1);

%%%%%%%%%%%%%%%%%%%%%%
%2
eta=0.5.*repmat(sample_locus_prop,I,1).*S;
Y=binornd(1,exp(eta)./(1+exp(eta)));
Y_binary_nocov(:,(1:100)+1*100)=(Y==1);

%%%%%%%%%%%%%%%%%%%%%%
%3
eta=0.6.*repmat(sample_locus_prop,I,1).*S;
Y=binornd(1,exp(eta)./(1+exp(eta)));
Y_binary_nocov(:,(1:100)+2*100)=(Y==1);

%%%%%%%%%%%%%%%%%%%%%%
%4
eta=0.7.*repmat(sample_locus_prop,I,1).*S;
Y=binornd(1,exp(eta)./(1+exp(eta)));
Y_binary_nocov(:,(1:100)+3*100)=(Y==1);


%%%%%%%%%%%%%%%%%%%%%%
%5
eta=0.8.*repmat(sample_locus_prop,I,1).*S;
Y=binornd(1,exp(eta)./(1+exp(eta)));
Y_binary_nocov(:,(1:100)+4*100)=(Y==1);

file_name=strcat(save_dir,'simu_nonnull_binary_nocov_',num2str(idnum),'.csv');
csvwrite(file_name,[Y_binary_nocov,S,E]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%single loci nonull models: binary, with E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y_binary_cov=zeros(I,5*J);

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
eta=0.4.*repmat(sample_locus_prop,I,1).*S+c.*E;
Y=binornd(1,exp(eta)./(1+exp(eta)));
Y_binary_cov(:,1:100)=(Y==1);

%%%%%%%%%%%%%%%%%%%%%%
%2
eta=0.5.*repmat(sample_locus_prop,I,1).*S+c.*E;
Y=binornd(1,exp(eta)./(1+exp(eta)));
Y_binary_cov(:,(1:100)+1*100)=(Y==1);

%%%%%%%%%%%%%%%%%%%%%%
%3
eta=0.6.*repmat(sample_locus_prop,I,1).*S+c.*E;
Y=binornd(1,exp(eta)./(1+exp(eta)));
Y_binary_cov(:,(1:100)+2*100)=(Y==1);

%%%%%%%%%%%%%%%%%%%%%%
%4
eta=0.7.*repmat(sample_locus_prop,I,1).*S+c.*E;
Y=binornd(1,exp(eta)./(1+exp(eta)));
Y_binary_cov(:,(1:100)+3*100)=(Y==1);


%%%%%%%%%%%%%%%%%%%%%%
%5
eta=0.8.*repmat(sample_locus_prop,I,1).*S+c.*E;
Y=binornd(1,exp(eta)./(1+exp(eta)));
Y_binary_cov(:,(1:100)+4*100)=(Y==1);

file_name=strcat(save_dir,'simu_nonnull_binary_cov_',num2str(idnum),'.csv');
csvwrite(file_name,[Y_binary_cov,S, E]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%single loci nonull models: continuous NO E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RandState=4000;
c=0;
rand('state',RandState+idnum*100+c*100000); % set arbitrary seed for uniform draws
randn('state',RandState+idnum*100+c*100000);
 
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
Y_con_nocov=zeros(I,5*J);
%%%%%%%%%%%%%%%%%%%%%%
%1
eta=0.2.*repmat(sample_locus_prop,I,1).*S;
Y_con_nocov(:,1:100)=normrnd(eta,1,I,J);

%%%%%%%%%%%%%%%%%%%%%%
%2
eta=0.25.*repmat(sample_locus_prop,I,1).*S;
Y_con_nocov(:,(1:100)+1*100)=normrnd(eta,1,I,J);

%%%%%%%%%%%%%%%%%%%%%%
%3
eta=0.3.*repmat(sample_locus_prop,I,1).*S;
Y_con_nocov(:,(1:100)+2*100)=normrnd(eta,1,I,J);

%%%%%%%%%%%%%%%%%%%%%%
%4
eta=0.35.*repmat(sample_locus_prop,I,1).*S;
Y_con_nocov(:,(1:100)+3*100)=normrnd(eta,1,I,J);

%%%%%%%%%%%%%%%%%%%%%%
%5
eta=0.4.*repmat(sample_locus_prop,I,1).*S;
Y_con_nocov(:,(1:100)+4*100)=normrnd(eta,1,I,J);

file_name=strcat(save_dir,'simu_nonnull_cont_nocov_',num2str(idnum),'.csv');
csvwrite(file_name,[Y_con_nocov,S,E]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%single loci nonull models: continuous with E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RandState=4000;
c=1;
rand('state',RandState+idnum*100+c*100000); % set arbitrary seed for uniform draws
randn('state',RandState+idnum*100+c*100000);
 
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
Y_con_cov=zeros(I,5*J);
%%%%%%%%%%%%%%%%%%%%%%
%1
eta=0.2.*repmat(sample_locus_prop,I,1).*S+c.*E;
Y_con_cov(:,1:100)=normrnd(eta,1,I,J);

%%%%%%%%%%%%%%%%%%%%%%
%2
eta=0.25.*repmat(sample_locus_prop,I,1).*S+c.*E;
Y_con_cov(:,(1:100)+1*100)=normrnd(eta,1,I,J);

%%%%%%%%%%%%%%%%%%%%%%
%3
eta=0.3.*repmat(sample_locus_prop,I,1).*S+c.*E;
Y_con_cov(:,(1:100)+2*100)=normrnd(eta,1,I,J);

%%%%%%%%%%%%%%%%%%%%%%
%4
eta=0.35.*repmat(sample_locus_prop,I,1).*S+c.*E;
Y_con_cov(:,(1:100)+3*100)=normrnd(eta,1,I,J);

%%%%%%%%%%%%%%%%%%%%%%
%5
eta=0.4.*repmat(sample_locus_prop,I,1).*S+c.*E;
Y_con_cov(:,(1:100)+4*100)=normrnd(eta,1,I,J);

file_name=strcat(save_dir,'simu_nonnull_cont_cov_',num2str(idnum),'.csv');
csvwrite(file_name,[Y_con_cov,S,E]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%two loci nonull models: binary/cont, no E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idnum = 1:1;

RandState=4000;
rand('state',RandState+idnum); % set arbitrary seed for uniform draws
randn('state',RandState+idnum);

file_name=strcat(data_dir,'Afr_prop.mat');
load(file_name);
sample_prop=Afr_prop(:,[15:65,266:316]);

[I,J]=size(sample_prop);
S=zeros(I,J);
E=zeros(I,1);

for j=1:J
S(:,j)=mnrnd(1,[(1-sample_prop(:,j)).^2,2.*sample_prop(:,j).*(1-sample_prop(:,j)),sample_prop(:,j).^2])*[0;1;2];
end;

S=S-repmat(mean(S),I,1);

%binary
eta=0.7.*S(:,[25,75])*locus_prop([39,289])';

Y=binornd(1,exp(eta)./(1+exp(eta)));
Y_binary=(Y==1);
S_binary = S;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=zeros(I,J);
E=zeros(I,1);

for j=1:J
S(:,j)=mnrnd(1,[(1-sample_prop(:,j)).^2,2.*sample_prop(:,j).*(1-sample_prop(:,j)),sample_prop(:,j).^2])*[0;1;2];
end;

S=S-repmat(mean(S),I,1);
%cont
eta=0.35.*S(:,[25,75])*locus_prop([39,289])';

Y_cont=normrnd(eta,1,I,1); % GAM use continuous.
S_cont=S;
file_name=strcat(save_dir,'simu_nonnull_2loci_',num2str(idnum),'.csv');
csvwrite(file_name,[Y_binary,Y_cont,S_binary,S_cont]);
end





