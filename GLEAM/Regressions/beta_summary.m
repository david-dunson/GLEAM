function [means uppers lowers]=beta_summary(save_dir, simu_type, save_name,lower_q, upper_q)

betas_out_all=zeros(1000*100,1299);

idnum = 1;
for idnum=1:100
file_name=strcat(save_dir, simu_type,save_name,'_',num2str(idnum+900),'.mat');
load(file_name,'betas_0_out');
betas_out_all((1+(idnum-1)*1000):(1000+(idnum-1)*1000),:)=betas_0_out(1001:2000,2:end);
clear betas_0_out; 
end

means=zeros(1,1299);
uppers=zeros(1,1299);
lowers=zeros(1,1299);

for j=1:1299
means(j)=mean(betas_out_all(:,j));
uppers(j)=quantile(betas_out_all(:,j),upper_q);
lowers(j)=quantile(betas_out_all(:,j),lower_q);
end

clear betas_out_all;
