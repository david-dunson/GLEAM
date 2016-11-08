function [means_lasso]=beta_lasso_summary(save_dir, simu_type, save_name)

betas_lasso_all=zeros(100,1299);

idnum = 1;
for idnum=1:100
file_name=strcat(save_dir, simu_type,save_name,'_lasso_',num2str(900+idnum),'.mat');
load(file_name,'beta_lasso');
betas_lasso_all(idnum,:)=beta_lasso';
end

means_lasso=mean(betas_lasso_all);

