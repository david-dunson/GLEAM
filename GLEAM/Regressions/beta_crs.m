function [probs_nonnull]=beta_crs(save_dir, simu_type, save_name)

betas_out_all=zeros(1000*100,1299);

idnum = 1;
for idnum=1:100
file_name=strcat(save_dir, simu_type,save_name,'_',num2str(idnum+900),'.mat');
load(file_name,'betas_0_out');
betas_out_all((1+(idnum-1)*1000):(1000+(idnum-1)*1000),:)=betas_0_out(1001:2000,2:end);
clear betas_0_out; 
end

epss=0.01:0.01:0.2;
probs_nonnull=zeros(length(epss),1299);
%[100, 200,201, 300, 301,302, 400,401,402,403, 500,501,502,503,504]
for j=1:1299
    probs_nonnull(:,j)=mean(repmat(betas_out_all(:,j),1,length(epss)) > repmat(epss,1000*100,1))'; 
end


