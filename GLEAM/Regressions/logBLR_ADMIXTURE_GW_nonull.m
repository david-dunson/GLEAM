function logBLR=logBLR_ADMIXTURE_GW_nonull(sample_prop, sample_sub_prop,casegroup)

[I,J]=size(sample_prop);

logBLR=zeros(1,J);
r=0.7:0.2:2.5;

global_prop=[(1-sample_sub_prop).^2,...
            2.*sample_sub_prop.*(1-sample_sub_prop),...
            sample_sub_prop.^2];
        
for j=1:J
local_prop=[(1-sample_prop(:,j)).^2,...
            2.*sample_prop(:,j).*(1-sample_prop(:,j)),...
            sample_prop(:,j).^2];

logBLR(j)=logBLR_ADMIXTURE(local_prop(casegroup(:,j),:), global_prop(casegroup(:,j),:), r);
end;