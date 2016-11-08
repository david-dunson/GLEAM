function logBF=logBF_QNM_gprior_GW(sample_prop,E,Y,response)

J=size(sample_prop,2);

logBF=zeros(1,J);

for j=1:J
S=mnrnd(1,[(1-sample_prop(:,j)).^2,2.*sample_prop(:,j).*(1-sample_prop(:,j)),sample_prop(:,j).^2])*[0;1;2];
S=S-mean(S); %sample the ancestry number 0 1 or 2;

[x,fval]=fminbnd(@(x)-logBF_QNM_gprior(Y,S,E,x,response),0.0001,0.2);
logBF(j)=-fval;
end   