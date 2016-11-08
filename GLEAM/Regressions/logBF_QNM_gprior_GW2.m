function logBF=logBF_QNM_gprior_GW2(S,E,Y,response)

J=size(S,2);

logBF=zeros(1,J);

for j=1:J
[x,fval]=fminbnd(@(x)-logBF_QNM_gprior(Y,S(:,j),E,x,response),0.0001,0.2);
logBF(j)=-fval;
end 