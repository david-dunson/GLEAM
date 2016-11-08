function logBF=logBF_QNM_gprior_GW_nonull(S,E,Y,response)

J=size(S,2);

logBF=zeros(1,J);

for j=1:J
[x,fval]=fminbnd(@(x)-logBF_QNM_gprior(Y(:,j),S(:,j),E(:,j),x,response),0.0001,0.2);
logBF(j)=-fval;
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [I,p] = size(S(:,j));
% [I,q] = size(E(:,j));
%  
%  if sum(E(:,j)==0)==I
%      X=S(:,j);
%  else
%      X=[S(:,j),E(:,j)];
%  end
%  
%   [b,dev,stats] = glmfit(X,Y(:,j),'binomial','link','logit',...
%                                 'constant','on',...'offset',repmat(log(prev/(1-prev)),I,1),...
%                                 'estdisp','on');
%                             
%  beta_hat=stats.beta(2:p+1);
%  sigma2=stats.s^2;
%  Sigma_hat=stats.covb(2:p+1,2:p+1)./sigma2;
%  
% 
%  T=beta_hat'*inv(Sigma_hat)*beta_hat*I*tau/(sigma2*(1+I*tau));
%  logbf=log10((p+T)/(p*(1+I*tau)^(p/2+1))*exp(T/2));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%