function logbf=logBF_QNM_gprior(Y,S,E,tau,response)
%Purpose: calculate Bayes factor for the QNM prior 
%Input:
% Y: I*1 outcome
% S: I*p AIMS
% E: I*q Environmental predictors.
% tau: 1*1 dispersion parameter
% response: 1*1 resposnse type
%Output:
% logbf: 1*1 log10(BF)

if nargin ~= 5
     error('Wrong number of input arguments')
end

if ~strcmp(response,'continuous') & ~strcmp(response,'binary') & ~strcmp(response,'count') 
    error ('response type not supported')
end

 [I,p] = size(S);
 [I,q] = size(E);
 
 if sum(E==0)==I
     X=S;
 else
     X=[S,E];
 end
 
 switch response
     case 'continuous'
         [b,dev,stats] = glmfit(X,Y,'normal','link','identity',...
                                'constant','on',... 'offset',repmat(mean(Y),I,1),...
                                'estdisp','on');
     case 'binary'
          % prev=mean(Y);
         [b,dev,stats] = glmfit(X,Y,'binomial','link','logit',...
                                'constant','on',...'offset',repmat(log(prev/(1-prev)),I,1),...
                                'estdisp','on');
     case 'count'
         [b,dev,stats] = glmfit(X,Y,'poisson','link','log',...
                                'constant','on',... 
                                'estdisp','on');
 end

 beta_hat=stats.beta(2:p+1);
 sigma2=stats.s^2;
 Sigma_hat=stats.covb(2:p+1,2:p+1)./sigma2;
  
%   beta_hat=stats.beta(1:p);
%   sigma2=stats.s^2;
%   Sigma_hat=stats.covb(1:p,1:p)./sigma2;


 T=beta_hat'*inv(Sigma_hat)*beta_hat*I*tau/(sigma2*(1+I*tau));
 logbf=log10((p+T)/(p*(1+I*tau)^(p/2+1))*exp(T/2));
