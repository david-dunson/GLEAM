%PURPOSE:
%Stochastic Search Variable Selection(SSVS) for logistic regression

%NOTES:
%SSVS is under various assumptions to make it as computational
%efficient as possible. Some assumptions may be relieved later.

%INPUT 
%Y: I*1 0/1 indicators
%S:  I*K covariates (not including intercept)
%Iter: iteration of MCMC
%Burnin: Burn in state of MCMC
%Thin: Number of thin
%tau2: configuration for SSVS
%eta2:  configuration for SSVS

%OUTPUT:
%betas_0_out:  (Iter-Burnin)/Thin*(K+1) coefficients
%deltas_out: (Iter-Burnin)/Thin*K prior class indicator

function [betas_0_out deltas_out]=logit_SSVS(Y, S, Iter, Burnin, Thin, tau2, eta2)

[I, K]=size(S);

X=[ones(I,1), S]; % add dummy covariate for the intercept;

% constants and latent vars for logistic approximation 
nu = 7.3;                       % df
sigma = pi*sqrt((nu-2)/(3*nu));    % standard deviation


%initial values
betas_0= zeros(K+1,1);
betas_0(1,1)=-7;
betas_0(2:K+1,1)=0;
phis = ones(I,1); % latent variable in normal-gamma mix
sds = sigma*phis.^(-0.5);
gs = zeros(I,1);
deltas=zeros(K,1);
lambdas=zeros(K+1,1);
lambdas(1,1)=100;
lambdas(2:K+1,1)=tau2+(eta2*tau2-tau2)*deltas;
ilambdas=diag(1./lambdas);

% saved outputs
betas_0_out=zeros((Iter-Burnin)/Thin,K+1);
deltas_out=zeros((Iter-Burnin)/Thin,K);  

I_K0=eye(K+1);

%profile on;
for g=1:Iter
    
    %Impute the latent variables gs
    etas = X*betas_0; % linear predictor in logistic model, inverse cdf methods;
    gs(Y==1) = norminv(unifrnd(normcdf(0,etas(Y==1),sds(Y==1)),ones(sum(Y),1)),etas(Y==1),sds(Y==1)); 
    gs(Y==0) = norminv(unifrnd(zeros(sum(1-Y),1),normcdf(0,etas(Y==0),sds(Y==0))),etas(Y==0),sds(Y==0));
    gs(isinf(gs)==1)=.000001*gs(isinf(gs)==1);  %???
 
    %update phis and sds
    %phis = gamrnd((nu + 1)/2, 2./(nu + sigma^(-2)*(gs-etas).^2));
    phis = randgamma(repmat((nu + 1)/2, I, 1)) .* 2./(nu + sigma^(-2)*(gs-etas).^2);
    sds = sigma*phis.^(-0.5); 
     
    % update betas
    %SD=repmat(1./sds,[1 K+1]);  
    %V_beta=inv((SD.*X)'*(SD.*X)+iLambda); % SD is not iGamma
    %E_beta=V_beta*((SD.*X)'*(gs./sds)+iLambda*mus_0);
    
    iSD2=repmat(1./sds.^2, 1, K+1); % I*(K+1)
    iSD2Xp=(iSD2.*X)';
    %V_beta=inv(iSD2Xp*X+ilambdas);
    V_beta=I_K0/(iSD2Xp*X+ilambdas); % faster version of inverse
    E_beta=V_beta*(iSD2Xp*gs);
    V_beta_chol=chol(V_beta);
    betas_0=randnorm(1,E_beta,V_beta_chol); %function from light speed package
    betas=betas_0(2:K+1);
  
    if sum(isnan(betas_0))>0
       break
    end
        
    % updata deltas
    tmp=exp(betas.^2./(2.*tau2).*(1-1/eta2))./sqrt(eta2);
    deltas=sample_vector([ones(K,1), tmp]')'-1;
    
    
    
%     [b_LS, sigma_b_LS, s_LS] = lscov(X,gs,sds.^2);
%     size(gs)
    
    % update lambdas
    lambdas(2:K+1,1)=tau2+(eta2*tau2-tau2)*deltas;
    ilambdas=diag(1./lambdas);
    
    %disp([g,sum(deltas)]);
    if (g>Burnin && rem(g/Thin,1)==0)
        
    betas_0_out((g-Burnin)/Thin,:)= betas_0';
    deltas_out((g-Burnin)/Thin,:) = deltas';

    end

end
%profile report;


