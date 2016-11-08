%PURPOSE:
%Generalized doubel parato shrinkage(GDP) for logistic regression
%NOTES:
%INPUT 
%Y: I*1 0/1 indicators
%S:  I*K covariates (not including intercept)
%Iter: iteration of MCMC
%Burnin: Burn in state of MCMC
%Thin: Number of thin
%gdp_alpha: configuration for GDP
%dgp_eta:  configuration for GDP
%OUTPUT:
%betas_0_out:  (Iter-Burnin)/Thin*(K+1) coefficients
%taus_out:     (Iter-Burnin)/Thin*(K) coefficients
%gdp_alpha_out (Iter-Burnin)/Thin*(1) coefficients
%gdp_eta_out  (Iter-Burnin)/Thin*(1) coefficients
function [betas_0_out taus_out gdp_alpha_out gdp_eta_out]=logit_GDP(Y, S, Iter, Burnin, Thin, gdp_alpha, gdp_eta)

[I, K]=size(S);

X=[ones(I,1), S]; % add dummy covariate for the intercept;

% constants and latent vars for logistic approximation 
nu = 7.3;                       % df
sigma = pi*sqrt((nu-2)/(3*nu));    % standard deviation


%initial values
betas_0= zeros(K+1,1); 
betas_0(1,1)=-7;
%betas_0(2:K+1,1)=0;
phis = ones(I,1); % latent variable in normal-gamma mix
sds = sigma*phis.^(-0.5);
gs = zeros(I,1);

taus=ones(K+1,1);
taus(1,1)=100;
itaus=diag(1./taus);
lambdas=ones(K,1);

% saved outputs
betas_0_out=zeros((Iter-Burnin)/Thin,K+1);
taus_out=zeros((Iter-Burnin)/Thin,K); 
gdp_alpha_out=zeros((Iter-Burnin)/Thin,1);  
gdp_eta_out=zeros((Iter-Burnin)/Thin,1);  

I_K0=eye(K+1);

as=(0.0001:0.002:0.9999)';
es=(0.0001:0.002:0.9999)';
L=length(es);

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
    V_beta=I_K0/(iSD2Xp*X+itaus); % faster version of inverse
    E_beta=V_beta*(iSD2Xp*gs);
    V_beta_chol=chol(V_beta);
    betas_0=randnorm(1,E_beta,V_beta_chol); %function from light speed package
    betas_0(abs(betas_0)<1e-6) = 1e-6;
    betas=betas_0(2:K+1);
  
    if sum(isnan(betas_0))>0
       break
    end
        
    % update taus
    tmp2=igrnd(abs(lambdas ./ betas),lambdas.^2,K);
    itaus=diag([1/100;tmp2]);
    
    % update lambdas
    lambdas = randgamma(repmat((gdp_alpha + 1), K, 1)) .* 1./(abs(betas) + repmat(gdp_eta, K,1));
      
    %update gdp_alpha
    tmp=K*log(1./as-1) - sum(log( 1 + abs(betas)./gdp_eta ))./as;
    gdp_a=as(sample(exp(tmp-max(tmp)),1));
    gdp_alpha=1/gdp_a-1;
    
    %update gdp_eta
   tmp=K*log(es./(1.-es))-(gdp_alpha+1)*sum( log(repmat((es./(1.-es)),1,K).*repmat(abs(betas)',L,1)+1),2);
   gdp_e=es(sample(exp(tmp-max(tmp)),1)); 
   gdp_eta=1/gdp_e-1;
%     [b_LS, sigma_b_LS, s_LS] = lscov(X,gs,sds.^2);
%     size(gs)
        
    %disp([g,sum(deltas)]);
    if (g>Burnin && rem(g/Thin,1)==0)
        
    betas_0_out((g-Burnin)/Thin,:)= betas_0';
    taus_out((g-Burnin)/Thin,:)= (1./tmp2)';
    gdp_alpha_out((g-Burnin)/Thin,:)=gdp_alpha;
    gdp_eta_out((g-Burnin)/Thin,:)=gdp_eta;
    end

end
%profile report;


