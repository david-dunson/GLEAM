function [betas_0_out thetas_out Sigmas_out zs_out taus_out gdp_alpha_out gdp_eta_out] = ASPR_GDP(Y,X_betas,zs,Gsim,burnin,thin,gdp_alpha, gdp_eta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PURPOSE: MCMC for two-component multivariate latent class model with weight depending on
%         high-deminsional covaraites and each component following
%         mulitivariate normal distributions.
%INPUT:   Y        n by s responses
%         X_betas  n by p covaraites
%         zs       n by 1 initial classification
%         Gsim     1 by 1 number of simulations
%         burnin   1 by 1 number of burn in 
%         thin     1 by 1 number of thin
%OUTPUT:  betas_0_out  (Gsim-burnin)/thin by p+1 betas
%         thetas_out   (Gsim-burnin)/thin by 2*s thetas
%         Sigmas_out   (Gsim-burnin)/thin by s*(s+1) Sigmas
%                       upper triangular by rows 
%         zs_out       (Gsim-burnin)/thin by n zs

%taus_out:     (Gsim-burnin)/thin*(p) coefficients
%gdp_alpha_out (Gsim-burnin)/thin*(1) coefficients
%gdp_eta_out  (Gsim-burnin)/thin*(1) coefficients

%         tout=[betas_0_out,thetas_out,Sigmas_out,zs_out]
%Writen by: Bin Zhu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[n,s]=size(Y);
[n,p]=size(X_betas);

X=[ones(n,1),X_betas]; % add dummy covariate for the intercept;

  c0=0;      %  prior mean for intercept
  cv0=100;    % prior variance for intercept
  %a=1;        % alpha value
  
  %construct prior for the latent class mean and covariance matrix;
  rho_1=6;
  rho_2=rho_1; 
  psi_1=1;
  psi_2=psi_1;
  
  thetas_0=mean(Y);
  Sigma_0=(rho_1-s-1)/2.*cov(Y);
  
  thetas_0_1=thetas_0;
  Sigma_0_1=Sigma_0;  
  thetas_0_2=thetas_0; 
  Sigma_0_2=Sigma_0;

  
  % constants and latent vars for logistic approximation 
  nu = 7.3;                       % df
  sigma = pi*sqrt((nu-2)/(3*nu));    % standard deviation
  phis = ones(n,1);                % latent variable in normal-gamma mix
  sds = sigma*phis.^(-0.5);
  gs = zeros(n,1);

  betas_0= zeros(p+1,1);
  %muG = 0;  
  %ks=ones(p,1); % Bin index for DP;
  %mus=zeros(p,1);
  %mus_0=[c0,mus']';
  %p_star=max(ks);              %no. grps
  %Ssum=(repmat(ks(:,1),[1 p_star])-repmat(1:p_star,[p 1])==0);
  %ms = ones(1,p)*Ssum;              %no in each grp.
  
  %MCMC saves
  
  betas_0_out=zeros((Gsim-burnin)/thin,p+1);
  %mus_out=zeros((Gsim-burnin)/thin,p);
  %taus_out=zeros((Gsim-burnin)/thin,p);
  
  thetas_out=zeros((Gsim-burnin)/thin,2*s);
  Sigmas_out=zeros((Gsim-burnin)/thin,s*(s+1));
  
  zs_out=zeros((Gsim-burnin)/thin,n);
  
  taus_out=zeros((Gsim-burnin)/thin,p);
  gdp_alpha_out=zeros((Gsim-burnin)/thin,1);
  gdp_eta_out=zeros((Gsim-burnin)/thin,1);
  
  %phis_out=zeros((Gsim-burnin)/thin,length(zs));
  %lambdas_out=zeros((Gsim-burnin)/thin,p);
  
  %ks_out=zeros((Gsim-burnin)/thin,p);
  
  %lambdas=ones(p,1);
  %iLambda=inv(diag([cv0,lambdas']));
  %XpX=X'*X;
  
  %taus_star=1;
  %mus_star=0;
  %taus=ones(p,1);
  
  taus=ones(p+1,1);
  taus(1,1)=100;
  itaus=diag(1./taus);
  lambdas=ones(p,1);
  
  I_K0=eye(p+1);

  as=(0.0001:0.002:0.9999)';
  es=(0.0001:0.002:0.9999)';
  L=length(es);

for g=1:Gsim
  % --------------------------------------------------------------------- %
  % update thetas_1/2 and Sigma_1/2;
  % -------------------------------------------------------------------- % 
  Y_bar_1=mean(Y(find(zs==1),:));
  n_1=sum(zs==1);
  tmp=Y(find(zs==1),:)- Y_bar_1(ones(1,n_1),:);
  S_1=tmp'*tmp;
  
  tmp=n_1/(n_1+psi_1);
  theta_hat_1=tmp*Y_bar_1+(1-tmp)*thetas_0_1;
  psi_hat_1=psi_1+n_1;
  rho_hat_1=rho_1+n_1;
  Sigma_hat_1=Sigma_0_1+S_1+n_1/(1+n_1/psi_1).*(Y_bar_1-thetas_0)'*(Y_bar_1-thetas_0);
  
  Y_bar_2=mean(Y(find(zs==0),:));
  n_2=sum(zs==0);
  tmp=Y(find(zs==0),:)- Y_bar_2(ones(1,n_2),:);
  S_2=tmp'*tmp;
  
  tmp=n_2/(n_2+psi_2);
  theta_hat_2=tmp*Y_bar_2+(1-tmp)*thetas_0_2;
  psi_hat_2=psi_2+n_2;
  rho_hat_2=rho_2+n_2;
  Sigma_hat_2=Sigma_0_2+S_2+n_2/(1+n_2/psi_2).*(Y_bar_2-thetas_0)'*(Y_bar_2-thetas_0);

  Sigma_1=iwishrnd(Sigma_hat_1,psi_hat_1);
  thetas_1=randnorm(1,theta_hat_1',chol(Sigma_1./psi_hat_1) )'; 
 
  Sigma_2=iwishrnd(Sigma_hat_2,psi_hat_2);
  thetas_2=randnorm(1,theta_hat_2',chol(Sigma_2./psi_hat_2) )';
  
  while(sum(thetas_1>0)~=s | sum(thetas_2>0)~=s | thetas_1(1) >= thetas_2(1) | thetas_1(2) >= thetas_2(2))
  
  Sigma_1=iwishrnd(Sigma_hat_1,psi_hat_1);
  thetas_1=randnorm(1,theta_hat_1',chol(Sigma_1./psi_hat_1) )'; 
 
  Sigma_2=iwishrnd(Sigma_hat_2,psi_hat_2);
  thetas_2=randnorm(1,theta_hat_2',chol(Sigma_2./psi_hat_2) )';
   end;
  
  % --------------------------------------------------------------------- %
  % impute latent variables                                               %
  % --------------------------------------------------------------------- %
  etas = X*betas_0; % linear predictor in logistic model, inverse cdf methods;
  gs(zs==1) = norminv(unifrnd(normcdf(0,etas(zs==1),sds(zs==1)),ones(sum(zs),1)),etas(zs==1),sds(zs==1)); 
  gs(zs==0) = norminv(unifrnd(zeros(sum(1-zs),1),normcdf(0,etas(zs==0),sds(zs==0))),etas(zs==0),sds(zs==0));
  gs(isinf(gs)==1)=.000001*zs(isinf(gs)==1);  %???
  
  % --------------------------------------------------------------------- %
  % Impute zs
  % -------------------------------------------------------------------- %
  prob  = 1-normcdf(0,etas,sds); % prob of z_i=1
  prop1 = prob.*mvnpdf(Y,thetas_1,Sigma_1);
  prop2 = (1-prob).*mvnpdf(Y,thetas_2,Sigma_2);
  post.prob = prop1./(prop1+prop2);
  
  zs=binornd(1,post.prob);% impute zs from bernoulli
    
  
  % --------------------------------------------------------------------- %
  % update phis and sds                                                   %
  % --------------------------------------------------------------------- %
  phis = gamrnd((nu + 1)/2, 2./(nu + sigma^(-2)*(gs-etas).^2));
  sds = sigma*phis.^(-0.5); 
  %pdi = sds.^(-2);

  if sum(isnan(phis))>0
      break
  end
  
  % --------------------------------------------------------------------- %
  % Update intercept and betas'
  % -------------------------------------------------------------------- %
   % update betas
    %SD=repmat(1./sds,[1 K+1]);  
    %V_beta=inv((SD.*X)'*(SD.*X)+iLambda); % SD is not iGamma
    %E_beta=V_beta*((SD.*X)'*(gs./sds)+iLambda*mus_0);
    
    iSD2=repmat(1./sds.^2, 1, p+1); % n*(p+1)
    iSD2Xp=(iSD2.*X)';
    %V_beta=inv(iSD2Xp*X+ilambdas);
    V_beta=I_K0/(iSD2Xp*X+itaus); % faster version of inverse
    E_beta=V_beta*(iSD2Xp*gs);
    V_beta_chol=chol(V_beta);
    betas_0=randnorm(1,E_beta,V_beta_chol); %function from light speed package
    betas_0(abs(betas_0)<1e-6) = 1e-6;
    betas=betas_0(2:p+1);
  
    if sum(isnan(betas_0))>0
       break
    end
        
    % update taus
    tmp3=igrnd(abs(lambdas ./ betas),lambdas.^2,p);
    itaus=diag([1/100;tmp3]);
    
    % update lambdas
    lambdas = randgamma(repmat((gdp_alpha + 1), p, 1)) .* 1./(abs(betas) + repmat(gdp_eta, p,1));
      
    %update gdp_alpha
    tmp=p*log(1./as-1) - sum(log( 1 + abs(betas)./gdp_eta ))./as;
    gdp_a=as(sample(exp(tmp-max(tmp)),1));
    gdp_alpha=1/gdp_a-1;
    
    %update gdp_eta
   tmp=p*log(es./(1.-es))-(gdp_alpha+1)*sum( log(repmat((es./(1.-es)),1,p).*repmat(abs(betas)',L,1)+1),2);
   gdp_e=es(sample(exp(tmp-max(tmp)),1)); 
   gdp_eta=1/gdp_e-1;
%     [b_LS, sigma_b_LS, s_LS] = lscov(X,gs,sds.^2);
%     size(gs)
  % --------------------------------------------------------------------- %
  % Save outputs
  % -------------------------------------------------------------------- %
  if (g>burnin & rem(g/thin,1)==0)
    
    
    betas_0_out((g-burnin)/thin,:)=betas_0';
    %mus_out((g-burnin)/thin,:)=mus';
    %taus_out((g-burnin)/thin,:)=taus';
    
    thetas_out((g-burnin)/thin,:)=[thetas_1,thetas_2];
    
    tmp1=reshape(triu(Sigma_1)',1,s*s);
    tmp2=reshape(triu(Sigma_2)',1,s*s);
    Sigmas_out((g-burnin)/thin,:)=[tmp1(find(tmp1~=0)),tmp2(find(tmp2~=0))];
    
    zs_out((g-burnin)/thin,:)=zs';
    taus_out((g-burnin)/thin,:)=(1./tmp3)';
    gdp_alpha_out((g-burnin)/thin,:)=gdp_alpha;
    gdp_eta_out((g-burnin)/thin,:)=gdp_eta;
  
  end
  
  
  %
  
end



