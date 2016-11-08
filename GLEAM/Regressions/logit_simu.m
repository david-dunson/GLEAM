function zs = logit_simu(betas, X_betas, RandState, intercept)

rand('state',RandState); % set arbitrary seed for uniform draws
randn('state',RandState); % set arbitrary seed for normal draws

[n p]=size(X_betas);

betas_0=[intercept, betas']';
X=[ones(n,1) X_betas];

omegas_1=exp(X*betas_0)./(1+exp(X*betas_0));
zs=binornd(1,omegas_1);

