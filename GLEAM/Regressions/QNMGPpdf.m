function f=QNMGPpdf(beta, tau,sigma2,X)

% purpose:
% pdf for the quadratic normal moment g prior;
% input: 
% beta: q*1;
% tau, sigma2: scaler;
% X: n*q;

% ouput: 
% f: scaler density

[n,q]=size(X);

MU=zeros(q,1);
SIGMA=n*tau*sigma2*inv(X'*X);

f=beta'*(X'*X)*beta*mvnpdf(beta,MU,SIGMA)/(n*tau*sigma2*q);




