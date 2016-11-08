function [betalist siglist] = GDP_gibbs_p_large(X,y,alpha,eta,it,burnin,thin,part)

% This function if appropriate to use when p>>n. If p is not a multiple of 
% the input variable "part", some adjustments need be done in the update of 
% the regression coefficients. eta can be fixed at 1 and alpha should be 
% chosen to adjust for multiplicity.
%
% Created: 
%   April 2011
%
% Copyright Owners:
%   Artin Armagan, David Dunson and Jaeyong Lee
%   All rights reserved

%Constants, starting values, etc...
[n p] = size(X);
tau = repmat(1e3,p,1);
a = 1;
tauinv = 1./tau;
be = zeros(p,1);
la = (alpha+1)./(abs(be)+eta);
sig = 1;

betalist = zeros((it-burnin)/thin,p);
siglist = zeros((it-burnin)/thin,1);

tic
k = 0;
for i = 1:it;
    
    beold = be;
    %Update of regression coefficients
    yrest = X*beold;
    for j = 1:(p/part)
        index = ((j-1)*part+1):((j-1)*part+part);
        Xindex = X(:,index);
        XXinv = (Xindex'*Xindex+diag(tauinv(index)))\eye(part);
        %Xdiag = Xindex.*repmat((tau(index))',n,1);
        %XXinv = diag(tau(index))-Xdiag'*((eye(n)+Xdiag*Xindex')\eye(n))*Xdiag;
        bmean = XXinv*Xindex'*(y-(yrest-Xindex*be(index)));
        bCov = sig.*XXinv;
        bCov = (bCov + bCov')/2;
        be(index) = mvnrnd(bmean,bCov);
        yrest = yrest+Xindex*(be(index)-beold(index));
    end  
    be(abs(be)<1e-6) = 1e-6;
    
    %Update of error variance
    residual = (y-X*be); 
    dhat = (residual'*residual+sum((be.^2).*tauinv))./2;
    chat = (n+p)./2;
    sig = 1./gamrnd(chat, 1./dhat);
    
    %Update of first level scale parameter
    invgauss_mean = sqrt((la.^2).*sig./(be.^2));
    invgauss_scale = la.^2;
    nu = randn(p,1);
    yy = nu.^2;
    xx = invgauss_mean+(invgauss_mean.^2).*yy./(2.*invgauss_scale) ...
        - (invgauss_mean./(2.*invgauss_scale)).*sqrt(4.*invgauss_mean ...
        .*invgauss_scale.*yy+(invgauss_mean.^2).*(yy.^2));
    zz = rand(p,1); 
    cons = invgauss_mean./(invgauss_mean+xx);
    index_1 = (zz<=cons);
    index_2 = (zz>cons);
    tauinv(index_1) = xx(index_1);
    tauinv(index_2) = (invgauss_mean(index_2).^2)./xx(index_2);

    %Update of second level scale parameter
    la = gamrnd(alpha+1, 1./(abs(be)./sqrt(sig)+eta));

    if (i>burnin && sum((i-burnin)==(1:thin:(it-burnin)))==1)
        k = k+1;
        betalist(k,:) = be;
        siglist(k) = sig;
    end

end