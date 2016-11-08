function uni_pMOM = uni_pMOMpdf(theta, sigma2, tau, r)

uni_pMOM=normpdf(theta, 0, sqrt(tau*sigma2))*theta^(2*r)/((tau*sigma2)^r*fact2(2*1-1));