function pr=QNM_int_prob(a,tau)

pr=2*(normcdf(a/sqrt(tau),0,1)-a*exp(-a^2/(2*tau))/sqrt(2*3.1415926*tau)-0.5);