%PURPOSE:
%Calcualte tau and eta given negative predictive value(npv) and 
%positive predictive value(ppv) 

function [tau eta] =tau_eta(beta, npv, ppv)

eta=npv/(1-npv);

tau=abs(beta)*sqrt((1-1/eta^2)/(2*log(eta*ppv/(1-ppv))));