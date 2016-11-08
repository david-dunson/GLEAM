% PURPOSE 
% update PaPb matrix 
% INPUT
% S: I*J state matrix
% SNPs: I*J SNPs matrix 
% pa_0: 1*J parental population A's proportion of variant alleles
% pb_0: 1*J parental population B's proportion of variant alleles
% pa)old: 1*J previous population A's proportion of variant alleles
% pb_old: 1*J previous population B's proportion of variant alleles
% tau_a: inflate parameter for population A
% tau_b: inflate parameter for population B 
% OUTPUT
% pa: 1*J current population A's proportion of variant alleles
% pb: 1*J current population B's proportion of variant alleles

% pa_0=(afrn_var/60)'; pb_0=(caucn_var/60)'; tau_a=300; tau_b=300;
% pa_old=pa_0; pb_old=pb_0;
% SNPs2=obj.SNPs;

function [pa,pb]=PaPb_update(S,SNPs2,pa_0,pb_0,pa_old,pb_old,tau_a,tau_b)

[I,J]=size(S);
pa=zeros(1,J);
pb=zeros(1,J);

 for j=1:J
   %sx_cnt = crosstab(S(:,j), SNPs2(:,j));
   sx_cnt = my_crosstab(S(:,j), SNPs2(:,j));
   va_rb  = randbinom((pa_old(j)*(1-pb_old(j)))/(pa_old(j)*(1-pb_old(j))+pb_old(j)*(1-pa_old(j))),sx_cnt(2,2));
   pa(j)  = randbeta(tau_a*pa_0(j)+sx_cnt(3,2)+2*sx_cnt(3,3)+va_rb, ...
                     tau_a*(1-pa_0(j))+2*sx_cnt(3,1)+ sx_cnt(3,2)+sx_cnt(2,2)-va_rb);
   pb(j)  = randbeta(tau_b*pb_0(j)+sx_cnt(1,2)+2*sx_cnt(1,3)+sx_cnt(2,2)-va_rb,...
                     tau_b*(1-pb_0(j))+2*sx_cnt(1,1)+sx_cnt(1,2)+va_rb);
end
%