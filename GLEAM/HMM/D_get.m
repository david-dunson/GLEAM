%PURPOSE: caculate genetic distance
%INPUT:
%chrom: 1*J chromosome number
%gen_pos: 1*J genetic position in CM
%OUTPUT:
%D: 1*J disntance in Morgan(M)
function [D] = D_get(chrom, gen_pos)

D=[];

for k=1:23
    idx=(chrom==k);
    d=diff([0;gen_pos(idx)]);
    D=[D d'./100];
end

