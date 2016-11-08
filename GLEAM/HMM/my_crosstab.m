function mat=my_crosstab(x1,x2)

idx_10=(x1==0);
idx_11=(x1==1);
idx_12=(x1==2);

idx_20=(x2==0);
idx_21=(x2==1);
idx_22=(x2==2);

mat=[sum(idx_10&idx_20),sum(idx_10&idx_21),sum(idx_10&idx_22);...
     sum(idx_11&idx_20),sum(idx_11&idx_21),sum(idx_11&idx_22);...
     sum(idx_12&idx_20),sum(idx_12&idx_21),sum(idx_12&idx_22)
    ];