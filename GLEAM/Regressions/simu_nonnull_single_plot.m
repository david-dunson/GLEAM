function simu_nonnull_single_plot(tmp_logBLR,tmp_logBF)

x=1:5;
y1=tmp_logBLR(1,:);
L1=y1-tmp_logBLR(2,:);
U1=tmp_logBLR(3,:)-y1;

y2=tmp_logBF(1,:);
L2=y2-tmp_logBF(2,:);
U2=tmp_logBF(3,:)-y2;

errorbar(x-0.1,y1,L1,U1,'.r','MarkerSize', 10)
hold on
errorbar(x+0.1,y2,L2,U2,'xb','MarkerSize', 5)
hold off