function simu_nonnull_single_plot2(tmp_logBLR,tmp_logBF,tmp_lassoEnet)

x=1:5;
y1=tmp_logBLR(1,:);
L1=y1-tmp_logBLR(2,:);
U1=tmp_logBLR(3,:)-y1;

y2=tmp_logBF(1,:);
L2=y2-tmp_logBF(2,:);
U2=tmp_logBF(3,:)-y2;

y3=tmp_lassoEnet(1,:);
L3=y3-tmp_lassoEnet(2,:);
U3=tmp_lassoEnet(3,:)-y3;

errorbar(x-0.2,y1,L1,U1,'.r','MarkerSize', 10)
hold on
errorbar(x-0.1,y2,L2,U2,'xb','MarkerSize', 5)
errorbar(x+0.1,y3(:,1:5),L3(:,1:5),U3(:,1:5),'ok','MarkerSize', 5)
errorbar(x+0.2,y3(:,6:10),L3(:,6:10),U3(:,6:10),'sg','MarkerSize', 5)
hold off