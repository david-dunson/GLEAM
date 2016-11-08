function tab=nonnull_summary(intab, cutoff)

idx_sig=intab>cutoff;

[I,K]=size(intab);
J=K/5;

part1=idx_sig(:,1:J);
part2=idx_sig(:,1+J:J+J);
part3=idx_sig(:,1+J+J:J+J+J);
part4=idx_sig(:,1+3*J:4*J);
part5=idx_sig(:,1+4*J:5*J);


summay1=[median(sum(part1,2))/J;min(sum(part1,2))/J;max(sum(part1,2))/J];
summay2=[median(sum(part2,2))/J;min(sum(part2,2))/J;max(sum(part2,2))/J];
summay3=[median(sum(part3,2))/J;min(sum(part3,2))/J;max(sum(part3,2))/J];
summay4=[median(sum(part4,2))/J;min(sum(part4,2))/J;max(sum(part4,2))/J];
summay5=[median(sum(part5,2))/J;min(sum(part5,2))/J;max(sum(part5,2))/J];

tab=[summay1,summay2,summay3,summay4,summay5];
