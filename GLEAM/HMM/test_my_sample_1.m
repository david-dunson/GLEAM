cd('E:\dropbox\DropBox\Dropbox\MATLAB\MatlabToolbox\MyFunctions\Admixture\HMM')
mex my_sample_1.c

draws=zeros(10000,1);

prop=[1,2,7];

tic;
for k=1:10000
draws(k)=my_sample_1(prop);
end
toc;

tabulate(draws)

tic;
for k=1:10000
draws(k)=my_sample_1_old(prop);
end
toc;


