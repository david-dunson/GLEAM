

la=[1;2;3]
ra=[4;5;6]
mm=magic(3)

tic;
for k=1:10000
tmp=repmat(la,1,3) .* mm .* repmat(ra',3,1);
end
toc;

tic;
for k=1:1000
tmp=FWkernal(la,mm,ra);
end
toc;

cd('E:\dropbox\DropBox\Dropbox\MATLAB\MatlabToolbox\MyFunctions\Admixture\HMM')
mex FWkernal.c