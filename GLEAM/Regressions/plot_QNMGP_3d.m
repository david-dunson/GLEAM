%for 2-dim
cut0=norminv((1-pa)^2,0,1);
cut1=norminv(1-pa^2,0,1);

mu = [0 0];
SIGMA = [1 0; 0 1];
Z = mvnrnd(mu,SIGMA,n);

X=ones(n,2);
X(Z<=cut0)=0;
X(Z>cut1)=2;

count(X(:,1));
count(X(:,2));

X=X-repmat(mean(X),n,1); %centered;

betas=(-2.5:0.05:2.5)';

fs=zeros(length(betas),length(betas));

for i = 1:length(betas)
    for j = 1:length(betas)
        fs(i,j)=QNMGPpdf([betas(i); betas(j)], tau,sigma2,X);
    end
end

surfc(betas,betas,fs)
shading interp
axis([-2.5 2.5 -2.5 2.5 -0.2 0.2])
colorbar;

xlabel('\beta_2');
ylabel('\beta_1');
zlabel('f(\beta_1,\beta_2)');

% see:
% http://www.mathworks.com/support/solutions/en/data/1-17AF2/index.html?product=ML&solution=1-17AF2

% be drawn
new_level = -0.2;
% Get the handle to each patch object
h = findobj('type','patch');
% Create a loop to change the height of each contour
zd = get(h,'ZData');
for i = 1:length(zd)
set(h(i),'ZData',new_level*ones(length(zd{i}),1))
end
view([-50,25])

file_name=strcat(fig_dir,'bi_QNMGPpdf_cor0','.eps');
print2eps(file_name);


%for 2-dim
mu = [0 0];
SIGMA = [1 0.25; 0.25 1];
Z = mvnrnd(mu,SIGMA,n);

X=ones(n,2);
X(Z<=cut0)=0;
X(Z>cut1)=2;

count(X(:,1));
count(X(:,2));

X=X-repmat(mean(X),n,1); %centered;

betas=(-2.5:0.05:2.5)';

fs=zeros(length(betas),length(betas));

for i = 1:length(betas)
    for j = 1:length(betas)
        fs(i,j)=QNMGPpdf([betas(i); betas(j)], tau,sigma2,X);
    end
end

surfc(betas,betas,fs)
shading interp
axis([-2.5 2.5 -2.5 2.5 -0.2 0.2])
%colorbar;

xlabel('\beta_2');
ylabel('\beta_1');
zlabel('f(\beta_1,\beta_2)');

% see:
% http://www.mathworks.com/support/solutions/en/data/1-17AF2/index.html?product=ML&solution=1-17AF2

% be drawn
new_level = -0.2;
% Get the handle to each patch object
h = findobj('type','patch');
% Create a loop to change the height of each contour
zd = get(h,'ZData');
for i = 1:length(zd)
set(h(i),'ZData',new_level*ones(length(zd{i}),1))
end
view([-50,25])

file_name=strcat(fig_dir,'bi_QNMGPpdf_cor025','.eps');
print2eps(file_name);


%for 2-dim
mu = [0 0];
SIGMA = [1 0.5; 0.5 1];
Z = mvnrnd(mu,SIGMA,n);

X=ones(n,2);
X(Z<=cut0)=0;
X(Z>cut1)=2;

count(X(:,1));
count(X(:,2));

X=X-repmat(mean(X),n,1); %centered;

betas=(-2.5:0.05:2.5)';

fs=zeros(length(betas),length(betas));

for i = 1:length(betas)
    for j = 1:length(betas)
        fs(i,j)=QNMGPpdf([betas(i); betas(j)], tau,sigma2,X);
    end
end

surfc(betas,betas,fs)
shading interp
axis([-2.5 2.5 -2.5 2.5 -0.2 0.2])
%colorbar;

xlabel('\beta_2');
ylabel('\beta_1');
zlabel('f(\beta_1,\beta_2)');

% see:
% http://www.mathworks.com/support/solutions/en/data/1-17AF2/index.html?product=ML&solution=1-17AF2

% be drawn
new_level = -0.2;
% Get the handle to each patch object
h = findobj('type','patch');
% Create a loop to change the height of each contour
zd = get(h,'ZData');
for i = 1:length(zd)
set(h(i),'ZData',new_level*ones(length(zd{i}),1))
end
view([-50,25])

file_name=strcat(fig_dir,'bi_QNMGPpdf_cor05','.eps');
print2eps(file_name);

%for 2-dim
mu = [0 0];
SIGMA = [1 0.75; 0.75 1];
Z = mvnrnd(mu,SIGMA,n);

X=ones(n,2);
X(Z<=cut0)=0;
X(Z>cut1)=2;

count(X(:,1));
count(X(:,2));

X=X-repmat(mean(X),n,1); %centered;

betas=(-2.5:0.05:2.5)';

fs=zeros(length(betas),length(betas));

for i = 1:length(betas)
    for j = 1:length(betas)
        fs(i,j)=QNMGPpdf([betas(i); betas(j)], tau,sigma2,X);
    end
end

surfc(betas,betas,fs)
shading interp
axis([-2.5 2.5 -2.5 2.5 -0.2 0.2])
%colorbar;

xlabel('\beta_2');
ylabel('\beta_1');
zlabel('f(\beta_1,\beta_2)');

% see:
% http://www.mathworks.com/support/solutions/en/data/1-17AF2/index.html?product=ML&solution=1-17AF2

% be drawn
new_level = -0.2;
% Get the handle to each patch object
h = findobj('type','patch');
% Create a loop to change the height of each contour
zd = get(h,'ZData');
for i = 1:length(zd)
set(h(i),'ZData',new_level*ones(length(zd{i}),1))
end
view([-50,25])

file_name=strcat(fig_dir,'bi_QNMGPpdf_cor075','.eps');
print2eps(file_name);