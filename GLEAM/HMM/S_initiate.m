%PURPOSE: Initiate S matrix
%INPUT: 
%Q 3*3*I*J matrix state transition matrix
%Q_0 I*3 matrix of initial states mass matrix 
%chrom: J*1 chromosome number
%rand_seed: random seed
%OUTPUT: 
%S: I*J state matrix 

function [S] = S_initiate(Q, Q_0,chrom,rand_seed)

RandState=rand_seed;
rand('state',RandState); % set arbitrary seed for uniform draws
randn('state',RandState); % set arbitrary seed for normal draws


[tmp1,tmp2, I, J]=size(Q);

S=zeros(I,J);

for i=1:I
    
    S_i=zeros(1,J);
    
    for k=1:23
        idx=(chrom==k); %for chromosome k;
        Q_ik=Q(:,:,i,idx);
        %S_ikj=sample(Q_0(i,:),1); % initial value at position zero
        S_ikj=my_sample_1(Q_0(i,:));
        
        S_ik=zeros(1,sum(idx));
        
        for j=1:sum(idx)
            %S_ikj=sample(Q_ik(S_ikj,:,1,j),1);
            S_ikj=my_sample_1(Q_ik(S_ikj,:,1,j));
            S_ik(j)=S_ikj-1;
            %S_i=[S_i S_ikj-1];
        end
        
        S_i(idx)=S_ik;
    end
    
    S(i,:)=S_i;
end
