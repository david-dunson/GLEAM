%PURPOSE:
%Forward Filtering Backward Sampling only; calculate log likelihood
%INPUT:
%P: 3*3*J observation probability mass function
%Q_cdn: 3*3*1*3 conditional state transition matrix
%Q_0: 3*1 initial states probability
%R: 1*J recombination counts
%X: 1*J observations in 0,1,2
%OUTPUT:
%S: 1*J states
%S_0: initial state;
%logL: 1*J log likelihood

%P=P(:,:,idx);Q_cdn=Q_cdn(:,:,i,:);Q_0=Q_0(i,:)'; R=R(i,idx); SNPs2=obj.SNPs;X=SNPs2(i,idx);
function [S,S_0, L]=FFBS2(P,Q_cdn,Q_0,R,X)

J=length(X);

Q_F=zeros(3,3,J); %Q_F: 3*3*J forward matrix
S=zeros(1,J);
L=zeros(1,J); %likelihodd


%forward filtering
Q_F_j_marg=Q_0;
for j=1:J   
    %Q_F(:,:,j) = repmat(Q_F_j_marg,1,3) .* Q_cdn(:,:,1,R(j)+1) .* repmat(P(:,X(j)+1,j)',3,1);

    Q_F(:,:,j)=FWkernal(Q_F_j_marg,Q_cdn(:,:,1,R(j)+1),P(:,X(j)+1,j)); % not be scaled, but o.k.
    Q_F_j_marg(1,1)=  Q_F(1,1,j)+ Q_F(2,1,j)+ Q_F(3,1,j); %faster
    Q_F_j_marg(2,1)=  Q_F(1,2,j)+ Q_F(2,2,j)+ Q_F(3,2,j);
    Q_F_j_marg(3,1)=  Q_F(1,3,j)+ Q_F(2,3,j)+ Q_F(3,3,j);
     
    L(j) = Q_F_j_marg(1)+Q_F_j_marg(2)+Q_F_j_marg(3);%faster
    Q_F_j_marg=Q_F_j_marg/L(j); % MUST be scaled, otherwise numberic round errors;
    
end


%backward sampling
S_next=my_sample_1(Q_F_j_marg);% faster version

S(J)=S_next-1;
for j=J:-1:2
    S_next=my_sample_1(Q_F(:,S_next,j));
    S(j-1)=S_next-1;
end;

S_0 = my_sample_1(Q_F(:,S_next,1))-1;



