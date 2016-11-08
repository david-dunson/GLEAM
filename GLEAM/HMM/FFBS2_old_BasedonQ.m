%PURPOSE:
%Forward Filtering Backward Sampling only; calculate log likelihood
%INPUT:
%P: 3*3*J observation probability mass function
%Q: 3*3*J state transition matrix
%Q_0: 3*1 initial states probability
%X: 1*J observations in 0,1,2
%OUTPUT:
%S: 1*J states
%logL: 1*J log likelihood


function [S, logL]=FFBS2(P,Q,Q_0,X)

J=length(X);

Q_F=zeros(3,3,J); %Q_F: 3*3*J forward matrix
%Q_F_marg=zeros(3,J);
S=zeros(1,J);
logL=zeros(1,J);

%Q_F=[]; %Q_F: 3*3*J forward matrix
%Q_F_marg=[]; %3*J
%S=[]; % 1*J
%logL=[]; % 1*J

%forward filtering
Q_F_j_marg_prev=Q_0;
j=1;
for j=1:J
    tmp = repmat(Q_F_j_marg_prev,1,3) .* Q(:,:,j) .* repmat(P(:,X(j)+1,j)',3,1);
    L_j=sum(sum(tmp));
    Q_F_j = tmp ./ L_j;
    Q_F_j_marg =  sum(Q_F_j,1); %column sum
    Q_F(:,:,j)=Q_F_j;
    %Q_F=cat(3,Q_F,Q_F_j);
    %Q_F_marg(:,j)=Q_F_j_marg';
    %Q_F_marg=[Q_F_marg Q_F_j_marg'];
    logL(j)=log(L_j);
    %logL=[logL log(L_j)];
    Q_F_j_marg_prev=Q_F_j_marg';
end


%backward sampling
%S_next=randsample(3,1,true,Q_F_marg(:,J));
%S_next=randsample(3,1,true,Q_F_j_marg');
S_next=sample(Q_F_j_marg',1);

S(J)=S_next-1;
%S=S_next-1;
i=1;
for i=1:(J-1)
    j=J-i+1;
    %S_next=randsample(3,1,true,Q_F(:,S_next,j));
    S_next=sample(Q_F(:,S_next,j),1);
    S(j-1)=S_next-1;
    %S=[S_next-1 S];
end;
