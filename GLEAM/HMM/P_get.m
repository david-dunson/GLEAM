%PURPOSE:
%calcualte P matrix
%INPUT:
%pa: 1*J population A variant allele frequncies
%pb: 1*J population B variant allele frequncies
%OUTPUT: 
%P 3*3*J matrix observation probability mass matrix

% pa=(afrn_var/60)';pb=(caucn_var/60)';
function [P] = P_get(pa,pb)

J=length(pa);

one_pa=1-pa;
one_pb=1-pb;

P=zeros(3,3,J);

P(1,1,:)=one_pb.^2;
P(1,2,:)=2.*pb.*one_pb;
P(1,3,:)=pb.^2;
P(2,1,:)=one_pa.*one_pb;
P(2,2,:)=pa.*one_pb+pb.*one_pa;
P(2,3,:)=pa.*pb;
P(3,1,:)=one_pa.^2;
P(3,2,:)=2.*pa.*one_pa;
P(3,3,:)=pa.^2;
