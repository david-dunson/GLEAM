function cnts2=nonnull_2loci_cnt(idx, locus1,locus2)

cnts=zeros(1,6);

cnts(1)=length(idx)- ismember(locus1,idx)- ismember(locus2,idx);
if ismember(locus1,idx)==1 & ismember(locus2,idx)~=1
cnts(2)=1;
elseif ismember(locus1,idx)~=1 & ismember(locus2,idx)==1
cnts(3)=1;
elseif ismember(locus1,idx)==1 & ismember(locus2,idx)==1
cnts(4)=1;            
end

region25=[3:8,10:24,26:46];
region75=[66:74,76,78:102];

cnts(5)=sum(ismember(idx,region25));
cnts(6)=sum(ismember(idx,region75));

cnts(1)=cnts(1)-cnts(5)-cnts(6);

cnts2=[cnts(5)/100,cnts(6)/100,cnts(1)/100,cnts(2:4)];