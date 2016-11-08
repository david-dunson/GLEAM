idnum = as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
cat("Starting run idnum = ",idnum,"\n")

#your fns start here:

library(glmnet)

#for the pc
#read_dir = 'C:/MyProjects/Admixture/Simu_data_additional/';
#idnum = 1;

#for the cluster
read_dir = '/data/zhub/Admixture/Simu_data_additional/'

nonnullcases=read.csv(file=paste(read_dir,'simu_nonnull_binary_nocov_',idnum,'.csv',sep=""),head=F,sep=",")
#[Y_binary_nocov,S,E]

results = matrix(-9, 100,15)

for (rep in 1:100){
  print(rep)
###########################
#Y_binary_E0
###########################
#standard
for(i in 0:4){
logit <- glm(nonnullcases[,rep+i*100] ~ nonnullcases[,500+rep]+ nonnullcases[,600+rep]-1, family = "binomial")

results[rep,1+i]=summary(logit)$coefficients[1,4]

#lasso
X = as.matrix(cbind(nonnullcases[,501:600],nonnullcases[,600+rep]))

tmp=cv.glmnet(X,nonnullcases[,rep+i*100],family="binomial",alpha=1,standardize=T,intercept=F)
fit_lasso=glmnet(X,nonnullcases[,rep+i*100],family="binomial",alpha=1,standardize=T,intercept=F)
results[rep,5+1+i]=coef(fit_lasso,s=tmp$lambda.min)[1+rep]

#Enet
tmp=cv.glmnet(X,nonnullcases[,rep+i*100],family="binomial",alpha=0.5,standardize=T,intercept=F)
fit=glmnet(X,nonnullcases[,rep+i*100],family="binomial",alpha=0.5,standardize=T,intercept=F)
results[rep,10+1+i]=coef(fit,s=tmp$lambda.min)[1+rep]
}
}

# function stops here
outputName=paste("task-",idnum,".RData",sep="")
outputPath=file.path("/data/zhub/Admixture/output_nonnull_binary_nocov",outputName)
save("results",file=outputPath)


