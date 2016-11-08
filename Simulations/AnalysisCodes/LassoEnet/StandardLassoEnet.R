idnum = as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
cat("Starting run idnum = ",idnum,"\n")

#your fns start here:

library(glmnet)

#for the pc
#read_dir = 'C:/MyProjects/Admixture/Simu_data_additional/';
#idnum = 1;

#for the cluster
read_dir = '/data/zhub/Admixture/Simu_data_additional/'

nullcases=read.csv(file=paste(read_dir,'simu_null_',idnum,'.csv',sep=""),head=F,sep=",")
#[Y_binary_E0,Y_binary_E1,Y_cont_E0,Y_cont_E1,E0,E1,S_all]

results = matrix(-9, dim(nullcases)[2]-6,12)

for (rep in 1:10000){
###########################
#Y_binary_E0
###########################
#standard
logit <- glm(nullcases[,1] ~ nullcases[,5]+ nullcases[,6+rep]-1, family = "binomial")

results[rep,1]=summary(logit)$coefficients[2,4]

###########################
#Y_binary_E1
###########################
#standard
logit <- glm(nullcases[,2] ~ nullcases[,6]+ nullcases[,6+rep]-1, family = "binomial")

results[rep,4]=summary(logit)$coefficients[2,4]


###########################
#Y_cont_E0
###########################
#standard
lm <- glm(nullcases[,3] ~ nullcases[,5]+ nullcases[,6+rep]-1, family = "gaussian")

results[rep,7]=summary(lm)$coefficients[2,4]


###########################
#Y_cont_E1
###########################
#standard
logit <- glm(nullcases[,4] ~ nullcases[,6]+ nullcases[,6+rep]-1, family = "gaussian")

results[rep,10]=summary(logit)$coefficients[2,4]

}

for (rep in 1:100){
  ###########################
  #Y_binary_E0
  ###########################
   
  #lasso
  X = as.matrix(cbind(nullcases[,5], nullcases[,6+1:100+(rep-1)*100]))
  
  tmp=cv.glmnet(X,nullcases[,1],family="binomial",alpha=1,standardize=T,intercept=F)
  fit_lasso=glmnet(X,nullcases[,1],family="binomial",alpha=1,standardize=T,intercept=F)
  results[1:100+(rep-1)*100,2]=coef(fit_lasso,s=tmp$lambda.min)[1:100+2]
  
  #Enet
  tmp=cv.glmnet(X,nullcases[,1],family="binomial",alpha=0.5,standardize=T,intercept=F)
  fit=glmnet(X,nullcases[,1],family="binomial",alpha=0.5,standardize=T,intercept=F)
  results[1:100+(rep-1)*100,3]=coef(fit,s=tmp$lambda.min)[1:100+2]
  
  ###########################
  #Y_binary_E1
  ###########################
   
  #lasso
  X = as.matrix(cbind(nullcases[,6], nullcases[,6+1:100+(rep-1)*100]))
  tmp=cv.glmnet(X,nullcases[,2],family="binomial",alpha=1,standardize=T,intercept=F)
  fit=glmnet(X,nullcases[,2],family="binomial",alpha=1,standardize=T,intercept=F)
  results[1:100+(rep-1)*100,5]=coef(fit,s=tmp$lambda.min)[1:100+2]
  
  #Enet
  tmp=cv.glmnet(X,nullcases[,2],family="binomial",alpha=0.5,standardize=T,intercept=F)
  fit=glmnet(X,nullcases[,2],family="binomial",alpha=0.5,standardize=T,intercept=F)
  results[1:100+(rep-1)*100,6]=coef(fit,s=tmp$lambda.min)[1:100+2]
  
  ###########################
  #Y_cont_E0
  ###########################
  #lasso
  X = as.matrix(cbind(nullcases[,5], nullcases[,6+1:100+(rep-1)*100]))
  
  tmp=cv.glmnet(X,nullcases[,3],family="gaussian",alpha=1,standardize=T,intercept=F)
  fit_lasso=glmnet(X,nullcases[,3],family="gaussian",alpha=1,standardize=T,intercept=F)
  results[1:100+(rep-1)*100,8]=coef(fit_lasso,s=tmp$lambda.min)[1:100+2]
  
  #Enet
  tmp=cv.glmnet(X,nullcases[,3],family="gaussian",alpha=0.5,standardize=T,intercept=F)
  fit=glmnet(X,nullcases[,3],family="gaussian",alpha=0.5,standardize=T,intercept=F)
  results[1:100+(rep-1)*100,9]=coef(fit,s=tmp$lambda.min)[1:100+2]
  
  ###########################
  #Y_cont_E1
  ###########################
    
  #lasso
  X = as.matrix(cbind(nullcases[,6], nullcases[,6+1:100+(rep-1)*100]))
  tmp=cv.glmnet(X,nullcases[,4],family="gaussian",alpha=1,standardize=T,intercept=F)
  fit=glmnet(X,nullcases[,4],family="gaussian",alpha=1,standardize=T,intercept=F)
  results[1:100+(rep-1)*100,11]=coef(fit,s=tmp$lambda.min)[1:100+2]
  
  #Enet
  tmp=cv.glmnet(X,nullcases[,4],family="gaussian",alpha=0.5,standardize=T,intercept=F)
  fit=glmnet(X,nullcases[,4],family="gaussian",alpha=0.5,standardize=T,intercept=F)
  results[1:100+(rep-1)*100,12]=coef(fit,s=tmp$lambda.min)[1:100+2]
}

# function stops here
outputName=paste("task-",idnum,".RData",sep="")
outputPath=file.path("/data/zhub/Admixture/output_null",outputName)
save("results",file=outputPath)

