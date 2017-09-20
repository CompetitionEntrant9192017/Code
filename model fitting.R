setwd('')
load('processed.RData')
library(Matrix)
library(foreach)
library(glmnet)
library(parallel)
library(doParallel)
x = c(append(rep(append(T, F), 23), T))
Y.train = Y.train[x]
X.train = X.train[x,]
p.fac2 = rep(0.6,(265+1063*265))
p.fac2[1:265] = 0
fit = glmnet(X.train,Y.train,family='binomial',intercept=F, alpha = 0)
cvfit = cv.glmnet(X.train,Y.train,family='binomial',intercept=F, alpha = 0, nfolds = 10)
save.image(file = 'fittedNoCancerType.RData')
save.image(file = 'fittedOriginal.RData')
load("fittedOriginal.RData")
submission = read.csv('submission.6M_drug_pair.csv',head=T)
## Coefficient estimation
p.fac = rep(0.6,(265+1250*265+31*265))
p.fac[1:265] = 0
fitLinear = glmnet(X.train,Y.train,penalty.factor=p.fac,intercept=F)
save.image("LinearIntercept.RData")
cvfit = cv.glmnet(X.train,Y.train,family='multinomial',penalty.factor=p.fac,intercept=T, parallel=TRUE, nfolds =10)
save.image("fittedFinalInterceptCV.RData")
load("fittedFinalIntercept.RData")
print(fit)
print(cvfit)
load("cvfit.RData")
#cvfit = cv.glmnet(X.train,Y.train,family='multinomial',penalty.factor=p.fac,intercept=F, parallel=TRUE, nfolds = 10)
save.image("floydFittedFinal.RData")
load("Final.RData")
p.fac2 = rep(0.6,(265+1250*265+31*265))
p.fac2[1:265] = 0
fit = glmnet(X.train,Y.train,family='multinomial',penalty.factor=p.fac2,intercept=F)
save.image(file = 'fitted.RData')
load("FittedFinalIntercept.RData")
load("floydFittedElasticNet.RData")
#save.image(file = "cvfit.RData")
##prediction
Y.test.logistic = predict(fit,newx=X.test,type = "response", s=c(cvfit$lambda.min))
Y.test.linear = predict(fitLinear,newx=X.test,type = "response", s=c(fit$lambda[90]))
#Y.test.linear = predict(fitLin,newx=X.test,s=c(fit$lambda[c(100,74)]))
Y.test.class = matrix(NA,length(Y.test.logistic[,1]),1)
Y.test.class[Y.test.logistic[,1]<0.5,1] = -1
Y.test.class[Y.test.logistic[,1]>0.5,1] = 1
index.linear = is.na(Y.test.class)
Y.test.class[index.linear]=Y.test.linear[index.linear]
submission[,2] = Y.test.class[,1]
submission[,3]=Y.test.class[,1]
write.csv(submission,file='submission_ElasticLog.csv',row.names=F)
##submission[,2] = Y.test.class[,2]
##write.csv(submission,file='submission_zwz_linear_aic.csv',row.names=F)
##finalData = data.frame(Y.test.class[,1:27],Y.test)
##submission[,2] = Y.test.class[,1]
print(index.linear[300:3000])
finalData = data.frame(Y.test.class[,1])
write.csv(finalData, file = 'commonRank2.csv', row.names=F)
#finalData = data.frame(Y.test14.class[,1],Y.test12.class[,1],Y.test34.class[,1],Y.test.class[,1],Y.test)
#write.csv(finalData, file = 'confidence.csv', row.names=F)
drug.names = t(drug.names)
commonCoef = data.frame(drug.names, X)
write.csv(commonCoef, file ='commoncoef.csv', row.names = F)
##save.image('dgData.RData') 60.5%
##print(fit)
common = fit$beta[1:265,90]
d1 = common[dgComp[ index.test, 2]]
d2 = common[dgComp[ index.test, 3]]
drugCommon = data.frame(d1, d2)
write.csv(drugCommon, file = 'commonTest.csv', row.names=F)
remove(featureMatrix)
coefficientMatrix = coef(fit)[2:290176,]
featureMatrix = matrix(NA,265,1095)
featureMatrix[,1]=coefficientMatrix[1:265]
foreach(i=1:265)%do%
{
  featureMatrix[i,2:1064]=coefficientMatrix[(1063*i-797):(1063*i+265),74]
  print(i)
}
foreach(i = 1:265)%do%
{
  featureMatrix[i, 1065:1095]=coefficientMatrix[(31*i+281930):(31*i+281960), 74]
  print(i)
}
drug.names = aucVals[2, 3:267]
foreach(i =1:265)%do%{
  repeatnum = c(drug.names[i])
  if (sum(drug.names[]==repeatnum)>1){
    repeatVal = drug.names[]==repeatnum
    repeatVal[i]=F
    drug.names[repeatVal]=paste(drug.names[i], "(alternative)")
  }
}
rownames(featureMatrix)=drug.names
cancer.names = unique(cell.cancer)
colnames(featureMatrix)=c("Common Ranking", gene.names, cancer.names)
featureVectors = data.frame(featureMatrix)
write.csv(featureVectors, file='featureCoefficientsOriginal.csv', row.names=T)