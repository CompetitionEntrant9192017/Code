setwd('')
library(Matrix)
library(foreach)
library(parallel)
library(doParallel)
library(Rfast)
cellfea.DF = read.csv('cell_line_features.csv',header=T,as.is=T) #Data frame and as.is keeps the characters  instead of converting into factor
cell.gene = t(apply(as.matrix(cellfea.DF[-1,-1]),2,as.numeric))
gene.names = cellfea.DF[-1,1]
colnames(cell.gene) = gene.names
cell.cancer = as.character(cellfea.DF[1,-1])


colnames(cell.gene) = 1:1250

#relabel names of cancer silghtly adjusted from origin data
cancer.names = unique(cell.cancer)
cancer.names = c(cancer.names[-3],'NOT CLASSIFIED') 
cell.cancer.l = factor(cell.cancer,levels=cancer.names,labels=1:31)

dgComp = read.csv('train3M_drug_pairMAIN.csv',header=T)
drug.names = unique(c(as.vector(dgComp$Drug1),as.vector(dgComp$Drug2)))
ll = factor(c(as.vector(dgComp$Drug1),as.vector(dgComp$Drug2)),levels=drug.names,labels=1:265)
dgComp$Drug1 = ll[1:3000000]
dgComp$Drug2 = ll[3000001:6000000]
cell.names = rownames(cell.gene)
sample.names = factor(as.character(dgComp[,1]),levels=cell.names,labels=1:990)
dgComp[,1] = as.numeric(sample.names)
floydMatrices = array(NA, dim = c(990, 265, 265))
dgComp2=dgComp[500001:1000000,]
dgComp3=dgComp[1000001:1500000,]
dgComp4=dgComp[1500001:2000000,]
dgComp5=dgComp[2000001:2500000,]
dgComp6=dgComp[2500001:3000000,]
foreach(i=1:500000)%do%{
#length(dgComp[,1])
  if(is.na(dgComp[i,4])){}
  else if (dgComp[i, 4]==1){floydMatrices[dgComp[i, 1], dgComp[i, 2], dgComp[i, 3]]=1}
  else if (dgComp[i, 4]==-1){floydMatrices[dgComp[i, 1], dgComp[i, 3], dgComp[i, 2]]=1}
  else if (dgComp[i, 4]==0){
    floydMatrices[dgComp[i, 1], dgComp[i, 3], dgComp[i, 2]]=0
    floydMatrices[dgComp[i, 1], dgComp[i, 2], dgComp[i, 3]]=0
  }
}
print("1")
foreach(i=1:500000)%do%{
  if(is.na(dgComp2[i,4])){}
  else if (dgComp2[i, 4]==1){floydMatrices[dgComp2[i, 1], dgComp2[i, 2], dgComp2[i, 3]]=1}
  else if (dgComp2[i, 4]==-1){floydMatrices[dgComp2[i, 1], dgComp2[i, 3], dgComp2[i, 2]]=1}
  else if (dgComp2[i, 4]==0){
    floydMatrices[dgComp2[i, 1], dgComp2[i, 3], dgComp2[i, 2]]=0
    floydMatrices[dgComp2[i, 1], dgComp2[i, 2], dgComp2[i, 3]]=0
  }
}
print("2")
foreach(i=1:500000)%do%{
  if(is.na(dgComp3[i,4])){}
  else if (dgComp3[i, 4]==1){floydMatrices[dgComp3[i, 1], dgComp3[i, 2], dgComp3[i, 3]]=1}
  else if (dgComp3[i, 4]==-1){floydMatrices[dgComp3[i, 1], dgComp3[i, 3], dgComp3[i, 2]]=1}
  else if (dgComp3[i, 4]==0){
    floydMatrices[dgComp3[i, 1], dgComp3[i, 3], dgComp3[i, 2]]=0
    floydMatrices[dgComp3[i, 1], dgComp3[i, 2], dgComp3[i, 3]]=0
  }
}
print("3")
foreach(i=1:500000)%do%{
  if(is.na(dgComp4[i,4])){}
  else if (dgComp4[i, 4]==1){floydMatrices[dgComp4[i, 1], dgComp4[i, 2], dgComp4[i, 3]]=1}
  else if (dgComp4[i, 4]==-1){floydMatrices[dgComp4[i, 1], dgComp4[i, 3], dgComp4[i, 2]]=1}
  else if (dgComp4[i, 4]==0){
    floydMatrices[dgComp4[i, 1], dgComp4[i, 3], dgComp4[i, 2]]=0
    floydMatrices[dgComp4[i, 1], dgComp4[i, 2], dgComp4[i, 3]]=0
  }
}
print("4")
foreach(i=1:500000)%do%{
  if(is.na(dgComp5[i,4])){}
  else if (dgComp5[i, 4]==1){floydMatrices[dgComp5[i, 1], dgComp5[i, 2], dgComp5[i, 3]]=1}
  else if (dgComp5[i, 4]==-1){floydMatrices[dgComp5[i, 1], dgComp5[i, 3], dgComp5[i, 2]]=1}
  else if (dgComp5[i, 4]==0){
    floydMatrices[dgComp5[i, 1], dgComp5[i, 3], dgComp5[i, 2]]=0
    floydMatrices[dgComp5[i, 1], dgComp5[i, 2], dgComp5[i, 3]]=0
  }
}
print("5")
foreach(i=1:500000)%do%{
  if(is.na(dgComp6[i,4])){}
  else if (dgComp6[i, 4]==1){floydMatrices[dgComp6[i, 1], dgComp6[i, 2], dgComp6[i, 3]]=1}
  else if (dgComp6[i, 4]==-1){floydMatrices[dgComp6[i, 1], dgComp6[i, 3], dgComp6[i, 2]]=1}
  else if (dgComp6[i, 4]==0){
    floydMatrices[dgComp6[i, 1], dgComp6[i, 3], dgComp6[i, 2]]=0
    floydMatrices[dgComp6[i, 1], dgComp6[i, 2], dgComp6[i, 3]]=0
  }
}
#May be changed due to Rfast being updated with the floyd algorithm after I corrected them
floyd_helper<-function(x,na_or_0 = NA){
  i4_huge<-2147483647
  x[is.na(x)]<-i4_huge
  y<-Rfast::floyd(x)
  y[is.na(y)]<-na_or_0
  y
}
foreach(i = 1:length(floydMatrices[,1,1]))%do%{
    print(i)
    floydMatrices[i,,]=floyd_helper(floydMatrices[i,,])
}
floydMatrices[floydMatrices==0]=NA
Comparison = c(rep(append(1, -1), trunc((265*265*990-sum(is.na(floydMatrices[,,]))))/2))
if ((265*265*990-sum(is.na(floydMatrices[,,])))%%2==1)
{
  Comparison = c(Comparison, 1)
}
length(Comparison)
(265*265*990-sum(is.na(floydMatrices[,,])))
dgDetails = which(!is.na(floydMatrices),arr.ind = T)
even_indexes=seq(2,(265*265*990-sum(is.na(floydMatrices[,,]))),2)
odd_indexes=seq(1,(265*265*990-sum(is.na(floydMatrices[,,]))),2)
Drug1 = matrix(NA,(265*265*990-sum(is.na(floydMatrices[,,]))),1)
Drug2 = matrix(NA,(265*265*990-sum(is.na(floydMatrices[,,]))),1)
Drug1[odd_indexes,1]=dgDetails[odd_indexes,2]
Drug1[even_indexes,1]=dgDetails[even_indexes,3]
Drug2[odd_indexes,1]=dgDetails[odd_indexes,3]
Drug2[even_indexes,1]=dgDetails[even_indexes,2]
CellLineID = dgDetails[,1]
dgCompNew = data.frame(CellLineID, Drug1, Drug2, Comparison)
gene.names = cellfea.DF[-1,1]
colnames(cell.gene) = gene.names
cell.cancer = as.character(cellfea.DF[1,-1])
colnames(cell.gene) = 1:1250

#relabel names of cancer silghtly adjusted from origin data
cancer.names = unique(cell.cancer)
cancer.names = c(cancer.names[-3],'NOT CLASSIFIED') 
cell.cancer.l = factor(cell.cancer,levels=cancer.names,labels=1:31)
save.image(file = "completeGraph0.RData")
load("completeGraph0.RData")
index.train = rep(append(rep(T,1),rep(F,10)),1706043)
#index.test = (is.na(dgComp[,4]) ==T)