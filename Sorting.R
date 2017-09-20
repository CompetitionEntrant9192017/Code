setwd('')
library(foreach)
library(Matrix)
#Sorting into acyclic array/matrix
aucVals = read.csv(file = 'IC50.csv', header = T, as.is=T)
drug.names = aucVals[2, 3:180]
fingerprints = read.csv(file = 'vougas_fingerprints.csv')
drug.names = unique(c(as.vector(aucVals[2, 3:267])))
compMatrix = array(NA, dim = c(990, 265, 265))
foreach(i=1:990)%do%
{
  sortMatrix=matrix(NA,2,265)
  sortMatrix[1,]=seq(1,265)#as.matrix(aucVals[1,3:267])
  sortMatrix[2,]=as.matrix(aucVals[(i+2),3:267])
  sortMatrix=t(sortMatrix)
  sortMatrix = sortMatrix[order(sortMatrix[,2], decreasing=FALSE), na.last=TRUE]
  sortMatrix[is.na(sortMatrix[,2]),1]=NA
  sortMatrix=t(sortMatrix)
  maxTerms = (265-sum(is.na(sortMatrix[1,])))
  foreach(j=1:maxTerms)%do%
  {
    foreach(k=j:maxTerms)%do%
    {
      drug1 = as.numeric(sortMatrix[1,j])
      drug2 = as.numeric(sortMatrix[1,k])
      if(sortMatrix[2, j]<sortMatrix[2, k])
      {
        compMatrix[i, drug1, drug2]=1
      }else{
        compMatrix[i, drug1, drug2]=0
      }
    }
  }
  foreach(j=1:265)%do%
  {
    compMatrix[i,j,j]=NA
  }
  print(i)#to keep track of progress, goes to 990
}
save.image(file = "matrixComparison.RData")
load("matrixComparison.RData")
sum = 265*265*990-sum(is.na(compMatrix[,,]))
dgDetails = which(!is.na(compMatrix),arr.ind = T)
even_indexes=seq(2,(265*265*990-sum(is.na(compMatrix[,,]))),2)
odd_indexes=seq(1,(265*265*990-sum(is.na(compMatrix[,,]))),2)
Drug1 = matrix(NA,(265*265*990-sum(is.na(compMatrix[,,]))),1)
Drug2 = matrix(NA,(265*265*990-sum(is.na(compMatrix[,,]))),1)
Drug1[odd_indexes,1]=dgDetails[odd_indexes,2]
Drug1[even_indexes,1]=dgDetails[even_indexes,3]
Drug2[odd_indexes,1]=dgDetails[odd_indexes,3]
Drug2[even_indexes,1]=dgDetails[even_indexes,2]
CellLineID = dgDetails[,1]
Comparison = c(rep(append(1, -1), trunc((265*265*990-sum(is.na(compMatrix[,,]))))/2))
if ((265*265*990-sum(is.na(compMatrix[,,])))%%2==1)
{
  Comparison = c(Comparison, 1)
}
compFinal = compMatrix[!is.na(compMatrix)]*Comparison
dgCompNew = data.frame(CellLineID, Drug1, Drug2, compFinal)
save.image("pairwise.RData")
#load("pairwise.RData")
cellfea.DF = read.csv('cell_line_features.csv',header=T,as.is=T) #Data frame and as.is keeps the characters  instead of converting into factor
gene.names = cellfea.DF[-1,1]
cell.gene = t(apply(as.matrix(cellfea.DF[-1,-1]),2,as.numeric))
colnames(cell.gene) = gene.names
cell.cancer = as.character(cellfea.DF[1,-1])
#relabel names of gene by 1-1250 corresponding to order in origin data
colnames(cell.gene) = 1:1250
gene.names = gene.names[apply(cell.gene, 2, function(x)sum(x!=0)>0)]
cell.gene = cell.gene[,apply(cell.gene, 2, function(x)sum(x!=0)>0)]

#relabel names of cancer silghtly adjusted from origin data
cancer.names = unique(cell.cancer)
cancer.names = c(cancer.names[-3],'NOT CLASSIFIED') 
cell.cancer.l = factor(cell.cancer,levels=cancer.names,labels=1:31)
save.image(file = "data.RData")
