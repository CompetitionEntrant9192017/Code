setwd('')
library(foreach)
library(Matrix)
#Sorting into acyclic array/matrix
aucVals = read.csv(file = 'IC50.csv', header = T, as.is=T)
fingerprints = read.csv(file = 'vougas_fingerprintsFixed.csv')
smiles = read.csv(file = 'vougas_smiles.csv')
drug.names = unique(c(as.vector(fingerprints$Name)))
ll = factor(aucVals[2, 3:267],levels=drug.names,labels=1:212)
ll = c(0, 0, ll)
aucVals[2,]=ll
aucVals = aucVals[,!is.na(ll)]
x = c(0)
dupes = c(0)
foreach(i = 1:212)%do%
{
  if(sum(aucVals[2,]==i)>1)
    dupes=c(dupes, i)
}
dupes = dupes[2:length(dupes)]
dupeCoord = c(0)
dupeCoordInitial = c(0)
foreach(i = 1:9)%do%
{
  drug.names = c(drug.names, paste(drug.names[dupes[i]], "(alt)"))
  dupeCoord = c(dupeCoord, which(aucVals[2,]==dupes[i])[2])
  dupeCoordInitial = c(dupeCoordInitial, which(aucVals[2,]==dupes[i])[1])
}
dupeCoord = dupeCoord[2:10]
dupeCoordInitial = dupeCoordInitial[2:10]
aucVals[2, dupeCoord] = seq(213,221)
compMatrix = array(NA, dim = c(990, 221, 221))
foreach(i=1:990)%do%
{
  sortMatrix=matrix(NA,2,221)
  sortMatrix[1,]=seq(1,221)#as.matrix(aucVals[1,3:267])
  sortMatrix[2,]=as.matrix(aucVals[(i+2),3:223])
  sortMatrix=t(sortMatrix)
  sortMatrix = sortMatrix[order(sortMatrix[,2], decreasing=FALSE), na.last=TRUE]
  sortMatrix[is.na(sortMatrix[,2]),1]=NA
  sortMatrix=t(sortMatrix)
  maxTerms = (221-sum(is.na(sortMatrix[1,])))
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
  foreach(j=1:221)%do%
  {
    compMatrix[i,j,j]=NA
  }
  print(i)#to keep track of progress, goes to 990
}
sum(compMatrix[!is.na(compMatrix)]==0)
save.image(file = "drugFeatureMatrixComparison.RData")
xx = factor(fingerprints[,2],levels=drug.names,labels=1:221)
fingerprints[,2]=xx
fingerprints[213:221,2] = seq(213, 221)
fingerprints[213:221,3:883]=fingerprints[dupes,3:883]
fingerprints[,884]=c(rep(0, 212), rep(1, 9))
sum = 221*221*990-sum(is.na(compMatrix[,,]))
dgDetails = which(!is.na(compMatrix),arr.ind = T)
even_indexes=seq(2,(221*221*990-sum(is.na(compMatrix[,,]))),2)
odd_indexes=seq(1,(221*221*990-sum(is.na(compMatrix[,,]))),2)
Drug1 = matrix(NA,(221*221*990-sum(is.na(compMatrix[,,]))),1)
Drug2 = matrix(NA,(221*221*990-sum(is.na(compMatrix[,,]))),1)
Drug1[odd_indexes,1]=dgDetails[odd_indexes,2]
Drug1[even_indexes,1]=dgDetails[even_indexes,3]
Drug2[odd_indexes,1]=dgDetails[odd_indexes,3]
Drug2[even_indexes,1]=dgDetails[even_indexes,2]
CellLineID = dgDetails[,1]
Comparison = c(rep(append(1, -1), trunc((221*221*990-sum(is.na(compMatrix[,,]))))/2))
if ((221*221*990-sum(is.na(compMatrix[,,])))%%2==1)
{
  Comparison = c(Comparison, 1)
}
compFinal = compMatrix[!is.na(compMatrix)]*Comparison
dgCompNew = data.frame(CellLineID, Drug1, Drug2, compFinal)
save.image("drugFeaturePairwise.RData")
load("drugFeaturepairwise.RData")
cellfea.DF = read.csv('cell_line_features.csv',header=T,as.is=T) #Data frame and as.is keeps the characters  instead of converting into factor
gene.names = cellfea.DF[-1,1]
cell.gene = t(apply(as.matrix(cellfea.DF[-1,-1]),2,as.numeric))
colnames(cell.gene) = gene.names
cell.cancer = as.character(cellfea.DF[1,-1])
#relabel names of gene by 1-1250 corresponding to order in origin data
colnames(cell.gene) = 1:1250
cell.gene = cell.gene[,apply(cell.gene, 2, function(x)sum(x!=0)>0)]
#relabel names of cancer silghtly adjusted from origin data
cancer.names = unique(cell.cancer)
cancer.names = c(cancer.names[-3],'NOT CLASSIFIED') 
cell.cancer.l = factor(cell.cancer,levels=cancer.names,labels=1:31)
fingerprints[213:221, 1]= seq(266:274)
fingerprints = fingerprints[,apply(fingerprints, 2, function(x)sum(x!=0)>0)]
save.image(file = "drugFeatureData.RData")
