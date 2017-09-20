setwd('')
load('data.RData')
library(Matrix)
## Sparse design matrix
# 1st (drug baseline)
ii = as.numeric(dgComp[,2])
jj = as.numeric(dgComp[,3])
X1 = sparseMatrix(rep(1:23828227,2),c(ii,jj),x=rep(c(1,-1),each=23828227),dims=c(23828227,265))
# 2nd (gene mutation 3,000,000 * (1250 * 265))
cell.gene.s = Matrix(cell.gene,sparse=T)
index.cell.dg1 = (as.numeric(dgComp[,2]) - 1) * 990 + dgComp[,1]
index.cell.dg2 = (as.numeric(dgComp[,3]) - 1) * 990 + dgComp[,1]
index.cell.1 = sparseMatrix(i=1:23828227,j=index.cell.dg1, x=rep(1,23828227),dims=c(23828227,990*265))
index.cell.2 = sparseMatrix(i=1:23828227,j=index.cell.dg2, x=rep(-1,23828227),dims=c(23828227,990*265))
index.cell = index.cell.1 + index.cell.2 #3,000,000 * (990*265)

cell.gene.block = bdiag(rep(list(A=cell.gene.s),265)) #(990*265) * (1250*265)
X2 = index.cell %*% cell.gene.block

# 3rd (cancer type 3,000,000 * (265*31))
index.cancer = rep(1,23828227)
for(ii in 1:23828227){
	index.cancer[ii] = cell.cancer.l[dgComp[ii,1]]
}

index.cancer.dg1 = (as.numeric(dgComp[,2]) - 1)*31 + index.cancer
index.cancer.dg2 = (as.numeric(dgComp[,3]) - 1)*31 + index.cancer
X3 = sparseMatrix(1:23828227,index.cancer.dg1,x=rep(1,23828227),dims=c(23828227,31*265)) - sparseMatrix(1:23828227, index.cancer.dg2,x=rep(1,23828227),dims=c(23828227,31*265))
# Train sample index and Test sample index
Y = dgComp[,4]
index.train = rep(T,23828227)
index.0 = (dgComp[,4] != 0)
index.notna = (is.na(dgComp[,4])==F) #not NA index
index.train = index.train & index.0 & index.notna
Xonetwo = cbind(X1[index.train,], X2[index.train,])
remove(X1)
remove(X2)
X.train = cbind(Xonetwo, X3[index.train,])
X.train = X1[index.train,]
remove(X3)
remove(Xonetwo)
X.test = cbind(X1[index.test,],X2[index.test,],X3[index.test,])
Y.train = Y[index.train]
save.image(file="processed.RData")
