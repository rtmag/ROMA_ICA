# MIN ICA

# X = A x S

# Matrix A is the meta sample for which
# y are the Independent components    x are the samples

# matrix S is the metagenes for which
# y are the Independent components x are the genes,cpg,dmr,etc testing

# you choose the number of IC usually 20 components is a good start for deconvolution

# the resolts include a projection file that tells you which are the most contributing genes, 
# A clinical data one clinical 
# variable can be qualitative such as  ER status (positive, negative mutation), Tumor grade, 
# or Quantitative like age

biocLite("breastCancerMAINZ")
library(breastCancerMAINZ)
data(mainz)
show(mainz)

library(MineICA)

library(JADE)
 library("Biobase")
exprs(mainz) <- t(apply(exprs(mainz),1,scale,scale=F))
colnames(exprs(mainz)) <- sampleNames(mainz)
resJade <- runICA(X=exprs(mainz),nbComp=5,method="JADE", maxit=10000 )

library(fastICA)
res = clusterFastICARuns(X=exprs(mainz),nbComp=5,alg.type="deflation",nbIt=10,funClus="hclust",method="average")
