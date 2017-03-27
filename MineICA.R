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


source("https://bioconductor.org/biocLite.R")
biocLite("KEGG.db")
###########

### Fragment of R code from vignette source 'MineICA.Rnw'

###################################################
### code chunk number 1: morelib
###################################################
library(Biobase)
library(plyr)
library(ggplot2)
library(foreach)
library(xtable)
library(biomaRt)
library(GOstats)
library(cluster)
library(marray)
library(mclust)
library(RColorBrewer)
library(igraph)
library(Rgraphviz)
library(graph)
library(colorspace)
library(annotate)
library(scales)
library(gtools)
library(KEGG.db)


###################################################
### code chunk number 2: lib
###################################################
library(MineICA)


###################################################
### code chunk number 4: loadMainz
###################################################
## load Mainz expression data and sample annotations.
library(breastCancerMAINZ)
data(mainz)
show(mainz)
## we restrict the data to the 10,000 probe sets with the highest IQR
mainz <- selectFeatures_IQR(mainz,10000)


###################################################
### code chunk number 5: runJade
###################################################
library(JADE)
## Features are mean-centered before ICA computation
exprs(mainz) <- t(apply(exprs(mainz),1,scale,scale=FALSE))
colnames(exprs(mainz)) <- sampleNames(mainz)
## run ICA-JADE
resJade <- runICA(X=exprs(mainz), nbComp=5, method = "JADE", maxit=10000) 


###################################################
### code chunk number 6: fastica (eval = FALSE)
###################################################
## library(fastICA)
## ## Random initializations are used for each iteration of FastICA
## ## Estimates are clustered using hierarchical clustering with average linkage
## res <- clusterFastICARuns(X=exprs(mainz), nbComp=5, alg.type="deflation", nbIt=10, 
##                           funClus="hclust", method="average")


###################################################
### code chunk number 7: buildParams
###################################################
## build params
params <- buildMineICAParams(resPath="mineICA/", selCutoff=3, pvalCutoff=0.05)


###################################################
### code chunk number 8: libannot
###################################################
## load annotation package
biocLite('hgu133a.db')

###################################################
### code chunk number 9: lspackage
###################################################
#ls("package:hgu133a.db")


###################################################
### code chunk number 10: mart
###################################################
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")



###################################################
### code chunk number 13: buildIcaSet
###################################################
## Define typeID, Mainz data originate from affymetrix HG-U133a  microarray 
## and are indexed by probe sets.
## The probe sets are annotated into Gene Symbols
typeIDmainz <-  c(geneID_annotation="SYMBOL", geneID_biomart="hgnc_symbol", 
                    featureID_biomart="affy_hg_u133a")

## define the reference samples if any, here no normal sample is available
refSamplesMainz <- character(0)

resBuild <- buildIcaSet(params=params, A=data.frame(resJade$A), S=data.frame(resJade$S),
                        dat=exprs(mainz), pData=pData(mainz), refSamples=refSamplesMainz,
                        annotation="hgu133a.db", typeID= typeIDmainz,
                        chipManu = "affymetrix", mart=mart)

icaSetMainz <- resBuild$icaSet
params <- resBuild$params


###################################################
### code chunk number 14: showIcaSet
###################################################
icaSetMainz


###################################################
### code chunk number 15: defannot
###################################################
annot <- pData(icaSetMainz)


###################################################
### code chunk number 16: lookvar
###################################################
#varLabels(icaSetMainz)[1:5]
#icaSetMainz$grade[1:5]


###################################################
### code chunk number 17: lookfeatures
###################################################
#featureNames(icaSetMainz)[1:5] # probe set ids
#geneNames(icaSetMainz)[1:5] #gene symbols
#sampleNames(icaSetMainz)[1:5] 


###################################################
### code chunk number 18: lookdat (eval = FALSE)
###################################################
## head(dat(icaSetMainz)) #probe set level
## head(datByGene(icaSetMainz)) #gene level


###################################################
### code chunk number 19: lookAS (eval = FALSE)
###################################################
## A(icaSetMainz)
## S(icaSetMainz)
## SByGene(icaSetMainz)


###################################################
### code chunk number 20: lookcomp (eval = FALSE)
###################################################
## nbComp(icaSetMainz)
## compNames(icaSetMainz)
## indComp(icaSetMainz)


###################################################
### code chunk number 21: witG
###################################################
#witGenes(icaSetMainz)[1:5]
## We can for example modify the second contributing gene 
#witGenes(icaSetMainz)[2] <- "KRT16"


###################################################
### code chunk number 22: MineICA.Rnw:362-363 (eval = FALSE)
###################################################
## compNames(icaSetMainz) <- paste("IC",1:nbComp(icaSetMainz),sep="")


###################################################
### code chunk number 23: MineICA.Rnw:371-379 (eval = FALSE)
###################################################
## ## select tumor samples of grade 3
## keepSamples <- sampleNames(icaSetMainz)[icaSetMainz$grade=="3"]
## ## Subset icaSetMainz to the grade-3 samples 
## icaSetMainz[,keepSamples]
## ## Subset icaSetMainz to the grade-3 samples and the first five components
## icaSetMainz[,keepSamples,1:5]
## ## Subset icaSetMainz to the first 10 features
## icaSetMainz[featureNames(icaSetMainz)[1:10],keepSamples]


###################################################
### code chunk number 24: histProj6
###################################################
#hist(S(icaSetMainz)[,1], breaks=50, main="Distribution of feature projection on the first component", 
#     xlab="projection values")
#abline(v=c(3,-3), col="red", lty=2)


###################################################
### code chunk number 25: selectContribGenes
###################################################
## Extract the contributing genes
#contrib <- selectContrib(icaSetMainz, cutoff=3, level="genes")
## Show the first contributing genes of the first and third components
#sort(abs(contrib[[1]]),decreasing=TRUE)[1:10]
#sort(abs(contrib[[3]]),decreasing=TRUE)[1:10]


###################################################
### code chunk number 26: selectContribGenes (eval = FALSE)
###################################################
## ## One can also want to apply different cutoffs depending on the components
## ## for example using the first 4 components:
## contrib <- selectContrib(icaSetMainz[,,1:4], cutoff=c(4,4,4,3), level="genes")


###################################################
### code chunk number 27: extractComp
###################################################
## extract sample contributions and gene projections of the second component
#comp2 <- getComp(icaSetMainz, level="genes", ind=2)
## access the sample contributions 
#comp2$contrib[1:5]
## access the gene projections
#comp2$proj[1:5]


###################################################
### code chunk number 28: runAn
###################################################
## select the annotations of interest
varLabels(icaSetMainz)
# restrict the phenotype data to the variables of interest
keepVar <- c("age","er","grade")
# specify the variables that should be treated as character
icaSetMainz$er <- c("0"="ER-","1"="ER+")[as.character(icaSetMainz$er)]
icaSetMainz$grade <- as.character(icaSetMainz$grade)


###################################################
### code chunk number 29: runAn (eval = FALSE)
###################################################
## ## Run the analysis of the ICA decomposition
## # only enrichment in KEGG gene sets are tested
runAn(params=params, icaSet=icaSetMainz, writeGenesByComp = TRUE, 
       keepVar=keepVar, dbGOstats = "KEGG")
