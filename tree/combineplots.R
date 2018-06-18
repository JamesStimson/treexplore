# Run testforSA.R and phybreak_vegard first time.

# something about indicating number sampled/unsampled?

# comment in/out as applicable when updating. scotti and beastlier are quick, as they are just loading from file.
# pick correct filenames for matList and pbMatList
# warning: hard-coded number of sims may need review

library(TransPhylo)
library(ape)
#####library(treescape)
library(treespace)
library(phybreak)
library(shinyBS)
library(shinythemes)
library(scatterD3)
library(markdown) # for the navbar code
library(visNetwork)

source("tree/transphylo_extras.R")

#source("testforSA.R")
#source("phybreak_vegard.R")
source("tree/scotti_load.R")
source("tree/beastlier_load.R")
source("tree/scotti_roetz.R")
source("tree/beastlier_roetz.R")

#source("tp_roetz.R")

nTransPhylo = 100
nMcmc = 100#nrow(mcmc)/step
#nScotti set in scotti load
nBeastlier = 100

numSamp = 13
numSampRoetz = 86

matList = readRDS("tree/TransPhylo100000Pi99DoubleShape.RData")
pbMatList = readRDS("tree/PhyBreak20000DoubleShape.RData") #20000 rather than 10000 moves to transphylo, away from beastlier

cMatList = c(matList,pbMatList,scMatList,bMatList)

names(cMatList)[1:nTransPhylo] <- paste0("TransPhylo",1:nTransPhylo)
names(cMatList)[(nTransPhylo+1):(nTransPhylo+nMcmc)] <- paste0("PhyBreak",1:nMcmc)
names(cMatList)[(nTransPhylo+nMcmc+1):(nTransPhylo+nMcmc+nScotti)] <- paste0("SCOTTI",1:nScotti)
names(cMatList)[(nTransPhylo+nMcmc+nScotti+1):(nTransPhylo+nMcmc+nScotti+nBeastlier)] <- paste0("Beastlier",1:nBeastlier)

Dtype <- c(rep("TransPhylo",nTransPhylo),rep("PhyBreak",nMcmc),rep("SCOTTI",nScotti),rep("Beastlier",nBeastlier)) # see https://cran.r-project.org/web/packages/treescape/vignettes/TransmissionTreesVignette.html
colours <- c(rep("#4682B4",1),rep("#DC143C",1),rep("#228B22",1),rep("#EE9900",1))

newdist <- wiwTreeDist(cMatList, sampled=1:numSamp)
newwiwMDS <- dudi.pco(newdist, scannf=FALSE, nf=3) # why nf=3? 3->5 seemed to make no difference

treeNames = c(seq(from=1, to=nTransPhylo), seq(from=1, to=nMcmc), seq(from=1, to=nScotti), seq(from=1, to=nBeastlier))

#plotGrovesD3(newwiwMDS, groups=Dtype, point_size=32, colors=colours, col_lab="Method", ellipses=TRUE, ellipses_level=0.95, treeNames=treeNames) # names to sanity check plot



# PLUS Roetzer data
matListRoetz = readRDS("tree/TransPhyloRoetz100000.RData")
pbMatListRoetz = readRDS("tree/PhyBreak4000Roetz.RData")
cMatListRoetz = c(matListRoetz, pbMatListRoetz, scMatListRoetz, bMatListRoetz)  # INCOMPLETE

DtypeRoetz <- c(rep("TransPhylo",nTransPhylo),rep("PhyBreak",nMcmc),rep("SCOTTI",nScottiRoetz),rep("Beastlier",nBeastlier)) # see https://cran.r-project.org/web/packages/treescape/vignettes/TransmissionTreesVignette.html
treeNamesRoetz = c(seq(from=1, to=nTransPhylo), seq(from=1, to=nMcmc), seq(from=1, to=nScottiRoetz), seq(from=1, to=nBeastlier))

newdistRoetz <- wiwTreeDist(cMatListRoetz, sampled=1:numSampRoetz) # 
newwiwRoetz <- dudi.pco(newdistRoetz, scannf=FALSE, nf=3)
#plotGrovesD3(newwiwRoetz, groups=DtypeRoetz, point_size=32, colors=colours, col_lab="Method", ellipses=TRUE, ellipses_level=0.95, treeNames=treeNamesRoetz)

# these rely on individual runs for transphylo and phybreak. Should save to file.
cTreeList = readRDS("tree/TransPhylo100000_cTreeList.RData")
cTreeListRoetz = readRDS("tree/TransPhylo100000_cTreeListRoetz.RData")

tTreeList = readRDS("tree/TransPhylo100000_tTreeList.RData")
pbTTreeList = readRDS("tree/PhyBreak20000_TTreeList.RData")
cTTreeList = c(tTreeList, pbTTreeList, scTTreeList, bTTreeList)

tTreeListRoetz = readRDS("tree/TransPhylo100000_tTreeListRoetz.RData")
pbTTreeListRoetz = readRDS("tree/PhyBreak2000_TTreeListRoetz.RData")
cTTreeListRoetz = c(tTreeListRoetz, pbTTreeListRoetz, scTTreeListRoetz, bTTreeListRoetz)

MedianTree1 = extractTTree(readRDS("tree/MedianTPv_Trunc.RData"))
MedianTree2 = readRDS("tree/MedianPBv_Trunc.RData")
MedianTree3 = readRDS("tree/MedianSCv_Trunc.RData")
MedianTree4 = readRDS("tree/MedianBv_Trunc.RData")

myMCMCstate = readRDS("tree/MCMCstate20000.RData")
roetzMCMCstate = readRDS("tree/MCMCstate2000Roetz.RData")

butterfly = readRDS("tree/Caterpillar.RData")


