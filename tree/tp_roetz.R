#
# Run AND save Transphylo inference
#
# Set mcmcit, the number of iterations
# Set the filename in saveRDS(...) call

# see roetzerload.R for orginal. Merged that and testforSA.R to get this file

library(TransPhylo)
library(ape)
library(treescape)

set.seed(2310)

roetree_nwk=read.tree(file="tree/roetzer2013.nwk") #equivalent I believe
roetree=read.nexus(file="roetzer2013.nex")
rfasta=read.dna("tree/roetzer2013.fasta",format="fasta")
dates=gsub("^(.*?):","",rownames(rfasta))
LastDate=sort(dates)[length(dates)] # 2010-11-01 # JS use length
LastDate=2010+(as.numeric(as.Date(LastDate)-as.Date("2010-01-01")))/365; # now in format for ptreeFromPhylo

roestart = ptreeFromPhylo(roetree, dateLastSample = LastDate)

mcmcit = 100000 # CHANGE
thinn = 10
skip = mcmcit/(200*thinn) # so we get 100 plots in second half
shift = 1 # for scaling the generation and sampling time params.
# dateT set to dateLastSample + at least one year

roerecord=inferTTree(roestart,w.shape=1.3*shift, w.scale=1/0.3, ws.shape=1.45*shift, ws.scale=1/0.3, mcmcIterations=mcmcit, thinning=thinn, 
                     startNeg=1.5, startOff.r=1, startOff.p=0.5, startPi=0.99, updateNeg=F, updatePi=F, dateT=(LastDate+1))
#plotCTree(roerecord[[1000]]$ctree)
#plotTTree(roerecord[[1000]]$ttree)
#roecontree=consTTree(roerecord[1000:1000])
#plotTTree(roecontree,w.shape=1.3, w.scale=1/0.3)
# alter from=1 as a way to implement burn-in 
# alter skipnum=1 in calling code in order to achieve 100 samples
getTTreeDistInfo3 <- function(record,skipnum=1,start_point=1) {
  ind=seq(from=start_point, by=skipnum,to=length(record))
  record=record[ind]
  matList <- lapply(1:length(record), function(x) {
    info <- extractTTree(record[[x]]$ctree)$ttree
    wiw <- cbind(info[,3],1:length(info[,1]))
    findMRCIs(wiw)$mrciDepths
  })
  return(matList)
}
# Save relevant tree depth data to file to avoid unnecessary re-runs
numSampRoetz=max(which(!is.na(extractTTree(roerecord[[1]]$ctree)$ttree[,2])))
matListRoetz = getTTreeDistInfo3(roerecord,skip,(skip+mcmcit/(thinn*2))) 
saveRDS(matListRoetz,file="TransPhyloRoetz100000_jun23.RData")
dist <- wiwTreeDist(matListRoetz, sampled=1:numSampRoetz)
wiwMDS <- dudi.pco(dist, scannf=FALSE, nf=3)
treeNames = seq(from=1, by=1, to=100)
plotGrovesD3(wiwMDS, treeNames=treeNames) # just the sim index numbers so I can see that it has burnt in OK

# Find a median:
med <- wiwMedTree(matListRoetz,sampled=1:numSampRoetz)
# medInfo <- extractTTree(record[[samp[[med$median]]]]$ctree)$ttree # JS altered: 
medInfo <- extractTTree(roerecord[[med$median]]$ctree)$ttree  # JS correct?

# Plot the blob-arrow "bean-bag" transmission tree for the median using visNetwork
# networkPlot(record[[med$median]]$ctree)
# Tidy up the names
medianTree = roerecord[[med$median]]$ctree
for (i in seq(from=1, to=numSampRoetz)){
  pos = regexpr('_', medianTree$nam[i])
  if (pos>1){medianTree$nam[i] <- substr(medianTree$nam[i],1,pos-1)}
}
networkPlot(medianTree)
plotTTree(extractTTree(medianTree),w.shape=1.3, w.scale=1/0.3)


cTreeListRoetz = getCTrees(roerecord,skip) # access eg plotCTree(cTreeList[[5]])
tTreeListRoetz <- lapply(1:length(treeNames), function(x) {extractTTree(cTreeListRoetz[[x]])})

saveRDS(cTreeListRoetz,file="tree/TransPhylo100000_cTreeListRoetz.RData")
saveRDS(tTreeListRoetz,file="tree/TransPhylo100000_tTreeListRoetz.RData")

# TRY extractPTree, hmm there is no PlotPTree()

