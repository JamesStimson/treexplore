#
# Run AND save Transphylo inference
#
# Set mcmcit, the number of iterations
# Set the filename in saveRDS(...) call

# NEED TO RUN transphylo_extras.R

library(TransPhylo)
library(ape)
library(treescape)

# JS 30/05 set to 100,000 first. Up to 1 million when all looks OK.

set.seed(2310)
sclade=read.tree(file="tree/sclade-newlabs.nwk") # plot(sclade) gives basic phylo tree
startpoint=ptreeFromPhylo(sclade, dateLastSample = 2012.5) # is this accurate enough? 2012.53

mcmcit = 100000
thinn = 10
skip = mcmcit/(200*thinn) # so we get 100 plots in second half
shift = 1 # for scaling the generation and sampling time params.
# dateT set to dateLastSample + at least one year
# thinning - change from 10? Not for now, use skipnum below
#record=inferTTree(startpoint,w.shape=1.3, w.scale=1/0.3, ws.shape=1.45,ws.scale=1/0.3, mcmcIterations=mcmcit, thinning=thinn, startNeg=1.5, startOff.r=1, startOff.p=0.5, startPi=0.99, updateNeg=F, updatePi=F, dateT=2013.6)
record=inferTTree(startpoint,w.shape=1.3*shift, w.scale=1/0.3, ws.shape=1.45*shift, ws.scale=1/0.3, mcmcIterations=mcmcit, thinning=thinn, 
                  startNeg=1.5, startOff.r=1, startOff.p=0.5, startPi=0.99, updateNeg=F, updatePi=F, dateT=2013.6)

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
numSamp=max(which(!is.na(extractTTree(record[[1]]$ctree)$ttree[,2])))
matList = getTTreeDistInfo3(record,skip,(skip+mcmcit/(thinn*2))) #added by JS 16/05, amended 14/06
saveRDS(matList,file="TransPhylo100000Pi99_jun23.RData")
dist <- wiwTreeDist(matList, sampled=1:numSamp)
wiwMDS <- dudi.pco(dist, scannf=FALSE, nf=3)
treeNames = seq(from=1, by=1, to=100)
plotGrovesD3(wiwMDS, treeNames=treeNames) # just the sim index numbers so I can see that it has burnt in OK

# Find a median:
med <- wiwMedTree(matList,sampled=1:numSamp)
# medInfo <- extractTTree(record[[samp[[med$median]]]]$ctree)$ttree # JS altered: 
medInfo <- extractTTree(record[[med$median]]$ctree)$ttree  
# Plot the blob-arrow "bean-bag" transmission tree for the median using visNetwork
# Tidy up the names
medianTree = record[[med$median]]$ctree
for (i in seq(from=1, to=numSamp)){
  pos = regexpr('_', medianTree$nam[i])
  if (pos>1){medianTree$nam[i] <- substr(medianTree$nam[i],5,pos-1)}
}
networkTPlot(extractTTree(medianTree))
plotTTree(extractTTree(medianTree),w.shape=1.3, w.scale=1/0.3)
saveRDS(medianTree,file="tree/MedianTPv_Trunc.RData")

# Extract all the trees for all the individual graphs
getCTrees <- function(record,skipnum=1) {
  ind=seq(from=skip+mcmcit/(thinn*2), by=skipnum,to=length(record))
  record=record[ind]
  cTreeList <- lapply(1:length(record), function(x) {
    for (i in seq(from=1, to=numSamp)){
      pos = regexpr('_', record[[x]]$ctree$nam[i])
      if (pos>1){record[[x]]$ctree$nam[i] <- substr(record[[x]]$ctree$nam[i],1,pos-1)}
    }
    record[[x]]$ctree
  })
  return(cTreeList)
}
cTreeList = getCTrees(record,skip) # access eg plotCTree(cTreeList[[5]])
tTreeList <- lapply(1:length(treeNames), function(x) {extractTTree(cTreeList[[x]])})

saveRDS(cTreeList,file="tree/TransPhylo100000_cTreeList_jun23.RData")
saveRDS(tTreeList,file="tree/TransPhylo100000_tTreeList_jun23.RData")


