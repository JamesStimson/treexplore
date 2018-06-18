# phybreak attempt JS May 8
library(phybreak)
set.seed(2310)
# use phybreak() to create a phybreak object
# input sequences and sampling times
seq_data = read.phyDat("tree/sclade_sequence_alignment.fa", format="fasta", type="DNA")
# assume times are years from start so scotti data/365.0
#times = c(5.5, 0.5, 99.5, 267.5, 736.5, 861.5, 963.5, 1244.5, 1780.5, 1813.5, 1838.5, 2129.5, 2401.5) 
time_years = c(0.015, 0.001, 0.273, 0.733, 2.018, 2.36, 2.64, 3.41, 4.878, 4.968, 5.037, 5.834, 6.579)
first_date = 2005.956
names(time_years) = c("Case27_2005-12-20", "Case28_2005-12-15", "Case29_2006-03-24", "Case30_2006-09-08", "Case37_2007-12-21", 
                      "Case38_2008-04-24", "Case40_2008-08-04", "Case47_2009-05-12", "Case52_2010-10-30", "Case54_2010-12-02", 
                      "Case55_2010-12-27", "Case68_2011-10-14", "Case81_2012-07-12")
short_names = c("Case27", "Case28", "Case29", "Case30", "Case37", "Case38", "Case40", "Case47", "Case52", "Case54", "Case55", "Case68", "Case81")
very_short_names= c("27", "28", "29", "30", "37", "38", "40", "47", "52", "54", "55", "68", "81")
# Gamma: The mean and variance are E(X) = a*s and Var(X) = a*s^2, where a is shape (gen 1.3) and s is scale (gen 1/0.3) 
# changing slope from 1 to 3 does make a diff


shift = 1 # 1 is deafault

myMCMCstate = phybreak(data=seq_data, times=time_years, mu=1e-04, 
                       gen.shape=shift*1.3, gen.mean=shift*1.3/0.3, sample.shape=shift*1.45, sample.mean=shift*1.45/0.3, 
                       wh.model=3, wh.slope=3, est.gen.mean=TRUE, est.sample.mean=TRUE, 
                       prior.mean.gen.mean=shift*1.3/0.3, prior.mean.gen.sd=sqrt(shift*1.3)/0.3, 
                       prior.mean.sample.mean=shift*1.45/0.3, prior.mean.sample.sd=sqrt(shift*1.45)/0.3, 
                       est.wh.slope=FALSE, use.tree=FALSE)

# 20,000 takes about 5 mins
totSims = 20000 # 2000 for quick re-runs
thinn = 10

myMCMCstate <- burnin.phybreak(myMCMCstate, ncycles=totSims/2)
myMCMCstate <- sample.phybreak(myMCMCstate, nsample=totSims/(2*thinn), thin=thinn) # is num sims nsample*thin? I think so, based on times


# consensus transition tree
transtree(myMCMCstate, method = "edmonds") # don't know if needed - in documentation
transtree(myMCMCstate, method = "mpc", infection.times = "infector.sd")
plot(myMCMCstate, plot.which = "mpc")

# use get functions
#get.mcmc has all the infector information
mcmc = get.mcmc(myMCMCstate)
#mcmc[i,] is data from ith sample
#eg mcmc[3,"mS"] eg mcmc[3,17:29]

step = totSims/(200*thinn) #10
rowindex = seq(from=step, by=step, to=nrow(mcmc))
pbMatList <- vector("list", length(rowindex)) #nrow(mcmc))
pbTTreeList <- vector("list", length(rowindex)) 


for (i in rowindex){
  pbwiw <- cbind(mcmc[i,(length(mcmc[i,])-length(time_years)):(length(mcmc[i,])-1)],1:length(time_years))
  #print(pbwiw)
  depths = findMRCIs(pbwiw)$mrciDepths
  pbMatList[[i/step]] <- depths # warning, this better divide evenly! no error catching here
  
  rawTTimes <- cbind(mcmc[i,4:(3+length(time_years))])
  tTree <- matrix(nrow=(length(time_years)), ncol=3) # transPhylo ttree format, ?an extra row for the index case, infected from unsampled
  for (row in seq(from=1, to=(length(time_years)))){
    tTree[row,1] <- (first_date + rawTTimes[row])
    tTree[row,2] <- (first_date + time_years[row])
    tTree[row,3] <- pbwiw[row,1]
  }
  tempObj = NULL
  tempObj$nam <- very_short_names #added very_ 24/06
  tempObj$ttree <- tTree
  pbTTreeList[[i/step]] <- tempObj
}
pbdist <- wiwTreeDist(pbMatList, sampled=1:length(time_years))
pbwiwMDS <- dudi.pco(pbdist, scannf=FALSE, nf=3)
saveRDS(pbMatList,file="tree/PhyBreak20000.RData")
saveRDS(pbTTreeList, file="tree/PhyBreak20000_TTreeList.RData")
saveRDS(myMCMCstate, file='tree/MCMCstate20000.RData')

treeNames = seq(from=1, by=1, to=100)
plotGrovesD3(pbwiwMDS, treeNames=treeNames)

pbmed <- wiwMedTree(pbMatList,sampled=1:length(time_years)) # works, pbmed$median is just an index num
pbmedtree <- pbTTreeList[[pbmed$median]]
#pbmedInfo <- extractTTree(record[[pbmed$median]]$ctree)$ttree 
for (i in seq(from=1, to=numSamp)){
    {pbmedtree$nam[i] <- substr(pbmedtree$nam[i],5,6)}#hacky
}
plotTTree(pbTTreeList[[pbmed$median]],w.shape=1.3, w.scale=1/0.3)
networkTPlot(pbTTreeList[[pbmed$median]])
saveRDS(pbTTreeList[[pbmed$median]],file="tree/MedianPBv_Trunc.RData")

# get the phylo tress
# plot(myMCMCstate, plot.which = "sample", samplenr = 34) # this is plain phylogeny

