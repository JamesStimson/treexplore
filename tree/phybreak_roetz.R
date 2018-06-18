# phybreak for roetzer data
library(phybreak)
set.seed(2310)
# use phybreak() to create a phybreak object
# input sequences and sampling times
seq_data = read.phyDat("tree/roetzer2013.fasta", format="fasta", type="DNA")

time_years = c(1.834,10.091,2.338,9.505,10.757,0.253,3.922,0.001,0.338,0.582,0.168,0.667,1.001,1.168,0.086,1.168,1.582,1.582,1.667,1.919,1.834,2.086,1.834,2.253,3.256,2.752,2.001,1.919,2.752,3.171,2.585,2.669,3.752,3.752,4.502,4.669,4.502,4.922,4.922,4.922,5.752,5.752,6.42,6.587,7.006,7.672,7.924,10.508,10.259,10.423,10.675,10.842,10.927,11.346,11.423,11.423,11.842,11.261,3.417,9.091,2.752,8.672,8.754,9.006,9.091,9.343,9.587,9.259,10.927,11.59,12.508,12.261,11.59,12.261,8.42,8.505,11.346,10.006,12.59,12.675,12.757,12.842,12.675,13.009,13.094,12.346)
first_date = 1997.75

names(time_years) = names(seq_data)
short_names = names(seq_data)
very_short_names = names(seq_data)
for (row in seq(from=1, to=(length(time_years)))){
  pos = regexpr(':', short_names[row])
  if (pos>1){ very_short_names[row] <- substr(short_names[row],1,pos-1)}
}

# Gamma: The mean and variance are E(X) = a*s and Var(X) = a*s^2, where a is shape (gen 1.3) and s is scale (gen 1/0.3) 

shift = 1 # 1 is deafault

roetzMCMCstate = phybreak(data=seq_data, times=time_years, mu=1e-04, 
                       gen.shape=shift*1.3, gen.mean=shift*1.3/0.3, sample.shape=shift*1.45, sample.mean=shift*1.45/0.3, 
                       wh.model=3, wh.slope=3, est.gen.mean=TRUE, est.sample.mean=TRUE, 
                       prior.mean.gen.mean=shift*1.3/0.3, prior.mean.gen.sd=sqrt(shift*1.3)/0.3, 
                       prior.mean.sample.mean=shift*1.45/0.3, prior.mean.sample.sd=sqrt(shift*1.45)/0.3, 
                       est.wh.slope=FALSE, use.tree=FALSE)

# 2000 takes about 5 mins for roetzer
totSims = 2000#4000 # 2000 for quick re-runs
thinn = 10

roetzMCMCstate <- burnin.phybreak(roetzMCMCstate, ncycles=totSims/2)
roetzMCMCstate <- sample.phybreak(roetzMCMCstate, nsample=totSims/(2*thinn), thin=thinn) # is num sims nsample*thin? I think so, based on times


# consensus transition tree
transtree(roetzMCMCstate, method = "edmonds") # don't know if needed - in documentation
transtree(roetzMCMCstate, method = "mpc", infection.times = "infector.sd")
plot(roetzMCMCstate, plot.which = "mpc")

# use get functions
#get.mcmc has all the infector information
mcmc = get.mcmc(roetzMCMCstate)
#mcmc[i,] is data from ith sample
#eg mcmc[3,"mS"] eg mcmc[3,17:29]

step = totSims/(200*thinn) #10
rowindex = seq(from=step, by=step, to=nrow(mcmc))
pbMatListRoetz <- vector("list", length(rowindex)) #nrow(mcmc))
pbTTreeListRoetz <- vector("list", length(rowindex)) 


for (i in rowindex){
  pbwiw <- cbind(mcmc[i,(length(mcmc[i,])-length(time_years)):(length(mcmc[i,])-1)],1:length(time_years))
  #print(pbwiw)
  depths = findMRCIs(pbwiw)$mrciDepths
  pbMatListRoetz[[i/step]] <- depths # warning, this better divide evenly! no error catching here
  
  rawTTimes <- cbind(mcmc[i,4:(3+length(time_years))])
  tTree <- matrix(nrow=(length(time_years)), ncol=3) # transPhylo ttree format, ?an extra row for the index case, infected from unsampled
  for (row in seq(from=1, to=(length(time_years)))){
    tTree[row,1] <- (first_date + rawTTimes[row])
    tTree[row,2] <- (first_date + time_years[row])
    tTree[row,3] <- pbwiw[row,1]
  }
  tempObj = NULL
  tempObj$nam <- very_short_names
  tempObj$ttree <- tTree
  pbTTreeListRoetz[[i/step]] <- tempObj
}
pbdistRoetz <- wiwTreeDist(pbMatListRoetz, sampled=1:length(time_years))
pbwiwMDSRoetz <- dudi.pco(pbdistRoetz, scannf=FALSE, nf=3)
saveRDS(pbMatListRoetz,file="tree/PhyBreak2000Roetz.RData")
saveRDS(pbTTreeListRoetz, file="tree/PhyBreak2000_TTreeListRoetz.RData")
saveRDS(roetzMCMCstate, file='tree/MCMCstate2000Roetz.RData')

treeNames = seq(from=1, by=1, to=100)
plotGrovesD3(pbwiwMDSRoetz, treeNames=treeNames)

pbmed <- wiwMedTree(pbMatListRoetz,sampled=1:length(time_years)) # works, pbmed$median is just an index num
#pbmedInfo <- extractTTree(record[[pbmed$median]]$ctree)$ttree 
#plotTTree(extractTTree(record[[pbmed$median]]$ctree),w.shape=1.3, w.scale=1/0.3)
#networkPlot(record[[pbmed$median]]$ctree)

# get the phylo tress
# plot(roetzMCMCstate, plot.which = "sample", samplenr = 34) # this is plain phylogeny

