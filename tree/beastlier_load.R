# process beastlier output

# log file does have infection dates in it, write python code to get

time_years = c(0.015, 0.001, 0.273, 0.733, 2.018, 2.36, 2.64, 3.41, 4.878, 4.968, 5.037, 5.834, 6.579)
first_date = 2005.956
names(time_years) = c("Case27_2005-12-20", "Case28_2005-12-15", "Case29_2006-03-24", "Case30_2006-09-08", "Case37_2007-12-21", 
                      "Case38_2008-04-24", "Case40_2008-08-04", "Case47_2009-05-12", "Case52_2010-10-30", "Case54_2010-12-02", 
                      "Case55_2010-12-27", "Case68_2011-10-14", "Case81_2012-07-12")
short_names = c("Case27", "Case28", "Case29", "Case30", "Case37", "Case38", "Case40", "Case47", "Case52", "Case54", "Case55", "Case68", "Case81")
very_short_names= c("27", "28", "29", "30", "37", "38", "40", "47", "52", "54", "55", "68", "81")


filename = 'tree/Beastlier_VegardJun16c.csv' 
logfilename = 'tree/vegardjun16c.log.txt'
treefilename = 'tree/vegardjun16c.trees.txt'

beastlier_raw_list = scan(file=filename, what=character(), nmax=1000000, sep=',')
numCases = length(time_years)

binfectorList <- vector("numeric", numCases)
binfectedList <- vector("numeric", numCases)

bInfectionTimes <-  vector("numeric", numCases)

numRows = length(beastlier_raw_list)/3
numSims = numRows/numCases
indexi = seq(from=1,to=numSims)
indexj = seq(from=1,to=numCases)

beastlier_raw_log = scan(file=logfilename, what=character(), nmax=1000000, sep='\t') # what default is numeric

raw_offset = 2 # two single-cell header lines
raw_ncols = 128 # the standard CAN CHANGE
col_offset = 17 # 18th col is where infection times start
burnin_offset = 101 # a bit hacky

bPosterior <-  vector("numeric", burnin_offset+numSims)

for (bPost in seq(from=1, to=burnin_offset+numSims)){
  bPosterior[bPost] <- beastlier_raw_log[raw_offset + (bPost*raw_ncols) + 2]
}

# plot hairy caterpillar, exp(-posterior) check vs TransPhylo?
plot(bPosterior, ylab='Posterior probability', xlab='MCMC iterations', type='l') 
# plot post burn-in
plot(bPosterior[101:201], ylab='Posterior probability', xlab='MCMC iterations', type='l')

bMatList <- vector("list", numSims)
bTTreeList <- vector("list", numSims)
bPTreeList <- vector("list", numSims)

for (i in indexi){
  for (j in indexj){
    element = ((i-1)*3*numCases) + (j*3 -2)
    infection_time = raw_offset + (i+burnin_offset)*raw_ncols + col_offset + j
    bInfectionTimes[j] = beastlier_raw_log[infection_time]
    
    infector = beastlier_raw_list[element+1]
    infected = beastlier_raw_list[element+2]
    if (nchar(infector)==6) {infector = as.numeric(substr(infector,5,6))}
    else {infector = as.integer(infector)}
    if (nchar(infected)==6) {infected = as.numeric(substr(infected,5,6))}
    else {infected = as.integer(infected)}
    binfectorList[j] = infector
    binfectedList[j] = infected
  }
  #print(infectorList)
  bwiw <- cbind(Infector=binfectorList, Infectee=binfectedList)
  depths = findMRCIs(bwiw)$mrciDepths
  #print(scwiw)
  bMatList[[i]] <- depths
  
  tTree <- matrix(nrow=(length(time_years)), ncol=3) # transPhylo ttree format, ?an extra row for the index case, infected from unsampled
  for (row in seq(from=1, to=numCases)){
    #tTree[row,1] <- (first_date + time_years[row] - 0.1) # -1 to give us something to draw
    tTree[row,1] <- (first_date + (as.numeric(bInfectionTimes[row])/365)) 
    tTree[row,2] <- (first_date + time_years[row])
    tTree[row,3] <- bwiw[row,1]
  }
  tempObj = NULL
  tempObj$nam <- short_names #use very_short_names for small plots
  tempObj$ttree <- tTree
  bTTreeList[[i]] <- tempObj
}

halltrees = read.nexus(treefilename)
for (tnum in seq(from=1, to=numSims)){
  bPTreeList[[tnum]] <- halltrees[[length(halltrees)-numSims-1+tnum]]
}

bdist <- wiwTreeDist(bMatList, sampled=1:numCases)
bwiwMDS <- dudi.pco(bdist, scannf=FALSE, nf=3)

treeNames = seq(from=1, by=1, to=numSims)
#plotGrovesD3(bwiwMDS, treeNames=treeNames)

#median tree
bMed <- wiwMedTree(bMatList,sampled=1:numCases) # works, pbmed$median is just an index num
#plotTTree(bTTreeList[[bMed$median[1]]],w.shape=1.3, w.scale=1/0.3)
#networkTPlot(bTTreeList[[bMed$median[1]]])
#saveRDS(bTTreeList[[bMed$median[1]]],file="tree/MedianBv_Trunc.RData")

