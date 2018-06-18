# process beastlier output roetzer


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

filename = 'tree/Beastlier_Roetz_Jan19.csv' 
logfilename = 'tree/roetzerjan19.log.txt'
treefilename = 'tree/roetzerjan19.trees.txt'

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
raw_ncols = 639 # num cols in lof file
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

# Column 2: posterior, Column 4: likelihood
# for (i in indexi){
#   for (j in indexj){
#     infection_time = raw_offset + (i+burnin_offset)*raw_ncols + col_offset + j
#     if (j==3) {print(beastlier_raw_log[infection_time])}
#   }
# }


bMatListRoetz <- vector("list", numSims)
bTTreeListRoetz <- vector("list", numSims)
bPTreeListRoetz <- vector("list", numSims)

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
  bMatListRoetz[[i]] <- depths
  
  tTree <- matrix(nrow=(length(time_years)), ncol=3) # transPhylo ttree format, ?an extra row for the index case, infected from unsampled
  for (row in seq(from=1, to=numCases)){
    #tTree[row,1] <- (first_date + time_years[row] - 0.1) # -1 to give us something to draw
    tTree[row,1] <- (first_date + (as.numeric(bInfectionTimes[row])/365)) # TESTED 12/06/17
    tTree[row,2] <- (first_date + time_years[row])
    tTree[row,3] <- bwiw[row,1]
  }
  tempObj = NULL
  tempObj$nam <- very_short_names
  tempObj$ttree <- tTree
  bTTreeListRoetz[[i]] <- tempObj
}

halltrees = read.nexus(treefilename)
for (tnum in seq(from=1, to=numSims)){
  bPTreeListRoetz[[tnum]] <- halltrees[[length(halltrees)-numSims-1+tnum]]
}

bdist <- wiwTreeDist(bMatListRoetz, sampled=1:numCases)
bwiwMDS <- dudi.pco(bdist, scannf=FALSE, nf=3)

treeNames = seq(from=1, by=1, to=numSims)
#plotGrovesD3(bwiwMDS, treeNames=treeNames)



