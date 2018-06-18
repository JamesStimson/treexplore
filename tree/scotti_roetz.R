# process scotti output

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

offset_years = c(9.9529585,1.6865961,4.5327562,3.6839701,2.470218,5.7481875,6.2198605,6.1036488,4.5091653,11.425832,5.8435252,3.9179613,2.8659211,2.1009223,0.2205783,1.3118954,5.7339976,3.8272818,3.6034003,7.7337611,4.2518909,3.2843367,4.0030684,0.2235413,4.0182371,2.5638372,3.5393764,5.0645723,1.0609842,1.0663169,9.7075246,8.512788,0.3327234,2.6308973,
3.9309829,1.4121784,3.9377321,1.0628362,1.4729052,4.5622304,2.6003722,4.4407747,3.9668094,1.5091882,1.1283084,2.5653427,5.9247179,1.0281529,7.3444355,1.4147986,7.9887674,6.8908888,1.7529496,4.278636,1.5702841,7.7981982,3.2499441,0.7033934,7.9573503,1.1201402,3.2846819,5.4268216,4.0112322,5.623314,2.250157,10.088073,11.465869,5.089577,5.1206253,1.3254505,0.9140015,0.3320062,1.3040585,3.6847144,5.3789205,12.458289,6.6379322,4.1540157,5.4753854,0.4874735,1.1293379,6.5095194,9.6672208,4.0645888,2.6435934,8.7168257)

filename = 'tree/roetzer_trans_wiw.txt' 
treefilename = 'tree/roetzer_output.trees' # 10,000 every 2 making 5000 again

scot_raw_list = scan(file=filename, what=character(), nmax=1000000, sep=',')

numCases = length(time_years)
infectorList <- vector("numeric", numCases)
infectedList <- vector("numeric", numCases)

numRows = length(scot_raw_list)/3
numSims = numRows/numCases
indexi = seq(from=1,to=numSims)
indexj = seq(from=1,to=numCases)

scMatListRoetz <- vector("list", 1)
scTTreeListRoetz <- vector("list", 75) # WARNING 
scPTreeListRoetz <- vector("list", 75) # WARNING 
count = 0 
for (i in indexi){
  for (j in indexj){
    element = ((i-1)*3*numCases) + (j*3 -2)
    #print(scot_raw_list[element])
    infector = scot_raw_list[element+1]
    infected = scot_raw_list[element+2]
    if (nchar(infector)==6) {infector = as.numeric(substr(infector,5,6))}
    else {infector = as.integer(infector)}
    if (nchar(infected)==6) {infected = as.numeric(substr(infected,5,6))}
    else {infected = as.integer(infected)}
    infectorList[j] = infector
    infectedList[j] = infected
  }
  #print(infectorList)
  
  scwiw <- cbind(Infector=infectorList, Infectee=infectedList)
  depths = try(findMRCIs(scwiw)$mrciDepths, TRUE)
  if (class(depths)!="try-error"){
    count = count+1
    scMatListRoetz[[count]] <- depths
    
    tTree <- matrix(nrow=(length(time_years)), ncol=3) # transPhylo ttree format, ?an extra row for the index case, infected from unsampled
    for (row in seq(from=1, to=(length(time_years)))){
      tTree[row,1] <- (first_date + time_years[row] - offset_years[row]) # -1 to give us something to draw
      tTree[row,2] <- (first_date + time_years[row])
      tTree[row,3] <- scwiw[row,1]
    }
    # go back through the tree and fix any anomalous dates (since we aren't given them, just a range from the inputs)
    for (pass in seq(from=1, to=10)){
      for (row in seq(from=1, to=(numCases))){
        if (tTree[row,3] > 0){    #not the index case
          if ((tTree[row,1] < tTree[tTree[row,3],1])&((tTree[tTree[row,3],1]+(0.01*row)) < tTree[row,2])) {
            tTree[row,1] = tTree[tTree[row,3],1] + (0.01*row) #must have unique times
          }
        }
      }
     } #pass
    tempObj = NULL
    tempObj$nam <- very_short_names
    tempObj$ttree <- tTree
    scTTreeListRoetz[[count]] <- tempObj # corrected i to count
  }
}

sctrees = read.nexus(treefilename) # every 20 for 100,000, so 5,000
for (tnum in seq(from=1, to=count)){
  scPTreeListRoetz[[tnum]] <- sctrees[[2500+tnum*25]]
}

scdist <- wiwTreeDist(scMatListRoetz, sampled=1:numCases)
scwiwMDS <- dudi.pco(scdist, scannf=FALSE, nf=3)

treeNames = seq(from=1, by=1, to=count)
#plotGrovesD3(scwiwMDS, treeNames=treeNames)
nScottiRoetz = count

scMed <- wiwMedTree(scMatListRoetz,sampled=1:numCases) # works, pbmed$median is just an index num

scMedRoetz <- wiwMedTree(scMatListRoetz,sampled=1:numCases) # works, pbmed$median is just an index num
#plotTTree(scTTreeListRoetz[[scMedRoetz$median]],w.shape=1.3, w.scale=1/0.3)
#networkTPlot(scTTreeListRoetz[[scMedRoetz$median]])
#saveRDS(scTTreeListRoetz[[scMedRoetz$median]],file="tree/MedianSCr.RData")
