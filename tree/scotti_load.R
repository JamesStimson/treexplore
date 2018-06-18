# process scotti output

time_years = c(0.015, 0.001, 0.273, 0.733, 2.018, 2.36, 2.64, 3.41, 4.878, 4.968, 5.037, 5.834, 6.579)
first_date = 2005.956
names(time_years) = c("Case27_2005-12-20", "Case28_2005-12-15", "Case29_2006-03-24", "Case30_2006-09-08", "Case37_2007-12-21", 
                      "Case38_2008-04-24", "Case40_2008-08-04", "Case47_2009-05-12", "Case52_2010-10-30", "Case54_2010-12-02", 
                      "Case55_2010-12-27", "Case68_2011-10-14", "Case81_2012-07-12")
short_names = c("Case27", "Case28", "Case29", "Case30", "Case37", "Case38", "Case40", "Case47", "Case52", "Case54", "Case55", "Case68", "Case81")
very_short_names= c("27", "28", "29", "30", "37", "38", "40", "47", "52", "54", "55", "68", "81")


offset_years = c(9.346, 4.07, 13.82, 17.951, 24.777, 3.866, 9.351, 1.878, 0.744, 5.287, 8.192, 2.796, 16.802)

filename = 'tree/vegard15double_wiw.txt' #2451,2576,2601,2851,2876,2976,3076 etc
treefilename = 'tree/vegard15double.trees'

scot_raw_list = scan(file=filename, what=character(), nmax=1000000, sep=',')

numCases = 13
infectorList <- vector("numeric", 13)
infectedList <- vector("numeric", 13)

numRows = length(scot_raw_list)/3
numSims = numRows/numCases
indexi = seq(from=1,to=numSims)
indexj = seq(from=1,to=numCases)

scMatList <- vector("list", 1)
scTTreeList <- vector("list", 75) # WARNING was 88
scPTreeList <- vector("list", 75) # WARNING 
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
  if ((class(depths)!="try-error") & (count<75)){  #75 to get in line with roetzer
    count = count+1
    scMatList[[count]] <- depths
    
    tTree <- matrix(nrow=(numCases), ncol=3) # transPhylo ttree format, ?an extra row for the index case, infected from unsampled
    for (row in seq(from=1, to=(numCases))){
      tTree[row,1] <- (first_date + time_years[row] - offset_years[row]) # -1 to give us something to draw
      tTree[row,2] <- (first_date + time_years[row])
      tTree[row,3] <- scwiw[row,1]
    }
    # go back through the tree and fix any anomalous dates (since we aren't given them, just a range from the inputs)
    for (pass in seq(from=1, to=12)){
    for (row in seq(from=1, to=(numCases))){
      if (tTree[row,3] > 0){    #not the index case
        if (tTree[row,1] < tTree[tTree[row,3],1]){
          tTree[row,1] = tTree[tTree[row,3],1] + (0.01*row) #must have unique times
        }
      }
    }
    } #pass
    tempObj = NULL
    tempObj$nam <- short_names
    tempObj$ttree <- tTree
    scTTreeList[[count]] <- tempObj
  }
}

sctrees = read.nexus(treefilename) # every 20 for 100,000, so 5,000
for (tnum in seq(from=1, to=count) ){
  temp <- sctrees[[2500+tnum*25]]
#   for (i in seq(from=1, to=numCases)){
#     temp$tip.label[i] <- substr(temp$tip.label[i],1,6)
#   }
  scPTreeList[[tnum]] <- temp
}

scdist <- wiwTreeDist(scMatList, sampled=1:numCases)
scwiwMDS <- dudi.pco(scdist, scannf=FALSE, nf=3)

treeNames = seq(from=1, by=1, to=count)
#plotGrovesD3(scwiwMDS, treeNames=treeNames)
nScotti = count

scMed <- wiwMedTree(scMatList,sampled=1:numCases) # works, pbmed$median is just an index num
#plotTTree(scTTreeList[[scMed$median]],w.shape=1.3, w.scale=1/0.3)
#networkTPlot(scTTreeList[[scMed$median]])
#saveRDS(scTTreeList[[scMed$median]],file="tree/MedianSCv_Trunc.RData")
