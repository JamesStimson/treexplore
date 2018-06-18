#' @param record MCMC output produced by inferTTree
#' @param k Case whose posterior infection times are to be extracted. Either an integer or a string matching one of the case names in the data
#' @return A vector of posterior infection times for case k 
#' @examples
#' gettimes(record[501:length(record)],"Case1") # posterior infection times for case "Case1", disregarding the first 500 record entries
#' gettimes(record,2) # all posterior infection times for case 2
getInfectionTimes <- function(record,k) {
  if (is.numeric(k)) {
    mytimes= vapply(1:length(record), function(x) { tt=extractTTree(record[[x]]$ctree); return(tt$ttree[k,1])},FUN.VALUE=1);
    return(mytimes)}
  else {
    mytimes= vapply(1:length(record), function(x) { tt=extractTTree(record[[x]]$ctree); ii=which(tt$nam==k); return(tt$ttree[ii,1])},FUN.VALUE=1);
    return(mytimes)}
}

#' @param ctree Combined tree object 'ctree' from the record produced by inferTTree, for example record[[1]]$ctree where record is the output of inferTTree
#' @return Vector of times between becoming infected and infecting others (generation times) in the posterior
#' @examples 
#' getGenerationTimes(record[[1]]$ctree) JS checked OK
getGenerationTimes <-  function(ctree) { tt=extractTTree(ctree)$ttree;
# 3rd column of ttree lists the infectors; exclude source
infectors=tt[,3]; infectors=infectors[infectors!=0];
# times at which each infector infected others:
infothers=tt[tt[,3]!=0,1]; 
# times at which each infector was herself infected:
gotinfd=tt[infectors,1];
return(infothers-gotinfd);}

# number unsampled cases 
#' @param ctree Combined tree object 'ctree' from the record produced by inferTTree, for example record[[1]]$ctree where record is the output of inferTTree
#' @return The number of unsampled cases in ctree
#' @examples 
#' getNumberUnsampled(ctree)
#' JS add $ttree checked OK
getNumberUnsampled <- function(ctree) {tt=extractTTree(ctree)$ttree; return(sum(is.na(tt[,2])))}

# simple visualisation of transmission tree with visNetwork JS checked OK
library(visNetwork)
#' @param ctree Combined tree object 'ctree' from the record produced by inferTTree, for example record[[1]]$ctree where record is the output of inferTTree
networkPlot <- function(ctree) {
    tt=extractTTree(ctree)
    info1=tt$ttree  
    SimpleWiw <- cbind(info1[,3],1:length(info1[,1])) # infector, infectee
  nodes <- data.frame(id = 1:nrow(info1), label=c(tt$nam, rep("Unsampled",nrow(info1)-length(tt$nam)))) # was (length(tt$nam)+1) : nrow(info1))
  edges <- data.frame(from = SimpleWiw[,1], to = SimpleWiw[,2],
                      arrows="to")
  visNetwork(nodes, edges)
  }

#JS let's have a ttree version of this
networkTPlot <- function(tt) {
  info1=tt$ttree  
  SimpleWiw <- cbind(info1[,3],1:length(info1[,1])) # infector, infectee
  #color = c("darkred", "grey", "orange", "darkblue", "purple") USE THIS!!
  nodes <- data.frame(id = 1:nrow(info1), label=c(tt$nam, rep("Unsampled",nrow(info1)-length(tt$nam))), 
                      color=c(rep("lightblue",length(tt$nam)),rep("orange",nrow(info1)-length(tt$nam))),
                      group=c(rep("Sampled",length(tt$nam)),rep("Unsampled",nrow(info1)-length(tt$nam)))) # was (length(tt$nam)+1) : nrow(info1))
  edges <- data.frame(from = SimpleWiw[,1], to = SimpleWiw[,2],
                      arrows="to")
  if(nrow(info1)!=length(tt$nam)){
  visNetwork(nodes, edges) %>%
    visGroups(groupname = "Sampled", color = "lightblue") %>%
    visGroups(groupname = "Unsampled", color = "orange") %>% 
    visLegend(width = 0.1, position = "left", main = "Legend")%>% 
    visNodes(font = list(size = 24))
  }
  else{# doesnt work well with just one group!
    visNetwork(nodes, edges)%>% 
      visNodes(font = list(size = 24))#%>%
      #visEvents(doubleClick = "itplot_dblclick") # THIS MAKES THE IMAGE DISAPPEAR!
  }
}

# requires treespace. 
# create list of wiw information in order to compute transmission tree distances
#' @param record  MCMC output produced by inferTTree
#' @param skipnum Number of record entries to skip, ie skipnum=10 uses the 1st, 11th, 21st ... elements of the record
#' @return list of MRCI information required by wiwTreeDist, one entry for each transmission tree that is included
#' #JS checked OK, use skipnum =100, say
 getTTreeDistInfo <- function(record,skipnum=1) {
   ind=seq(from=1, by=skipnum,to=length(record))
   record=record[ind]
  matList <- lapply(1:length(record), function(x) {
  info <- extractTTree(record[[x]]$ctree)$ttree
  wiw <- cbind(info[,3],1:length(info[,1]))
  findMRCIs(wiw)$mrciDepths
})
return(matList)
}

##################
## NOTES BELOW: using treespace to find distances for these objects
# this function requires knowledge of how many cases are sampled; this must be the same across all entries of the record (or the matList)

# in the ttree format, cases with NA in the second column are unsampled, and they are all listed after all the sampled cases. so the number sampled is: 
####numSamp=max(which(!is.na(extractTTree(record[[1]]$ctree)$ttree[,2])))
#dist <- wiwTreeDist(matList, sampled=1:numSamp) JS comment out

##JS code
####distinfo = getTTreeDistInfo(record,100)
####dist <- wiwTreeDist(distinfo, sampled=1:numSamp)

# use treescape functions to plot the MDS:
#wiwMDS <- dudi.pco(dist, scannf=FALSE, nf=3)
#plotGrovesD3(wiwMDS, tooltip_text=samp)#JS errors

# find a median:
#med <- wiwMedian(matList,sampled=1:86) JS removed
####med <- wiwMedTree(distinfo,sampled=1:numSamp)
#medInfo <- extractTTree(record[[samp[[med$median]]]]$ctree)$ttree
####medInfo <- extractTTree(record[[med$median]]$ctree)$ttree  # JS makes sense?


# matList <- getTTreeDistInfo(record)
  # see computeMatWIW from Transphylo
  #ttt[ttt < 2] <- 2t
  #ttt[1:4,1:4]
  #Like this: M[,c(1,7,9,11,13,15)] but use the NAMES
  #matMean <- lapply(1:numElements, function(x) {
  #  something <- matList[[x]]
  #}
 
 
 
 
 
 
