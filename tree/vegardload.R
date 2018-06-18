library(TransPhylo)
library(ape)

#ITERATIONS *10 (may 5) back may8 may9: try 1million (=*100)

# 15/05/17 JS add set seed
set.seed(2310)

# read input phylogenetic tree THIS IS THE 13 leaf new version
sclade=read.tree(file="sclade-newlabs.nwk")

# convert starting tree to Xavier's format:
startpoint=ptreeFromPhylo(sclade, dateLastSample = 2012.5) # is this accurate enough? 2012.53

# run the inference! WHAT DOES dateT do? Ans: sets end time of outreak. Make its last sample time + 1 year at least, so 2013.6
record=inferTTree(startpoint,w.shape=1.3, w.scale=1/0.3, ws.shape=1.45,ws.scale=1/0.3, mcmcIterations=1000, thinning=100, startNeg=1.5, startOff.r=1, startOff.p=0.5, startPi=0.9, updateNeg=F, updatePi=F, dateT=2013.6)

# plot one of the trees
plotCTree(record[[1000]]$ctree)

# consensus tree from part of the posterior: 
contree=consTTree(record[500:1000])
plotTTree(contree,w.shape=1.3, w.scale=1/0.3) # actually we should be using posterior rather than prior shape scale here but 
# anyway this gives a visualisation... 

# ploting some output requires two simple functions called gettimes and gentimes:  
# Xavier's ttree format: (row # is case ID); time of infectoin; time of sampling or NA; infector. 
# times of infection for each case 
gettimes= function(record,k) {
  if (is.numeric(k)) {
   mytimes= vapply(1:length(record), function(x) { tt=extractTTree(record[[x]]$ctree); return(tt$ttree[k,1])},FUN.VALUE=1);
  return(mytimes)}
  else {
    mytimes= vapply(1:length(record), function(x) { tt=extractTTree(record[[x]]$ctree); ii=which(tt$nam==k); return(tt$ttree[ii,1])},FUN.VALUE=1);
  return(mytimes)}
}

# generation times:
gentimes <-  function(ftree) { tt=extractTTree(ftree)$ttree;
infectors=tt[,3]; infectors=infectors[infectors!=0]; infothers=tt[tt[,3]!=0,1];  gotinfd=tt[infectors,1]; return(infothers-gotinfd);}

# number unsampled cases 
unsampnum <- function(ftree) {tt=extractTTree(ftree)$ttree; return(sum(is.na(tt[,2])))}

# sets up the 2x2 plot. see also: layout "R function of the day"
par(mfrow=c(2,2)) 

# plot 1
plot(sapply(record,function(x) x$pTTree+x$pPTree),ylab='Posterior probability',xlab='MCMC iterations',type='l')
# plot 2
plot(sapply(record,function(x) x$pTTree),ylab='Ttree probability',xlab='MCMC iterations',type='l')

# plot 3
hist(unlist(sapply(1:length(record), function(x) gentimes(record[[x]]$ctree))),breaks=50,freq=FALSE,xlab="generation times", main="Histogram of generation times, with prior");
lines(x,dgamma(x,shape=1.3,scale=1/0.3))

# plot 4
hist(unlist(sapply(1:length(record), function(x) unsampnum(record[[x]]$ctree))),breaks=10,freq=FALSE,xlab="number unsampled", main="Histogram of number unsampled");

# dev.off()


