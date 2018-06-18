#from caroline 03/05 adapted from 10/05 by JS
library(TransPhylo)
library(ape)

# 15/05/17 JS add set seed
set.seed(2310)

roetree_nwk=read.tree(file="roetzer2013.nwk") #equivalent I believe
roetree=read.nexus(file="roetzer2013.nex")
rfasta=read.dna("roetzer2013.fasta",format="fasta")

dates=gsub("^(.*?):","",rownames(rfasta))

LastDate=sort(dates)[length(dates)] # 2010-11-01 # JS use length

LastDate=2010+(as.numeric(as.Date(LastDate)-as.Date("2010-01-01")))/365; # now in format for ptreeFromPhylo

# convert starting tree to Xavier's format:
roestart = ptreeFromPhylo(roetree, dateLastSample = LastDate)

#last argument correct?
roerecord=inferTTree(roestart,w.shape=1.3, w.scale=1/0.3, ws.shape=1.45,ws.scale=1/0.3, mcmcIterations=10000, thinning=10, startNeg=1.5, startOff.r=1, startOff.p=0.5, startPi=0.9, updateNeg=F, updatePi=F, dateT=LastDate)
plotCTree(roerecord[[1000]]$ctree)
#plotTTree(roerecord[[1000]]$ttree)
roecontree=consTTree(roerecord[1000:1000])
plotTTree(roecontree,w.shape=1.3, w.scale=1/0.3)
