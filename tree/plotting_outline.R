
# simple visualisation of transmission tree with visNetwork
library(visNetwork)
#' @param ctree Combined tree object 'ctree' from the record produced by inferTTree, for example record[[1]]$ctree where record is the output of inferTTree
#' @author Caroline Colijn
#' @example 
#' vninfo=networkPlot(record[[10]]$ctree)
#' visNetwork(vninfo$nodes,vninfo$edges) %>% visLegend(width=0.2,addNodes=vninfo$lnodes,useGroups=F)
networkPlotColours <- function(ctree,showTimes=T,shapefactor=3) {
  tt=extractTTree(ctree)
  info1=tt$ttree  
  numCases=nrow(info1)
  numSamp=length(tt$nam)  # number sampled; first group
  numUnsamp=nrow(info1)-numSamp; # number unsampled 
  SimpleWiw <- cbind(info1[,3],1:numCases) # infector, infectee
  if (!showTimes)  nodes <- data.frame(id = 1:nrow(info1), label=c(tt$nam, (numSamp+1):numCases) )
  if (showTimes) {
    infTimes=tt$ttree[,1]
    labs=paste(c(tt$nam, (numSamp+1):numCases)," (",0.1*round(10*infTimes),") ",sep="") # should really convert to dd-mm-yy
    nodes <- data.frame(id = 1:numCases,label=labs)
  }
  nodes$groups=c(rep("sampled",numSamp), rep("unsampled",numUnsamp))
  nodes$value=1; nodes$value[nodes$groups=="sampled"]=shapefactor
  colors=brewer.pal(4,"Spectral")
  pal<-colorRampPalette(colors) 
  # early cases have higher values in this ordering: 
  nodes$color=pal(numCases)[numCases-rank(infTimes)+1]
  nodes$shape="circle"
  
  
  # 3 demonstrative colours for the colour key: 
  lnodes <- data.frame(label = c(round(min(infTimes)), round(median(infTimes)),round(max(infTimes))),
                       shape = c( "cirle"), color = pal(3)[c(3,2,1)],size=1)
  edges <- data.frame(from = SimpleWiw[,1], to = SimpleWiw[,2],
                      arrows="to")
  thesource=edges$to[which(edges$from==0)]
  #   nodes$shape[thesource]="star" # thought this made the source look too different
  nodes$color[thesource]="darkgrey"
  # more demonst cols
  lnodes<- data.frame(label=round(quantile(infTimes, seq(0,by=0.1,1))),shape=c("circle"), color=pal(11)[seq(11,1)],size=1);
  #   visNetwork(nodes, edges) %>% visLegend(width=0.3,addNodes=lnodes,useGroups = F)
  return(list(nodes=nodes,edges=edges,lnodes=lnodes))
}

# another function :

## getting weights 
getWeightsFromWIW <-  function(mywiw, edges,numCases) {
  # compute weight for each edge in the edgelist, only using wiw for those transmissions from sampled to sampled cases
  IND=which(edges$from <= numCases & edges$to <= numCases & edges$from>0)
  ws=vapply(1:length(IND), function(x) mywiw[edges$from[IND[x]], edges$to[IND[x]]], FUN.VALUE = 1)
  return(list(IND=IND, widths=ws))
}



#########################################3
# now use this to make the colourful plot 
##########################################

# set up your median entry for your record: 
#medInd= YourMedianEntry

# make info for the visnetwork plot
vninfo=networkPlot(record[[medInd]]$ctree,shapefactor = 1)



# you can specify widths (posterior probs)

# copmuteMatWIW gets a matrix, donors as rows, receivers as columns. i think computeMatWIW is in transphylo

mywiw=computeMatWIW(record)

# need to make a weights vector taking weights from mywiw where both cases are sampled (Roetzer in 1:86, though we could go through row and col names too)
w=getWeightsFromWIW(mywiw,vninfo$edges,86) # CHANGE THE 86 as needed 
vninfo$edges$width=1
vninfo$edges$width[w$IND]=1+round(20*w$widths) # empirical - what looked good

# unsampled to grey

vninfo$nodes$color[vninfo$nodes$groups=="unsampled"]="grey"  # I believe the other cols are set in networkPlot
vninfo$nodes$label=paste("  ", c(extractTTree(record[[medInd]]$ctree)$nam,87:92), "  ",sep="")

# plot and save . lnodes is already created by visNetwork
visNetwork(vninfo$nodes,vninfo$edges) %>% visLegend(width=0.2,addNodes=vninfo$lnodes,useGroups=F,stepY=50) %>% visInteraction(hover=T) 
visSave(file="myvisplot2.html")






