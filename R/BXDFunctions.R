#'
#' remove rows with NAs from matrix
#'
#'@param plev a matrix
#'
#'@export
#'
rmNarows <- function(plev )
{
  x <- apply(plev,1,function(x){sum(is.na(x))})
  plev <- plev[-which(x>0),]
}

#reads data rm na and log2
prepareData <- function(plev){
  pln <- log2(plev)
  return(pln)
}

##
selectSignificant2 <- function(dd, siglevel = 0.7, minPcor=0)
{
  tmpq5 <- apply(dd,c(2,3),quantile,1-siglevel)
  tmpq95 <- apply(dd,c(2,3),quantile,siglevel)
  diag(tmpq5) <- 0
  diag(tmpq95) <- 0

  indic <- (tmpq5 > minPcor | tmpq95 < -minPcor)

  tmpmedian <- apply(dd,c(2,3),median)
  tmpmedian[!indic] <- 0
  diag(tmpmedian)=0
  return(tmpmedian)
}

selectSignificant <- function(dd, siglevel = 0.8, minPcor=0)
{
  tmpq5 <- apply(dd,c(2,3),quantile,1-siglevel)
  tmpq95 <- apply(dd,c(2,3),quantile,siglevel)
  diag(tmpq5) <- 0
  diag(tmpq95) <- 0

  indic <- (tmpq5 > 0 | tmpq95 < 0)

  tmpmedian <- apply(dd,c(2,3),median)
  tmpmedian[!indic] <- 0
  diag(tmpmedian)=0
  indic <- (tmpmedian > minPcor | tmpmedian < -minPcor)
  tmpmedian[!indic] <- 0
  return(tmpmedian)

}

#'
#' create an igraph
#'
#'@param adjMatrix adjacency matrix
#'@param protnam protein names
#'
#'@export
#'
createGraph <- function(adjMatrix,protnam){
  library(igraph)
  graph.speed=graph.adjacency(adjmatrix=adjMatrix, mode="undirected", weighted=TRUE)
  V(graph.speed)$name = protnam
  temp.degree=degree(graph.speed)
  pname <- protnam
  pname[temp.degree==0]=NA
  V(graph.speed)$label=pname
  return(graph.speed)
}

#'
#' plot a graph - facade to igraphs plot function
#'
#'@param graph an igraph
#'@param myl my layout
#'
#'@export
#'
plotGraphSimple <- function(graph,myl=NULL){
  vdeg <- degree(graph)
  library(RColorBrewer)
  mypalette<-brewer.pal(9,"BuGn")
  vcols <- mypalette[round(rank(vdeg-max(vdeg))/length(vdeg) * 7)+1]

  mypalette <- brewer.pal(9,"YlGnBu")
  weights <- E(graph)$weight+1
  ecols <- round(((weights/max(weights)*8)+1))
  range(ecols)
  ecols <- mypalette[ecols]

  x <- plot(graph, layout=myl,  vertex.size=3, vertex.frame.color="white",  edge.color= ecols,
            vertex.color=vcols,
            vertex.label.font=1,
            vertex.label.cex=0.7,
            vertex.label.color=1)
  return(x)
}

#'
#' plot a graph - facade to igraphs plot function
#'
#'@param graph an igraph
#'@param myl my layout
#'@param cex size of vertex label
#'@param main main title
#'@param vertcol  vertex color
#'
#'@export
#'
plotGraph2 <- function(graph,myl,cex=0.5,main="",vertcol=NULL){
  vdeg <- degree(graph)
  library(RColorBrewer)
  mypalette<-brewer.pal(9,"BuGn")
  #vcols <- mypalette[round(rank(vdeg-max(vdeg))/length(vdeg) * 7)+1]

  mypalette <- brewer.pal(9,"RdBu")
  weights <- E(graph)$color
  weights[weights>0] <- sqrt(weights[weights>0])
  weights[weights<0] <- -1* sqrt(-1*weights[weights<0])
  weights <- weights+1
  ecols <- round(((weights/2 * 8)+1))
  ecols <- mypalette[ecols]

  x <- plot(graph, layout=myl,  vertex.size=7, vertex.frame.color="white",
            edge.color= ecols,
            edge.width = (weights-1)^6*4,
            vertex.label.font=1,
            vertex.label.cex=cex,
            vertex.label.color=1,
            vertex.color = V(graph)$color , main=main)
  return(x)
}



computeBoost <- function(pln,ffunc,nr=250)
{
  result1=ffunc(pln)
  bla <- dim(result1)
  dd <- array(0,dim=c(nr,bla[1],bla[2]))
  for(i in 1:nr){
    print(i)
    ps <- pln[sample(dim(pln)[1],replace=TRUE),]
    dd[i,,] <- ffunc(as.matrix(ps))
  }
  return(dd)
}

computeBoost2 <- function(pln,g1,g2,ffunc,nr=250)
{
  result1 <- ffunc(pln)
  bla <- dim(result1)
  dd <- array(0,dim=c(nr,bla[1],bla[2]))
  for(i in 1:nr){
    print(i)
    ps <- pln[c(
      g1[sample(length(g1),replace=TRUE)]
      ,g2[sample(length(g2),replace=TRUE)]
    ),]
    dd[i,,] <- ffunc(as.matrix(ps))
  }
  return(dd)
}

#'
#' plot a graph - facade to igraphs plot function
#'
#'@param graph an igraph
#'@param myl my layout
#'@param vfc vertex.frame.color
#'@param family vertex.label.family
#'@export
#'
plotGraphAnn <- function(graph,myl,vfc="white",family=1){
  library(RColorBrewer)

  mypalette <- brewer.pal(9,"RdBu")
  weights <- E(graph)$color
  weights[weights>0] <- (weights[weights>0])^0.2
  weights[weights<0] <- -1* (-1*weights[weights<0])^0.2
  weights <- weights+1
  ecols <- round(((weights/2 * 8)+1))
  ecols <- mypalette[ecols]

  x <- plot(graph, layout=myl,  vertex.size=6, vertex.frame.color=vfc, edge.width=4, edge.color= ecols,vertex.label=V(graph)$protname ,vertex.label.family=family,vertex.label.font=1 ,vertex.label.dist=0.4, vertex.label.cex=1 , vertex.label.color=1)
  return(x)
}

#'
#' plot a graph - facade to igraphs plot function
#'
#'@param graph an igraph
#'@param myl my layout
#'@param vfc vertex frame color
#'@param family vertex.label.family
#'@param edge.width edge.width
#'@param vertex.label.cex vertex.label.cex
#'@param vertex.label.color vertex.label.color
#'@export
#'
plotGraphAnn2 <- function(graph,myl,vfc="white",family=1,edge.width=NULL,vertex.label.cex=1,vertex.label.color=1){
  weights <- E(graph)$weight
  cols = rep("pink",length(weights))
  cols[weights > 0] = "cornflowerblue"
  x <- plot(graph, layout=myl,  vertex.size=6, vertex.frame.color=vfc, edge.width=edge.width,
            edge.color= cols,vertex.label=V(graph)$protname ,vertex.label.family=family,
            vertex.label.font=1 ,vertex.label.dist=0.4, vertex.label.cex=vertex.label.cex ,
            vertex.label.color=vertex.label.color)
  return(x)
}


#'
#' quantile normalization
#'
#'@param pln matrix which columns will be normalized
#'
#'@export
#'
myqnormalize <- function(pln){
  tmp <- normalize.quantiles(as.matrix(log(pln)))
  rownames(tmp) <- rownames(pln)
  colnames(tmp) <- colnames(pln)
  tmp <- scale(tmp)
  plns1 <- t(scale(t(tmp)))
}

#get Colors for sample
getColorsforSample <- function(pln){
  samplenames <- colnames(pln)
  cols <- rep("",length(samplenames))
  cols[grep("HFD",samplenames)] <- "gray"
  cols[grep("CD",samplenames)] <- "white"
  return(cols)
}

giveCols2names <- function(protnam,tt,ml){
  mcol <- rep("",length(protnam))
  for(i in 1:length(protnam)){
    x <- (tt[which(tt$ids==protnam[i]),2])
    if(length(x)==0){
      mcol[i] <- "white"
    }
    else if(is.na(x)){
      mcol[i] <- "white"
    }
    else{
      #cat(protnam[i]," " ,ml[[x]],"\n")
      if(length(ml[[x]])==0){
        mcol[i] <- "lightgray"
      }
      else{
        mcol[i] <- ml[[x]]
      }
    }
  }
  return(mcol)
}



getMatchingRows = function(uniprotnam , rnams){
  xx <- NULL
  for(i in 1:length(uniprotnam)){
    dd =  which(uniprotnam[i]  == rnams[,2])
    xx <-c(xx,(dd))
  }
  xx <- xx[!is.na(xx)]
  return(xx)
}


#'
#'
#' plots histograms from upper left, lower left, and lower right quarter of matrix
#'
#'@param tmp matrix
#'@param main main title of figure
#'
#'@export
#'
plotHistCor2Cond = function(tmp,main=""){
  dow = 1:(dim(tmp)[1]/2)
  up = (dim(tmp)[1]/2+1):dim(tmp)[1]
  hist(uppertriang(tmp[dow,dow]),angle=-45,col=1,density=10,xlim=c(-1,1),main=main)
  hist(uppertriang(tmp[up,up]),add=TRUE,col=2,border=T,density=10)
  hist(uppertriang(tmp[dow,up]),add=TRUE,col=3,border=T,density=5,angle=30)
  legend("topleft",c("HFD","CD","HFD/CD"),col=c(1,2,3),lty=c(1,1,1))
}

selectDifferenceInteractions = function(tmp, diffthresh = 0.54){
  dow = 1:(dim(tmp)[1]/2)
  up = (dim(tmp)[1]/2+1):dim(tmp)[1]
  HDFg = tmp[up,up]
  CDg = tmp[dow,dow]
  R = cor(c(HDFg),c(CDg),method="spearman")
  plot(c(HDFg),c(CDg),xlab= "cor(HDF)", ylab="cor(CD)")
  legend("topleft",paste("R=",round(R,digits=2)))

  diffHDFCD = (HDFg - tmp[dow,dow])
  idx=diffHDFCD > diffthresh | diffHDFCD < -diffthresh
  points(HDFg[idx],CDg[idx],col=2)
  diffHDFCD[!idx] = 0
  return(diffHDFCD)
}

#'
#'
#' plot igraph facade
#'
#'@param matr an adjacency matrix
#'@param names protein names
#'@param main main title
#'@param vertcol color of vertices
#'@param repulserad repulsion
#'@param myl my layout - pass your own i graph layout
#'@export
#'
plotGraph <- function(matr, names, main = "cor(HDF)-cor(CD)",vertcol = NULL,repulserad=10e3,myl = NULL){
  graph.speed <- createGraph((matr),names)
  ww = E(graph.speed)$weight
  E(graph.speed)$color=E(graph.speed)$weight
  if(!is.null(vertcol)){
    V(graph.speed)$color=vertcol
  }
  graph.speed <- delete.vertices(graph.speed,which(degree(graph.speed)==0))
  if(is.null(myl)){
    myl <- layout.fruchterman.reingold(graph.speed,repulserad=repulserad,weights=abs(ww))
  }
  plotGraph2(graph.speed,myl,cex=0.6,main=main)
  return(list(graph = graph.speed, myl = myl))
}
#'
#'
#'gets upper left quarter of a matrix
#'
#'@param tmp matrix
#'
#'@export
#'
getDown = function(tmp){
  dow = 1:(dim(tmp)[1]/2)
  return(tmp[dow,dow])
}
#'
#'
#'gets lower right quarter of a matrix
#'
#'@param tmp matrix
#'
#'@export
#'
getUp = function(tmp){
  up = (dim(tmp)[1]/2+1):dim(tmp)[1]
  return(tmp[up,up])
}

#'
#'
#'gets lower left quarter of a matrix
#'
#'@param tmp matrix
#'
#'@export
#'
getDownUp = function(tmp){
  down = 1:(dim(tmp)[1]/2)
  up = (dim(tmp)[1]/2+1):dim(tmp)[1]
  return(tmp[down,up])
}
#'
#'
#' computes graph measures degree, closeness, betweeness
#'
#'@param gdiff an igraph object
#'
#'@export
#'
graphMeasures = function(gdiff){
  res = list(
    degree = centralization.degree(gdiff)$centralization,
    centralization = centralization.closeness(gdiff)$centralization,
    betweeness = centralization.betweenness(gdiff)$centralization
  )
  return(res)
}
