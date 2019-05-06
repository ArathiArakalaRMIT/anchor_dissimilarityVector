############################################################################
#aim: Statistical Summary of network topological features
#     features investigated are 
#     num of vertices, num of edges, pl2, eyes, isolated nodes, max degree, num of components, size of largest comp, size of 2 largest comp,
#     var of degree distribution, avg degree, norm vcount, norm ecount, entropy, spectral radius, spectral gap.
#    

#   
#written by: Arathi Arakala
#written on: 17 June 2012
#input: Biometric Graphs
#output: Summary stats and histogram.
##############################################################################

#required packages
#=================
#setwd("I:/Academic/Research/Projects and Grants/Biometric Databases/Graph Databases")

#setwd("D:/Arathi/Graph Databases") #copied graphs to local machine to increase speed
library(splancs)
library(spatial)
library(network)
library(clue)
library(sna)
source("basicFunctions.R")
source("graphMatchingFunc_lib.R")

#Functions defined in this file
################################

#getTopology(g) # for the input graph g, we output the topological features of this graph - #nodes, #edges, degree distribution, and so on

####################################

getTopologySpectral<-function(g, g1, g2) {
  
 
  
  A<-as.matrix.network(g, matrix.type="adjacency")
  I<-as.matrix.network(g, matrix.type="incidence")
  E<-as.matrix.network(g, matrix.type="edgelist")
  
  
  
  # Connection matrix
  C<-matrix(0,dim(A)[1],dim(A)[2])
  if(dim(E)[1]>0){   
      for(e in 1:dim(E)[1]){
        row<-E[e,1]
        col<-E[e,2]
        C[row,col]<-C[row,col]+1
        C[col,row]<-C[col,row]+1 #due to symmetry of C
        }
  }
  #Degree matrix
  D<-matrix(0,dim(A)[1],dim(A)[2])
  for(i in 1: dim(I)[1]){
    D[i,i]<-sum(I[i,])
  }
  
  #Laplacian matrix
  L=C-D
  
  #---------------------------------------
  #1
  vertexCount<-network.size(g)
  
  
  #2
  edgeCount<-network.edgecount(g)
  
  #3
  A2<-A%*%A
  A2_nocircuit<-A2
  diag(A2_nocircuit)<-0
  pl2Count<-0.5*sum(A2_nocircuit)
  
  #4
  eyesCount<-0.5*length(which(C==2))
  
  #5
  isolatedCount<-component.dist(g)$cdist[1]
  
 
  
  #6
  degreeList<-numeric()
  for(i in 1:vertexCount){
    degreeList<-c(degreeList,sum(A[i,]))
  }
  maxDegreeCount<-max(degreeList)
  
  #7
 componentsCount<-components(g)
  
  #8
  largestComponentSize<-max(component.dist(g)$csize)
  
  #9
  indicesOrdered<-order(component.dist(g)$csize,decreasing=TRUE)[1:2]  
  twoLargestComponentSize<-sum(component.dist(g)$csize[indicesOrdered])
  
  #10
  varDegreeDisbn<-var(degreeList)
  
  #11                               
  averageDegree<-sum(degreeList)/length(degreeList)
  
  #12
#   A3<-A2%*%A
#   triangleCount<-sum(diag(A3))/6
  vertexCount_norm<-vertexCount/sqrt(network.size(g1)*network.size(g2))
   #vertexCount_norm<-vertexCount/network.size(g2)
  
  
  #13
#   clusteringcoefft<-0
#   if(pl2Count>0){
#     clusteringcoefft<-triangleCount/pl2Count
#   }
#   
  if(sqrt(network.edgecount(g1)*network.edgecount(g2))==0){edgeCount_norm<-0}
  if(sqrt(network.edgecount(g1)*network.edgecount(g2))!=0){edgeCount_norm<-edgeCount/sqrt(network.edgecount(g1)*network.edgecount(g2))}
  
#   if(network.edgecount(g2)==0){edgeCount_norm<-0}
#   if(network.edgecount(g2)!=0){edgeCount_norm<-edgeCount/network.edgecount(g2)}
# #   
  
  
  #14
#   A3_nocircuit<-A3
#   diag(A3_nocircuit)<-0
#   A3_nocircuit<- A3_nocircuit-A #remove paths of length 1
#   pl3Count<-0.5*sum(A3_nocircuit)-(2*pl2Count)
  degreeDisbn<-array(0,c(1,maxDegreeCount+1))
  for(i in 1:(maxDegreeCount+1)){    
   degreeDisbn[i]<-length(which(degreeList==(i-1)))     
  }#end of i loop
  h<-degreeDisbn/vertexCount
  entropy<- -sum(h*log2(h))
    
  #15
  nonTrivIndices<-which(eigen(A)$values!=0,arr.ind=TRUE)
  spectralRadius<-0
  if(length( nonTrivIndices)>0) { spectralRadius<-max(eigen(A)$values[nonTrivIndices]) }
  
  #16
  nonTrivIndices<-which(eigen(L)$values!=0,arr.ind=TRUE)
  spectralGap<-0
  if(length( nonTrivIndices)>0) { spectralGap<-max(eigen(L)$values[nonTrivIndices]) }
  
  
  #distribution<-hist(degreeList, breaks=seq(0,max(degreeList)+1),right=FALSE, include.lowest=FALSE)
  
   #17
  largestComponentSize_norm<-largestComponentSize/sqrt(network.size(g1)*network.size(g2))
  
list(vertexCount=vertexCount,edgeCount=edgeCount, pl2Count=pl2Count, eyesCount=eyesCount, isolatedCount=isolatedCount, maxDegreeCount=maxDegreeCount, componentsCount=componentsCount, largestComponentSize=largestComponentSize, twoLargestComponentSize= twoLargestComponentSize, varDegreeDisbn=varDegreeDisbn, averageDegree=averageDegree, vertexCount_norm=vertexCount_norm, edgeCount_norm=edgeCount_norm, entropy=entropy,  spectralRadius= spectralRadius, spectralGap= spectralGap, largestComponentSize_norm=largestComponentSize_norm)
}



graphTopologyMore<-function(g){
  
  #1 clustering coefficient
  clcoeff<-grtans(g, diag=FALSE, mode="graph", measure = c("strong"), use.adjacency = TRUE)
  
  # path distribution
 kpath<-kpath.census(g, maxlen = 3, mode = "graph", tabulate.by.vertex = TRUE, path.comembership ="none", dyadic.tabulation = "none")                                           
  
  #cycle distribution
  kcycle<-kcycle.census(g, maxlen = 3, mode = "graph", tabulate.by.vertex = TRUE, cycle.comembership = "none")
  
  
  #connectedness
 connected<- connectedness(g)
  
  #return the ouptput measures.
  list()
  
  
}#end of graph TopologyMore


#test code here

# biometric<-"FP"
# sourceid<-"Manual"
# g1<-getGraph(biometric, sourceid, 1,2)
# windows()
# plotspatialgraph(g1)
# 
# g1_top<-getTopologySpectral(g1)
