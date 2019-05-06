### This file must have basic functions to read in files and construct the graph ###
### We must have functions to plot a graph as well ###

### call the libraries needed ###
library(splancs)
library(spatial)
library(network)
library(clue)
library(sna)
library(ggplot2)
# W_o001_L_S1_Nr4_bpoints
# P_o001_R_S1_Nr1_links
#1. 
##############################################################################
#This function creates a graph object with edge and vertex attributes.
# biometric->P FOR 'Palm Vein, AND W FOR Wrist Vein
# person-> person number
# session -> session number
# number -> the number
# LR -> left or right sample
##############################################################################

getGraph<-function(biometric,person,session,number, LR){
  folderName<-""
  if(biometric=='P')folderName<-"PUT PalmVein"
  if(biometric=='W')folderName<-"PUT WristVein"
  person_prefix<-""
  if(person<10) person_prefix<-"00"
  if(person>=10 && person<100) person_prefix<-"0"
  
  
  inputNodesFile<-""
  inputEdgesFile<-""

    inputNodesFile<-paste(getwd(),"/BiometricDatabases/",folderName,"/", biometric, "_o",person_prefix,person,"_", LR, "_S",session,"_Nr",number,"_bpoints.txt",sep="")
    inputEdgesFile<-paste(getwd(),"/BiometricDatabases/",folderName,"/", biometric, "_o",person_prefix,person,"_", LR, "_S",session,"_Nr",number,"_links.txt",sep="")

  #Create graph
  #-------------
  v <- read.table(inputNodesFile,header=FALSE, sep="\t", fill=TRUE)
  N<-length(v[,1])
  e<-read.table(inputEdgesFile,header=FALSE, sep="\t", fill=TRUE)
  e<-as.matrix(e)

  #arrange edges with lower vertex index first and higher vertex index second #added on 24/02/2015
  toswap<-which(e[,1]>e[,2])
  for(k in toswap){
    val<-e[k,]
    newval<-c(val[2], val[1])
    e[k,]<-newval
  }
  #remove repeated edges
  e<-unique(e)
  
  # create a network using the edgelist
  g <- as.network.matrix(e, matrix.type="edgelist",directed=FALSE)
  gsize <- network.size(g)
  
  #if there are isolated nodes, add them
  if (N-gsize>0) add.vertices(g,N-gsize)  #add isolated vertices if necessary
  set.vertex.attribute(g, "xcoord", v[,2])
  set.vertex.attribute(g, "ycoord", v[,3])
  set.vertex.attribute(g,"type", v[,4]) 
  
  #delete edges where start and end nodes are same.
  todelete<-which(e[,1]==e[,2])
  delete.edges(g,todelete)
  
  
  lengths<-rep(0,times=length(e[,1]) )
  slopes<-rep(0,times=length(e[,1]) )
  
  for (i in 1:length(e[,1])){
    
    start <- e[i,1]
    end <- e[i,2]
    
    lengths[i] <- sqrt((v[start,2]-v[end,2])**2+(v[start,3]-v[end,3])**2)
    
    if((v[end,2]-v[start,2])==0){slopes[i]=90}
    if((v[end,2]-v[start,2])!=0){
      if (v[start,2]<v[end,2]) slopes[i] <- (180/pi)*atan((v[end,3]-v[start,3])/(v[end,2]-v[start,2]))
      if (v[start,2]>=v[end,2]) slopes[i] <- (180/pi)*atan((v[start,3]-v[end,3])/(v[start,2]-v[end,2]))
      
    }
    
    lengths[i]<-round(lengths[i],2)
    slopes[i]<-round(slopes[i],2)
    #print(i)
  }
  
  set.edge.attribute(g,"edgelength",lengths)
  set.edge.attribute(g,"edgeslope",slopes)
  
 
    #return graph
  #-------------
  g
  
}#end of getGraph function
###########################################################################

#### 2.
#################################################################################
############# START plotspatialgraph - plots g as a spatial graph ###############
plotspatialgraph <- function(g){
  
  true_edges <- as.matrix(g, matrix.type="edgelist")
  xvals <- get.vertex.attribute(g, "xcoord")
  yvals <- get.vertex.attribute(g, "ycoord")
  true_vertices <- cbind(xvals, -yvals)
  #true_vertices <- cbind(xvals, yvals)
  #true_vertices <- cbind(yvals, xvals)
  
  pointmap(true_vertices, add=FALSE, pch=20, col="darkblue", cex=1 , xaxt="n",yaxt="n")
  #textxy(xvals,-yvals,seq(1:dim(true_vertices)[1]) )
  #textxy(xvals,yvals,seq(1:dim(true_vertices)[1]) )
  
  num_edges <- length(true_edges[,1])
  if(num_edges>0){
    for (i in 1:num_edges){
      start <- true_edges[i,1]
      end <- true_edges[i,2]
      added_edge <- matrix(c(true_vertices[start,1],true_vertices[end,1],true_vertices[start,2],true_vertices[end,2]), 2, 2)
      lines(added_edge, lwd=1.5, col="darkblue")
    }
  }#end of if
  
}
############# END plotspatialgraph - plots g as a spatial graph ###############

### 3.
##################################################################################################
## START plotreggraphs
###################################################################################################
plotreggraphs<-function(g1,g2){
  
  g1_vertices<-cbind(get.vertex.attribute(g1,"xcoord"), get.vertex.attribute(g1,"ycoord"))
  g2_vertices<-cbind(get.vertex.attribute(g2,"xcoord"), get.vertex.attribute(g2,"ycoord"))
  
  g1_edgelist<-as.matrix.network(g1,matrix.type="edgelist")
  g2_edgelist<-as.matrix.network(g2,matrix.type="edgelist")
  
  #First "trick" R into making the plot big enough to see both graphs
  #-----------------------------------------------------------------------
  all_vertices<-rbind(g1_vertices,g2_vertices)
  pointmap(all_vertices, add=FALSE, pch=20, col="white",cex=1 , xaxt="n",yaxt="n")
  
  #plot vertices of g1
  #====================
  
  pointmap(g1_vertices, add=TRUE, pch=20, col="blue")
  #textxy(g1_vertices[,1],g1_vertices[,2],seq(1:dim(g1_vertices)[1]),dcol="blue" )
  
  #plot edges of g1
  #===================
  
  for (i in 1:length(g1_edgelist[,1])){
    
    start <- g1_edgelist[i,1]
    end <- g1_edgelist[i,2]
    
    added_edge <- matrix(c(g1_vertices[start,1],g1_vertices[end,1],g1_vertices[start,2],g1_vertices[end,2]), 2, 2)
    
    lines(added_edge, lwd=2, col="blue", lty=1)
  }
  
  
  #plot vertices of g2
  #====================
  
  pointmap(g2_vertices, add=TRUE, pch=16, col="darkgreen")
  #textxy(g2_vertices[,1],g2_vertices[,2],seq(1:dim(g2_vertices)[1]),dcol="red")
  
  #plot edges of g2
  #===================
  
  for (i in 1:length(g2_edgelist[,1])){
    
    start <- g2_edgelist[i,1]
    end <- g2_edgelist[i,2]
    
    added_edge <- matrix(c(g2_vertices[start,1],g2_vertices[end,1],g2_vertices[start,2],g2_vertices[end,2]), 2, 2)
    
    lines(added_edge, lwd=1.5, col="darkgreen")
  }
  
 abline(h=0, col="red", lty=3, lwd=0.5)
 abline(v=0,col="red", lty=3, lwd=0.5 )
  
}#end of function

###4.
###########################################################################################################
############# START rapidscore - "rapid" **count** of "common" nodes for two graphs g1 & g2 ###############
rapidscore <- function(g1, g2, tol){
  
  #tol = size of bounding box to determine a vertex match
  
  xvals <- get.vertex.attribute(g1, "xcoord")
  yvals <- get.vertex.attribute(g1, "ycoord")
  vertices_enr <- cbind(xvals, yvals)
  num_vertices_enr <- length(vertices_enr[,1])
  
  xvals <- get.vertex.attribute(g2, "xcoord")
  yvals <- get.vertex.attribute(g2, "ycoord")
  vertices_que <- cbind(xvals, yvals)
  num_vertices_que <- length(vertices_que[,1])
  
  match_score <- 0
  vertices_enr_taken <- rep(FALSE, times=num_vertices_enr)
  vertices_que_taken <- rep(FALSE, times=num_vertices_que)
  
  for (i in 1:num_vertices_enr){   
    for (j in 1:num_vertices_que){
      
      if (vertices_enr_taken[i]) next
      if (vertices_que_taken[j]) next
      
      diff <- sqrt((vertices_enr[i,1]-vertices_que[j,1])**2+(vertices_enr[i,2]-vertices_que[j,2])**2)
      
      if (diff<=tol) {
        match_score <- match_score + 1
        vertices_enr_taken[i]<-TRUE
        vertices_que_taken[j]<-TRUE
      }
    }
  }
  match_score <- 1 - match_score/(sqrt(num_vertices_enr*num_vertices_que))
  match_score
}
############# END rapidscore - "rapid" **count** of "common" nodes for two graphs g1 & g2 ###############

###5. 
##############################################################################################################################################
############# START shiftrotate - returns a rotated & shifted version of a spatial graph, with respect to a particular edge where ############
#############                     the start of the edge is the origin and the +ve x-axis is in the direction of the edge.         ############
shiftrotate <- function(g, edge){
  
  edges <- as.matrix.network(g, matrix.type="edgelist")
  
  xvals <- get.vertex.attribute(g, "xcoord")
  yvals <- get.vertex.attribute(g, "ycoord")
  vertices <- cbind(xvals, yvals)
  num_vertices <- length(vertices[,1])
  
  #APPLY shift and rotation to g
  #=============================
  
  start <- edges[edge,1]
  end <- edges[edge,2]
  
  if (vertices[start,1]>vertices[end,1]) {
    #reverse labels of "start" and "end", so that start always refers to the node having the smallest x-coordinate in the original co-ordinate system
    start <- edges[edge,2]
    end <- edges[edge,1]
  }
  
  if (vertices[start,1]==vertices[end,1]) theta <- 90
  if (vertices[start,1]<vertices[end,1]) theta <- (180/pi)*atan((vertices[end,2]-vertices[start,2])/(vertices[end,1]-vertices[start,1]))
  if (vertices[start,1]>vertices[end,1]) theta <- (180/pi)*atan((vertices[start,2]-vertices[end,2])/(vertices[start,1]-vertices[end,1]))
  
  newvertices <- vertices
  newervertices <- vertices
  
  #shift origin to co-ordinates of start
  newvertices[,1] <- newvertices[,1]-vertices[start,1]
  newvertices[,2] <- newvertices[,2]-vertices[start,2]
  
  #rotate all points such that the edge (start, end) becomes a secant coinciding with the positive x-axis.
  theta <- -(pi/180)*theta  #switches to radians
  for (i in 1:num_vertices){
    newervertices[i,1] <- newvertices[i,1]*cos(theta)-newvertices[i,2]*sin(theta)
    newervertices[i,2] <- newvertices[i,2]*cos(theta)+newvertices[i,1]*sin(theta)
  }
  
  
  # g_aligned <- as.network.matrix(edges, matrix.type="edgelist", directed=FALSE)
  # gsize <- network.size(g_aligned)
  # if (num_vertices-gsize>0) add.vertices(g_aligned,num_vertices-gsize)  #add isolated vertices if necessary
  g_aligned<-g
  
  set.vertex.attribute(g_aligned, "xcoord", newervertices[,1])
  set.vertex.attribute(g_aligned, "ycoord", newervertices[,2])
  
  #update length and slope of g_aligned
  #------------------------------------
  
  lengths<-rep(0,times=length(edges[,1]) )
  slopes<-rep(0,times=length(edges[,1]) )
  
  for (i in 1:length(edges[,1])){
    
    start <- edges[i,1]
    end <- edges[i,2]
    
    lengths[i] <- sqrt((newervertices[start,1]-newervertices[end,1])**2+(newervertices[start,2]-newervertices[end,2])**2)
    
    if(newervertices[start,1]==newervertices[end,1]) slopes[i]<- 90
    if (newervertices[start,1]<newervertices[end,1]) slopes[i] <- (180/pi)*atan((newervertices[end,2]-newervertices[start,2])/(newervertices[end,1]-newervertices[start,1]))
    if (newervertices[start,1]>newervertices[end,1]) slopes[i] <- (180/pi)*atan((newervertices[start,2]-newervertices[end,2])/(newervertices[start,1]-newervertices[end,1]))
    
    lengths[i]<-round(lengths[i],2)
    slopes[i]<-round(slopes[i],2)
    
  }
  
  set.edge.attribute(g_aligned,"edgelength",lengths)
  set.edge.attribute(g_aligned,"edgeslope",slopes)
  
  
  list(g_aligned=g_aligned, theta=theta)
}
############# END shiftrotate - returns a rotated & shifted version of a spatial graph, with respect to a particular edge where ############
#############                   the start of the edge is the origin and the +ve x-axis is in the direction of the edge.         ############
############################################################################################################################################


### 6.
############################################################################################################################################
### This function returns a list of graphs from a person
getPersonGraphList<-function(biometric, person, LR){
  # get the list of graphs from a person
  gList<-numeric()
  for(si in session_range){
    for(ni in number_range){
      session<-session_range[si]
      number<-number_range[ni]
      g<-getGraph(biometric,person,session,number, LR)
      g_tmp<-list(graph=g)
      gList<-rbind(gList, g_tmp)
    }
  }
  gList<-as.data.frame(gList)
  gList
}




############################################################################################################################################
###7. rapidscore_anchor, to be used in register anchor function
### The match score is redefined as the proportion of anchor nodes that match with nodes in g.
## g1 =g, g2=g_anchor
rapidscore_anchor <- function(g1, g2, tol){
  
  #tol = size of bounding box to determine a vertex match
  
  xvals <- get.vertex.attribute(g1, "xcoord")
  yvals <- get.vertex.attribute(g1, "ycoord")
  vertices_enr <- cbind(xvals, yvals)
  num_vertices_enr <- length(vertices_enr[,1])
  
  xvals <- get.vertex.attribute(g2, "xcoord")
  yvals <- get.vertex.attribute(g2, "ycoord")
  vertices_que <- cbind(xvals, yvals)
  num_vertices_que <- length(vertices_que[,1])
  
  match_score <- 0
  vertices_enr_taken <- rep(FALSE, times=num_vertices_enr)
  vertices_que_taken <- rep(FALSE, times=num_vertices_que)
  
  for (i in 1:num_vertices_enr){   
    for (j in 1:num_vertices_que){
      
      if (vertices_enr_taken[i]) next
      if (vertices_que_taken[j]) next
      
      diff <- sqrt((vertices_enr[i,1]-vertices_que[j,1])**2+(vertices_enr[i,2]-vertices_que[j,2])**2)
      
      if (diff<=tol) {
        match_score <- match_score + 1
        vertices_enr_taken[i]<-TRUE
        vertices_que_taken[j]<-TRUE
      }
    }
  }
  match_score <- 1 - (match_score/num_vertices_que)
  match_score
}
############# END rapidscore_anchor - "rapid" **count** of "common" nodes for two graphs g1 & g2 ###############

##### getTopology(g)
#### function that returns basic topology reesults

getTopology<-function(g){
  Size<-network.size(g)
  Edgecount<-network.edgecount(g)
  Edgenoderatio<-network.edgecount(g)/network.size(g)
  cd<-component.dist(g)
  largest_comp<-which(cd$csize==max(cd$csize))[1]
  
  C1<-sort(cd$csize, decreasing=TRUE)[1]
  C2<-sort(cd$csize, decreasing=TRUE)[2]
  Isolated<-cd$cdist[1]
  
  n_c1<-which( cd$membership==largest_comp)
  n_todelete<-setdiff(1:network.size(g), n_c1)
  g_c1<-g
  delete.vertices(g_c1, n_todelete)
  edges_g_c1<-0
  if(network.edgecount(g_c1)>0) C1_l<-sum( get.edge.value(g_c1,"edgelength") )
  C1c2<-(sort(cd$csize, decreasing=TRUE)[1]+sort(cd$csize, decreasing=TRUE)[2])
  ENRc1<-network.edgecount(g_c1)/network.size(g_c1)
  
  deg<-degree(g, gmode="graph", cmode="indegree")
  Dmax<-max(deg)
  
  starSize<-dim(get_star(g))[1]
  twostarSize<-dim(get_twostar(g))[1]  
  
  degree_tmp<-rep(0, times=max(deg) )
  for(j in 1:max(deg)){
    degree_tmp[j]<-length(which(deg==j))
  }
  degree_tmp<-degree_tmp/sum(degree_tmp)

  list(Size=Size, Edgecount=Edgecount, Edgenoderatio=Edgenoderatio, C1=C1, C2=C2, starSize=starSize, twostarSize=twostarSize)
  
}#end of function
