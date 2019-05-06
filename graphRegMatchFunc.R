#Functions defined in this file
#================================
# registerGraphs<-function(g1,g2,tol,f_plot)
# getGraphEditPath<-function(g1,g2,cid, f_edgeinclude)
# getEdgeEditCost<-function(g1, i, g2, j, cid)
# getMCS<-function(g1,g2, edits)

## 1.
################################################################################
## START registerGraphs function.
## This function takes 2 graphs as input, aligns them on the best aligning edge
## Also requires tol, the maximum matching tolerance for graph vertices
## It returns a list of 4 objects: the shifted and rotated graphs corresponding to 
## matching edge pair, the best matching edge pair and the best match score 
## (lower is better as score measures distance)
## It can show the best alignment in a plot, if the flag(f_plot) is set to 1.
####################################################################################

registerGraphs<-function(g1,g2,tol,f_plot){
  
  g1_edgelength<-get.edge.value(g1,"edgelength")
  g2_edgelength<-get.edge.value(g2,"edgelength")
  g1_edgeslope<-get.edge.value(g1,"edgeslope")
  g2_edgeslope<-get.edge.value(g2,"edgeslope")
  
  g1_edgelist<-as.matrix.network(g1, matrix.type="edgelist")
  g2_edgelist<-as.matrix.network(g2, matrix.type="edgelist")
  
  
  #match edges by comparing lengths and slopes
  #===========================================
  
  edge_matching <- matrix(0, length(g1_edgelength), length(g2_edgelength) )
  
  
  for (i in 1:length(edge_matching[,1])){
    for (j in 1:length(edge_matching[1,])){
      
      #      d1<-length(get.neighborhood(g1,g1_edgelist[i,1], type="combined"))+length(get.neighborhood(g1,g1_edgelist[i,2], type="combined"))
      #      d2<-length(get.neighborhood(g2,g2_edgelist[j,1], type="combined"))+length(get.neighborhood(g2,g2_edgelist[j,2], type="combined"))
      #      
      
      #edge_matching[i,j] <- abs(g1_edgelength[i]-g2_edgelength[j])
      
      edge_matching[i,j] <- sqrt((g1_edgelength[i]-g2_edgelength[j])**2+(g1_edgeslope[i]-g2_edgeslope[j])**2)
      edge_matching[i,j] <- edge_matching[i,j]/(0.5*g1_edgelength[i]+0.5*g2_edgelength[j])
      #edge_matching[i,j] <- edge_matching[i,j]/(d1+d2)
      
      
    }
  }
  
  ranking <- order(edge_matching, decreasing=FALSE)
  
  
  
  #take first 10% closest edge pairs to realign and get best match score
  #================================================================
  NEdgePairs<-0.1*length(ranking)
  #NEdgePairs<-length(ranking)
  if(length(ranking)<NEdgePairs) { NEdgePairs<- length(ranking) }
  bestMatchScore<-1
  bestAlignEdge<-rep(0,2)
  g1_shiftrot<-g1
  g2_shiftrot<-g2
  scores<-numeric()
  alignedEdge<-numeric()
  bestRank<-0
  
  
  for(i in 1:NEdgePairs) {
    #for(i in 2:2){
    
    rank<-i
    
    column <- ceiling(ranking[rank]/length(edge_matching[,1]))
    row <- ranking[rank]%%length(edge_matching[,1])
    if (row==0) row <- length(edge_matching[,1])
    edgepair <- c(row, column)
    alignedEdge<-rbind(alignedEdge,edgepair)
    
    #print(edgepair)
    
    g1_temp<-shiftrotate(g1,row)$g_aligned #new, after shiftrot changed to return list
    g2_temp<-shiftrotate(g2,column)$g_aligned
    tempscore<-rapidscore(g1_temp,g2_temp,tol)
    
    scores<-c(scores,tempscore)
    
    
    
    if(tempscore<bestMatchScore){
      bestMatchScore<-tempscore
      bestAlignEdge<-c(row,column)
      g1_shiftrot<-g1_temp
      g2_shiftrot<-g2_temp
      bestRank<-rank
    }
    
    
  } #end of i loop
  
  if(f_plot==1){
    windows()
    plotreggraphs(g1_shiftrot,g2_shiftrot)
    title("registered graphs")
  }
  
  output<-list(graph1=g1_shiftrot,graph2=g2_shiftrot,bestEdgePair=bestAlignEdge,bestScore=bestMatchScore,bestRank=bestRank)
  
}# end of function
########################################################################################################

###################################################################################
########### BEGIN getGraphEditPath ################################################
## This function takes 2 ALIGNED graphs and cid, the insertion/deletion cost as input
## It has a fourth input, a flag that selects if edges must be included in deciding the edit path (1 include, 0 exculde)
## It outputs the cheapest graph edit path in going from g1 to g2
## The cost matrix uses euclidean distance between vertex attributes as the basis 
## for computing cost


getGraphEditPath<-function(g1,g2,cid,f_edgeinclude){
  g1_vertices<-cbind(get.vertex.attribute(g1,"xcoord"),get.vertex.attribute(g1,"ycoord"))
  g2_vertices<-cbind(get.vertex.attribute(g2,"xcoord"),get.vertex.attribute(g2,"ycoord"))
  
  N1<-network.size(g1)
  N2<-network.size(g2)
  
  ##Build the cost matrix
  ############################
  
  block1<-matrix(0,N1, N2)
  
  for (i in 1:N1){
    for(j in 1:N2){
      
      sij<-sqrt( (g1_vertices[i,1]-g2_vertices[j,1])**2+(g1_vertices[i,2]-g2_vertices[j,2])**2 )    
      #   if(sij<=5){sij=0}
      #   if(sij>5){sij=100000}
      
      if(f_edgeinclude){
        edgeedit_cost<-getEdgeEditCost(g1,i,g2,j,cid)
      }
      else{
        edgeedit_cost<-0
      }
      
      block1[i,j]<-sij+edgeedit_cost
    }#end of j loop
  }#end of i loop
  
  block2<-matrix(100000,N1,N1)
  for(i in 1:N1) 
    block2[i,i]<- cid+ (cid*length(get.neighborhood(g1,i)))
  
  block3<-matrix(100000,N2,N2)
  for(j in 1:N2)
    block3[j,j]<- cid + (cid*length(get.neighborhood(g2,j)))
  
  block4<-matrix(0,N2,N1)
  
  upper<-cbind(block1, block2)
  lower<-cbind(block3, block4)
  C<-rbind(upper,lower)
  
  x<-solve_LSAP(C, maximum=FALSE)
  edits <- cbind(seq_along(x), x)
  GED <- sum(C[edits])
  
  list(editPath=edits, editCost=GED)
  
}#end of graph edit path
###################################################################################

getGraphEditPath_twocids<-function(g1,g2,cid_n, cid_e,f_edgeinclude){
  g1_vertices<-cbind(get.vertex.attribute(g1,"xcoord"),get.vertex.attribute(g1,"ycoord"))
  g2_vertices<-cbind(get.vertex.attribute(g2,"xcoord"),get.vertex.attribute(g2,"ycoord"))
  
  N1<-network.size(g1)
  N2<-network.size(g2)
  
  ##Build the cost matrix
  ############################
  
  block1<-matrix(0,N1, N2)
  
  for (i in 1:N1){
    for(j in 1:N2){
      
      sij<-sqrt( (g1_vertices[i,1]-g2_vertices[j,1])**2+(g1_vertices[i,2]-g2_vertices[j,2])**2 )    
      #   if(sij<=5){sij=0}
      #   if(sij>5){sij=100000}
      
      if(f_edgeinclude){
        edgeedit_cost<-getEdgeEditCost(g1,i,g2,j,cid_e)
      }
      else{
        edgeedit_cost<-0
      }
      
      block1[i,j]<-sij+edgeedit_cost
    }#end of j loop
  }#end of i loop
  
  block2<-matrix(100000,N1,N1)
  for(i in 1:N1) 
    block2[i,i]<- cid_n+ (cid_e*length(get.neighborhood(g1,i)))
  
  block3<-matrix(100000,N2,N2)
  for(j in 1:N2)
    block3[j,j]<- cid_n + (cid_e*length(get.neighborhood(g2,j)))
  
  block4<-matrix(0,N2,N1)
  
  upper<-cbind(block1, block2)
  lower<-cbind(block3, block4)
  C<-rbind(upper,lower)
  
  x<-solve_LSAP(C, maximum=FALSE)
  edits <- cbind(seq_along(x), x)
  GED <- sum(C[edits])
  
  list(editPath=edits, editCost=GED)
  
}#end of graph edit path_twocids

########start ###################################################################
## This function calculates the edge cost incurred when a node in g1 is substituted 
## for a node in g2, node in g1 is deleted or a node in g2 is inserted.

getEdgeEditCost<-function(g1, i, g2, j, cid){
  
  xvals <- get.vertex.attribute(g1, "xcoord")
  yvals <- get.vertex.attribute(g1, "ycoord")
  vertices_enr <- cbind(xvals, yvals)
  xvals <- get.vertex.attribute(g2, "xcoord")
  yvals <- get.vertex.attribute(g2, "ycoord")
  vertices_que <- cbind(xvals, yvals)
  
  neighbours1 <- get.neighborhood(g1,i)
  neighbours2 <- get.neighborhood(g2,j)
  
  E <- length(get.neighborhood(g1,i))
  Q <- length(get.neighborhood(g2,j))
  
  extra_cost<-0 #initialise
  
  if (E==0 && Q==0) {extra_cost<-0}
  
  if (E==0 && Q>0) {extra_cost<-cid*Q} # cost of inserting Q edges
  
  if (E>0 && Q==0) {extra_cost<-cid*E} # cost of deleting these edges
  
  if (E>0 && Q>0){
    
    block1 <- array(0, c(E,Q))
    for (i in 1:E){   
      for (j in 1:Q){
        
        n1 <- neighbours1[i]
        n2 <- neighbours2[j]
        sij <- sqrt((vertices_enr[n1,1]-vertices_que[n2,1])**2+(vertices_enr[n1,2]-vertices_que[n2,2])**2)
        #    if(sij<=5){sij<-0}
        #    if(sij>5){sij<-100000}
        block1[i,j]<- sij
        
      }
    }
    
    block2 <- array(100000, c(E,E))
    diag(block2) <- cid
    
    block3 <- array(100000, c(Q,Q))
    diag(block3) <- cid
    
    block4 <- array(0, c(Q,E))
    
    upper <- cbind(block1, block2)
    lower <- cbind(block3, block4)
    C <- rbind(upper, lower)
    
    z <- solve_LSAP(C, maximum=FALSE)
    
    extra_cost <- sum(C[cbind(seq_along(z), z)])
  }
  
  extra_cost
  
  
  
}#end of function

#####################################################################################################
########### Edit Path Using bounding box cost function


getGraphEditPath_bb<-function(g1,g2,cid,tol,f_edgeinclude){
  g1_vertices<-cbind(get.vertex.attribute(g1,"xcoord"),get.vertex.attribute(g1,"ycoord"))
  g2_vertices<-cbind(get.vertex.attribute(g2,"xcoord"),get.vertex.attribute(g2,"ycoord"))
  
  N1<-network.size(g1)
  N2<-network.size(g2)
  
  ##Build the cost matrix
  ############################
  
  block1<-matrix(0,N1, N2)
  
  for (i in 1:N1){
    for(j in 1:N2){
      
      sij<-sqrt( (g1_vertices[i,1]-g2_vertices[j,1])**2+(g1_vertices[i,2]-g2_vertices[j,2])**2 )    
      if(sij<=tol){sij=0}
      if(sij>tol){sij=100000}
      
      if(f_edgeinclude){
        edgeedit_cost<-getEdgeEditCost_bb(g1,i,g2,j,cid,tol)
      }
      else{
        edgeedit_cost<-0
      }
      
      block1[i,j]<-sij+edgeedit_cost
    }#end of j loop
  }#end of i loop
  
  block2<-matrix(100000,N1,N1)
  for(i in 1:N1) 
    block2[i,i]<- cid+ (cid*length(get.neighborhood(g1,i)))
  
  block3<-matrix(100000,N2,N2)
  for(j in 1:N2)
    block3[j,j]<- cid + (cid*length(get.neighborhood(g2,j)))
  
  block4<-matrix(0,N2,N1)
  
  upper<-cbind(block1, block2)
  lower<-cbind(block3, block4)
  C<-rbind(upper,lower)
  
  x<-solve_LSAP(C, maximum=FALSE)
  edits <- cbind(seq_along(x), x)
  GED <- sum(C[edits])
  
  list(editPath=edits, editCost=GED)
  
}#end of graph edit path
###################################################################################
########start ###################################################################
## This function calculates the edge cost incurred when a node in g1 is substituted 
## for a node in g2, node in g1 is deleted or a node in g2 is inserted.

getEdgeEditCost_bb<-function(g1, i, g2, j, cid, tol){
  
  xvals <- get.vertex.attribute(g1, "xcoord")
  yvals <- get.vertex.attribute(g1, "ycoord")
  vertices_enr <- cbind(xvals, yvals)
  xvals <- get.vertex.attribute(g2, "xcoord")
  yvals <- get.vertex.attribute(g2, "ycoord")
  vertices_que <- cbind(xvals, yvals)
  
  neighbours1 <- get.neighborhood(g1,i)
  neighbours2 <- get.neighborhood(g2,j)
  
  E <- length(get.neighborhood(g1,i))
  Q <- length(get.neighborhood(g2,j))
  
  extra_cost<-0 #initialise
  
  if (E==0 && Q==0) {extra_cost<-0}
  
  if (E==0 && Q>0) {extra_cost<-cid*Q} # cost of inserting Q edges
  
  if (E>0 && Q==0) {extra_cost<-cid*E} # cost of deleting these edges
  
  if (E>0 && Q>0){
    
    block1 <- array(0, c(E,Q))
    for (i in 1:E){   
      for (j in 1:Q){
        
        n1 <- neighbours1[i]
        n2 <- neighbours2[j]
        sij <- sqrt((vertices_enr[n1,1]-vertices_que[n2,1])**2+(vertices_enr[n1,2]-vertices_que[n2,2])**2)
        if(sij<=tol){sij<-0}
        if(sij>tol){sij<-100000}
        block1[i,j]<- sij
        
      }
    }
    
    block2 <- array(100000, c(E,E))
    diag(block2) <- cid
    
    block3 <- array(100000, c(Q,Q))
    diag(block3) <- cid
    
    block4 <- array(0, c(Q,E))
    
    upper <- cbind(block1, block2)
    lower <- cbind(block3, block4)
    C <- rbind(upper, lower)
    
    z <- solve_LSAP(C, maximum=FALSE)
    
    extra_cost <- sum(C[cbind(seq_along(z), z)])
  }
  
  extra_cost
  
  
  
}#end of function







######################################################################################################
############# START getmcs - takes an edit path and returns an mcs from two graphs ###################
## MCS is defined as a node induced subgraph of the query graph g2. 
## The MCS will have edges between the 2 nodes if corresponding nodes in g1 had an edge between them.
## inputs: g1 and g2 are registered graphs, edits are the list of edit operations from the 
## function getGraphEditPath()
#######################################################################################################

getMCS <- function(g1, g2, edits){
  
  E <- network.size(g1)
  Q <- network.size(g2)
  
  
  #create a new network for the mcs, with no edges, and same size as query
  #========================================================================
  
  mcs <- network.initialize(Q, directed=FALSE)
  
  
  # now consider all edges in g1 and ADD an edge to the mcs if the corresponding nodes 
  # in g2 are substituted as per the optimal assignment
  #============================================================
  
  edges <- as.matrix(g1, matrix.type="edgelist")
  num_edges <- dim(edges)[1]
  starts <- numeric()
  ends <- numeric()
  
  for (i in 1:num_edges) {
    edge <- edges[i,]
    base <- as.integer(edits[edge[1],2])
    tip <- as.integer(edits[edge[2],2])
    
    if (base<Q & tip<Q){
      
      possibles <- get.neighborhood(g2,base)
      degree <- length(possibles)
      
      if (degree>0) {
        
        flag <- 0
        for (j in 1:degree){
          if (possibles[j]==tip) {add.edge(mcs, base, tip)}
        }
        
      }
    }
  }
  
  
  #paste in vertex attributes from g2 into mcs
  #================================================
  
  
  xvals <- get.vertex.attribute(g2, "xcoord")
  yvals <- get.vertex.attribute(g2, "ycoord")
  vertices <- cbind(xvals, yvals)
  
  set.vertex.attribute(mcs, "xcoord", vertices[,1])
  set.vertex.attribute(mcs, "ycoord", vertices[,2])
  
  
  
  # delete from the mcs, nodes of g2 that have to be 
  # inserted as per the optimal assignment
  #====================================================
  
  vertices_to_delete <- numeric()
  
  for (i in 1:Q){
    
    r <- E+i
    fate <-edits[r,2]
    if (fate<=Q) {vertices_to_delete <- c(vertices_to_delete, i)}
  }
  
  delete.vertices(mcs, vertices_to_delete)
  
  if(network.size(mcs) > 0 ){
    #add edge and node attributes - this part newly added on 10/07/2013
    mcs_edges<-as.matrix.network(mcs, matrix.type="edgelist")
    if(length(mcs_edges[,1])>0){
      
      lengths<-rep(0,times=length(mcs_edges[,1]) )
      slopes<-rep(0,times=length(mcs_edges[,1]) )
      v<-cbind(get.vertex.attribute(mcs, "xcoord") , get.vertex.attribute(mcs, "ycoord") )
      
      for (i in 1:length(mcs_edges[,1])){
        
        start <- mcs_edges[i,1]
        end <- mcs_edges[i,2]
        
        lengths[i] <- sqrt((v[start,1]-v[end,1])**2+(v[start,2]-v[end,2])**2)
        
        if((v[end,1]-v[start,1])==0){slopes[i]=90}
        if((v[end,1]-v[start,1])!=0){
          if (v[start,1]<v[end,1]) slopes[i] <- (180/pi)*atan((v[end,2]-v[start,2])/(v[end,1]-v[start,1]))
          if (v[start,1]>=v[end,1]) slopes[i] <- (180/pi)*atan((v[start,2]-v[end,2])/(v[start,1]-v[end,1]))
          
        }
        
        lengths[i]<-round(lengths[i],2)
        slopes[i]<-round(slopes[i],2)
        #print(i)
      }
      
      set.edge.attribute(mcs,"edgelength",lengths)
      set.edge.attribute(mcs,"edgeslope",slopes)
      
      
      
    }
  }
  
  
  mcs
  
}#end of function
##########################################################################################################
