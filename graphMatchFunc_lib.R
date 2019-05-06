
############################################################################################################
getGraphEditPath_nodeBased<-function(g1,g2,cid_n, cid_e,f_edgeinclude){
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
  
}

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
#############################################################################################################
getMCS_nodeBased<- function(g1, g2, edits){
  
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
getGraphEditPath_edgeBased<-function(g1,g2,cid_n, cid_e){
  g1_length<-get.edge.value(g1, "edgelength")
  g1_slope<-tan(get.edge.value(g1, "edgeslope")*pi/180)
  g2_length<-get.edge.value(g2, "edgelength")
  g2_slope<-tan(get.edge.value(g2, "edgeslope")*pi/180)
  EL_1<-as.matrix.network(g1, "edgelist")
  EL_2<-as.matrix.network(g2, "edgelist")
  g1_vertices<-cbind(get.vertex.attribute(g1,"xcoord"),get.vertex.attribute(g1,"ycoord"))
  g2_vertices<-cbind(get.vertex.attribute(g2,"xcoord"),get.vertex.attribute(g2,"ycoord"))
  
  deg1<-degree(g1,gmode="graph" )
  deg2<-degree(g2, gmode="graph")
  
  #create a vector containing sum of the degrees of the end vertices of each edge
  EL_1_deg<-rep(0, times=dim(EL_1)[1])
  for(i in 1:dim(EL_1)[1]){
    node_1<-EL_1[i,1]
    node_2<-EL_1[i,2]
    deg_sum<-deg1[node_1]+deg1[node_2]    
    EL_1_deg[i]<- deg_sum
  }
  
  EL_2_deg<-rep(0, times=dim(EL_2)[1])
  for(i in 1:dim(EL_2)[1]){
    node_1<-EL_2[i,1]
    node_2<-EL_2[i,2]
    deg_sum<-deg2[node_1]+deg2[node_2]   
    EL_2_deg[i]<-deg_sum
  }
  
  N1<-length(g1_length)
  N2<-length(g2_length) 
  
  ##Build the cost matrix
  ############################
  
  block1<-matrix(0,N1, N2)
  
  for (i in 1:N1){
    for(j in 1:N2){
      
     #sij<-sqrt( (g1_length[i]-g2_length[j])^2 + (g1_slope[i]-g2_slope[j])^2 )
      
      g1_n<-EL_1[i,]
      g2_n<-EL_2[j,]
      x1y1<-g1_vertices[g1_n,]
      x1y1_swap<-rbind(x1y1[2,], x1y1[1,])
      x2y2<-g2_vertices[g2_n,]
      
      d1<-sum( sqrt(rowSums((x1y1-x2y2)^2 ) ) )
      d2<-sum( sqrt(rowSums((x1y1_swap-x2y2)^2 ) ) )
      sij<-min(d1,d2)
      
      #w*EL_1_deg[i]-EL_2_deg[j]    # w is a weight to give importance to this parameter
     
      
      block1[i,j]<-sij
    }#end of j loop
  }#end of i loop
  
  block2<-matrix(100000,N1,N1)
  for(i in 1:N1) 
    block2[i,i]<- cid_n+ (cid_e*EL_1_deg[i])
  
  block3<-matrix(100000,N2,N2)
  for(j in 1:N2)
    block3[j,j]<- cid_n + (cid_e*EL_2_deg[j])
  
  block4<-matrix(0,N2,N1)
  
  upper<-cbind(block1, block2)
  lower<-cbind(block3, block4)
  C<-rbind(upper,lower)
  
  x<-solve_LSAP(C, maximum=FALSE)
  edits <- cbind(seq_along(x), x)
  GED <- sum(C[edits])
  
  list(editPath=edits, editCost=GED)

}
##########################################################################################
getMCS_edgeBased<- function(g1, g2, edits){ #edge induced subgraph of destination graph
  
  E <- network.edgecount(g1)
  Q <- network.edgecount(g2)
  
  #edits[1:E,]
  #create a new network for the mcs starting with Q deleting all edges not substituted
  #========================================================================
  
  mcs <- g2
  
  #find all substituted edge from the edits
  edgeToKeep<-numeric()
  for(i in 1:E){
    edge_sub<-edits[i,]
    if(edge_sub[2]<=Q) edgeToKeep<-c(edgeToKeep,edge_sub[2])
    
  }
  
  
  edgeToDel<-setdiff(1:Q,edgeToKeep)
  delete.edges(mcs,edgeToDel)
  
  if(length(edgeToKeep)>0){
    mcs_vert_deg<-degree(mcs, gmode="graph")
    vertToDel<-which(mcs_vert_deg==0)
    delete.vertices(mcs, vertToDel)
    if(network.size(mcs)==0) mcs<-network.initialize(1)
  }
  
  if(length(edgeToKeep)==0){
    mcs<-network.initialize(1)
  }
  
#   ############# optional plotting#################
#   windows()
#   par(mfrow=c(1,2))
#   plotspatialgraph(g2)
#   title("Destination graph")
#   plotspatialgraph(mcs)
#   title(paste("Edge based with cid=(",cid_n,",",cid_e,")", sep=""))
##############################################################################    
  mcs
  
}#end of function





###############################################################################
###############Star-based############################

##########################################################################################################
getGraphEditPath_starBased<-function(g1,g2,cid_n, cid_e){
  g1_vertices<-cbind(get.vertex.attribute(g1,"xcoord"),get.vertex.attribute(g1,"ycoord"))
  g2_vertices<-cbind(get.vertex.attribute(g2,"xcoord"),get.vertex.attribute(g2,"ycoord"))
  
  deg1<-degree(g1,gmode="graph" )
  deg2<-degree(g2, gmode="graph")
  
  g1_star_df<-get_star(g1)
  g2_star_df<-get_star(g2)
  
  N1<-dim(g1_star_df)[1]
  N2<-dim(g2_star_df)[1]
  
  ##Build the cost matrix
  ############################
  
  block1<-matrix(0,N1, N2)
  
  for (i in 1:N1){
    for(j in 1:N2){
      
      star1_c<-g1_star_df$center[[i]]
      star2_c<-g2_star_df$center[[j]]
      
      star1_neighb<-unlist(g1_star_df$neighb[[i]])
      star2_neighb<-unlist(g2_star_df$neighb[[j]])
      
      x1y1<-g1_vertices[c(star1_c,star1_neighb),]
      x2y2<-g2_vertices[c(star2_c,star2_neighb),]
      
      sij<-sum( sqrt(rowSums((x1y1-x2y2)^2 ) ) )
      
     
      block1[i,j]<-sij
    }#end of j loop
  }#end of i loop
  
  block2<-matrix(100000,N1,N1)
  for(i in 1:N1) 
  {
    star1_neighb<-unlist(g1_star_df$neighb[[i]])
    sum_neighb_deg<-sum(deg1[star1_neighb])
    block2[i,i]<- cid_n+ (cid_e*sum_neighb_deg)
  }
     
  block3<-matrix(100000,N2,N2)
  for(j in 1:N2)
  {
    star2_neighb<-unlist(g2_star_df$neighb[[j]])
    sum_neighb_deg<-sum(deg2[star2_neighb])
    block3[j,j]<- cid_n+ (cid_e*sum_neighb_deg)
  }
    
  block4<-matrix(0,N2,N1)
  
  upper<-cbind(block1, block2)
  lower<-cbind(block3, block4)
  C<-rbind(upper,lower)
  
  x<-solve_LSAP(C, maximum=FALSE)
  edits <- cbind(seq_along(x), x)
  GED <- sum(C[edits])
  
  list(editPath=edits, editCost=GED)    
  
}
##########################################################################################
getMCS_starBased<- function(g1, g2, edits){ #edge induced subgraph of destination graph
  
  g1_star_df<-get_star(g1)
  g2_star_df<-get_star(g2)
  
  E<-dim(g1_star_df)[1]
  Q<-dim(g2_star_df)[1]
  
  #edits[1:E,]
  #create a new network for the mcs starting with Q deleting all edges not substituted
  #========================================================================
  
  mcs <- g2
  
  #find all substituted stars from edits. Identify by the center
starToKeep<-numeric()
  for(i in 1:E){
   star_sub<-edits[i,]
    if(star_sub[2]<=Q) starToKeep<-c(starToKeep,star_sub[2])
    
  }
  length(starToKeep)
  
centers<-as.numeric(unlist(g2_star_df$center)[starToKeep] )
  neighbs<-as.numeric(unlist(g2_star_df$neighb[starToKeep])) 
  
nodesToKeep<-union(centers,neighbs)
  nodesToDel<-setdiff( 1:network.size(g2), nodesToKeep)
  delete.vertices(mcs,nodesToDel)
  
  if(network.size(mcs)==0) mcs<-network.initialize(1)
#   ############# optional plotting#################
#   windows()
#   par(mfrow=c(1,2))
#   plotspatialgraph(g2)
#   title("Destination graph")
#   plotspatialgraph(mcs)
#   title(paste("Star based with cid=(",cid_n,",",cid_e,")", sep=""))
#   ##############################################################################    
  mcs
  
}#end of function


#################################################################################

###############Two Star-based############################

##########################################################################################################
getGraphEditPath_twostarBased<-function(g1,g2,cid_n, cid_e){
  g1_vertices<-cbind(get.vertex.attribute(g1,"xcoord"),get.vertex.attribute(g1,"ycoord"))
  g2_vertices<-cbind(get.vertex.attribute(g2,"xcoord"),get.vertex.attribute(g2,"ycoord"))
  
  deg1<-degree(g1,gmode="graph" )
  deg2<-degree(g2, gmode="graph")
  
  g1_twostar_df<-get_twostar(g1)
  g2_twostar_df<-get_twostar(g2)
  
  N1<-dim(g1_twostar_df)[1]
  N2<-dim(g2_twostar_df)[1]
  
  ##Build the cost matrix
  ############################
  
  block1<-matrix(0,N1, N2)
  
  for (i in 1:N1){
    for(j in 1:N2){
      
      g1_star1<-g1_twostar_df$star1[[i]]
      g1_star2<-g1_twostar_df$star2[[i]]
      
      g2_star1<-g2_twostar_df$star1[[j]]
      g2_star2<-g2_twostar_df$star2[[j]]
      
      g1_twostar_vert<-c( as.numeric(g1_star1$center), as.numeric( unlist(g1_star1$neighb) ), as.numeric(g1_star2$center),as.numeric( unlist(g1_star2$neighb)) )
      g2_twostar_vert<-c( as.numeric(g2_star1$center), as.numeric( unlist(g2_star1$neighb) ), as.numeric(g2_star2$center),as.numeric( unlist(g2_star2$neighb)) )
      
          
      x1y1<-g1_vertices[g1_twostar_vert,]
      x2y2<-g2_vertices[g2_twostar_vert,]
      
      sij<-sum( sqrt(rowSums((x1y1-x2y2)^2 ) ) )
      
      
      block1[i,j]<-sij
    }#end of j loop
  }#end of i loop
  
  block2<-matrix(100000,N1,N1)
  for(i in 1:N1) 
  {
    g1_star1<-g1_twostar_df$star1[[i]]
    g1_star2<-g1_twostar_df$star2[[i]]
    
    star1<-c( as.numeric(g1_star1$center), as.numeric( unlist(g1_star1$neighb) ) )
    star2<-c(as.numeric(g1_star2$center),as.numeric( unlist(g1_star2$neighb)) )
    
    neighb_vert<-setdiff(  union(star1,star2), intersect( star1, star2 ) )
       
    sum_neighb_deg<-sum(deg1[neighb_vert])-length(neighb_vert)
    block2[i,i]<- cid_n+ (cid_e*sum_neighb_deg)
  }
  
  block3<-matrix(100000,N2,N2)
  for(j in 1:N2)
  {
    g2_star1<-g2_twostar_df$star1[[j]]
    g2_star2<-g2_twostar_df$star2[[j]]
    
    star1<-c( as.numeric(g2_star1$center), as.numeric( unlist(g2_star1$neighb) ) )
    star2<-c(as.numeric(g2_star2$center),as.numeric( unlist(g2_star2$neighb)) )
    
    neighb_vert<-setdiff(  union(star1,star2), intersect( star1, star2 ) )    
    sum_neighb_deg<-sum(deg2[neighb_vert])-length(neighb_vert)
                 
    block3[j,j]<- cid_n+ (cid_e*sum_neighb_deg)
  }
  
  block4<-matrix(0,N2,N1)
  
  upper<-cbind(block1, block2)
  lower<-cbind(block3, block4)
  C<-rbind(upper,lower)
  
  x<-solve_LSAP(C, maximum=FALSE)
  edits <- cbind(seq_along(x), x)
  GED <- sum(C[edits])
  
  list(editPath=edits, editCost=GED)    
  
}
##########################################################################################
getMCS_twostarBased<- function(g1, g2, edits){ #edge induced subgraph of destination graph
  
  g1_twostar_df<-get_twostar(g1)
  g2_twostar_df<-get_twostar(g2)
  
  E<-dim(g1_twostar_df)[1]
  Q<-dim(g2_twostar_df)[1]
  
  #edits[1:E,]
  #create a new network for the mcs starting with Q deleting all edges not substituted
  #========================================================================
  
  mcs <- g2
  
  #find all substituted twostars from edits. Identify by the center
  twostarToKeep<-numeric()
  for(i in 1:E){
    twostar_sub<-edits[i,]
    if(twostar_sub[2]<=Q) twostarToKeep<-c(twostarToKeep,twostar_sub[2])
    
  }
  if(length(twostarToKeep)>0){
    vertexList<-numeric()
    for(i in 1:length(twostarToKeep)){
      ind<-as.numeric(twostarToKeep[i])
      g2_star1<-g2_twostar_df$star1[[ind]]
      vertexList<-c(vertexList, as.numeric(g2_star1$center))
      vertexList<-c(vertexList, as.numeric(unlist(g2_star1$neighb) ))
      
      g2_star2<-g2_twostar_df$star2[[ind]]
      vertexList<-c(vertexList, as.numeric(g2_star2$center))
      vertexList<-c(vertexList, as.numeric(unlist(g2_star2$neighb) ))
      
    }
    vertexList<-unique(vertexList)
    
    nodesToKeep<-vertexList
    nodesToDel<-setdiff( 1:network.size(g2), nodesToKeep)
    delete.vertices(mcs,nodesToDel)
    
    if(network.size(mcs)==0) mcs<-network.initialize(1)
    ############# optional plotting
  }
  if(length(twostarToKeep)==0){
    mcs<-network.initialize(1)
  }
  
#   #################
#   windows()
#   par(mfrow=c(1,2))
#   plotspatialgraph(g2)
#   title("Destination graph")
#   plotspatialgraph(mcs)
#   title(paste("Two Star based with cid=(",cid_n,",",cid_e,")", sep=""))
#   ##############################################################################    
  mcs
  
}#end of function


#################################################################################