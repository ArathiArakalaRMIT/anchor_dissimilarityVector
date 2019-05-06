library(doMC)
## 1.
## This function finds the anchor inside a given graph
## g1= g and g2 = g_anchor
## 
## Edge pairs are categorised as long and short and compared.
## The best 'L' matching edge pairs are used to translate n rotate
## the graphs n count the proportion of anchor nodes that matched.

registerAnchor<-function(g1=g,g2=g_anchor,tol,maxtheta, f_plot, L){
  #L is the hard number of edgepairs to consider
  
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
      
       
      edge_matching[i,j] <- abs(g1_edgelength[i]-g2_edgelength[j])
      
    #  edge_matching[i,j] <- sqrt((g1_edgelength[i]-g2_edgelength[j])**2+(g1_edgeslope[i]-g2_edgeslope[j])**2)
    # edge_matching[i,j] <- abs( (g1_edgelength[i]-g2_edgelength[j]) ) + abs( tan( (g1_edgeslope[i]-g2_edgeslope[j])*pi/2  ) )
 
      
    }
  }
  
  ranking <- order(edge_matching, decreasing=FALSE)
  
  #########added code to test#############
  #r<-ranking[1: (0.02*length(ranking))]
  
  
  ##########################################
  
  #select all edge pairs with average length greater than 10. Put them in 1 structure
  #select all other edge pairs and put them in another structure
  edgePairs_long<-numeric() 
  edgePairs_short<-numeric()
  longTh<-mean(c(median(g1_edgelength),median(g2_edgelength) ) )
  
  for(i in 1:length(ranking)){
  #  for(i in r){  
    rank<-i
    
    edge2 <- ceiling(ranking[rank]/length(edge_matching[,1]))
    edge1 <- ranking[rank]%%length(edge_matching[,1])
    if (edge1==0) edge1 <- length(edge_matching[,1])
    tmp<-list(e1=edge1, e2=edge2, rank=rank)
    
    if( g1_edgelength[edge1]>= longTh && g2_edgelength[edge2] >= longTh ) edgePairs_long<-rbind(edgePairs_long, tmp)
    if( g1_edgelength[edge1] < longTh || g2_edgelength[edge2]< longTh ) edgePairs_short<-rbind(edgePairs_short, tmp)    
  }
  edgePairs_long<-as.data.frame(edgePairs_long)
  edgePairs_short<-as.data.frame(edgePairs_short)
    
#   NEdgePairs_long<-max(1, round(0.02*length(edgePairs_long$e1)) )
#   
#   NEdgePairs_short<-max(1, round(0.02*length(edgePairs_short$e1)) )
  
    NEdgePairs_long<-min(length(edgePairs_long$e1), L/2 )
    NEdgePairs_short<-min(length(edgePairs_short$e1), L/2 )
  
  bestMatchScore<-numeric()
  bestAlignEdge<-numeric()
  #   g1_shiftrot<-g1
  #   g2_shiftrot<-g2
  g1_shiftrot<-numeric()
  g2_shiftrot<-numeric()
  scores<-numeric()
  alignedEdge<-numeric()
  bestRank<-numeric()
  bestIndex<-numeric()
  bestAnchorScore<-numeric()
  outputlist_long<-numeric()
  outputlist_short<-numeric()
  
 
  if(NEdgePairs_long>0){
    #registerDoMC(8)
    #outputlist_long<-foreach(i=1:NEdgePairs_long, .combine="cbind", .packages=c("network","sna", "clue")  )%dopar%{
    outputlist_long<-foreach(i=1:NEdgePairs_long, .combine="cbind", .packages=c("network","sna", "clue")  )%do%{
      source("basicFunctions.R")
      
      # for(i in 1:length(edgePairs_long$e1)){
      
      row<-edgePairs_long$e1[[i]]
      column<-edgePairs_long$e2[[i]]
      rank<-edgePairs_long$rank[[i]]
      edgepair<-c(row,column)
      
      g1_temp<-shiftrotate(g1,row)
      g2_temp<-shiftrotate(g2,column)
      tempscore<-rapidscore_anchor(g1_temp$g_aligned,g2_temp$g_aligned,tol)
      
      #     windows()
      #     plotreggraphs(g1_temp$g_aligned,g2_temp$g_aligned)
      #     title(paste("Long edge Rank#",rank," with rapidscore =", round(tempscore,2)))
      #     
      op<-list(score=tempscore, alignedEdgePair=edgepair, rank=rank, g1_shiftrot=g1_temp, g2_shiftrot=g2_temp)
      op
      
    } # end of foreach loop
    
    outputlist_long<-as.data.frame(outputlist_long)
    
    
  }
   
 
  #outputlist_short<-foreach(i=1:NEdgePairs_short, .combine="cbind", .packages=c("network","sna", "clue")  )%dopar%{
    outputlist_short<-foreach(i=1:NEdgePairs_short, .combine="cbind", .packages=c("network","sna", "clue")  )%do%{
      
          source("basicFunctions.R")
    
    row<-edgePairs_short$e1[[i]]
    column<-edgePairs_short$e2[[i]]
    rank<-edgePairs_short$rank[[i]]
    edgepair<-c(row, column)
      
    g1_temp<-shiftrotate(g1,row)
    g2_temp<-shiftrotate(g2,column)
    tempscore<-rapidscore_anchor(g1_temp$g_aligned,g2_temp$g_aligned,tol)
#         windows()
#         plotreggraphs(g1_temp$g_aligned,g2_temp$g_aligned)
#         title(paste("Short edge Rank#",rank," with rapidscore =", round(tempscore,2)))
        
        op<-list(score=tempscore, alignedEdgePair=edgepair, rank=rank, g1_shiftrot=g1_temp, g2_shiftrot=g2_temp)
        op
    
  } # end of foreach loop
  
  outputlist_short<-as.data.frame(outputlist_short)
  
  if(NEdgePairs_long>0)
    outputlist<-cbind(outputlist_long,outputlist_short)
  if(NEdgePairs_long==0)
    outputlist<-outputlist_short
  
  #find lowest score here by converting output to dataframe, geting lowest score and corresponding features
  scorelist<-numeric()
  
  for(i in 1:dim(outputlist)[2]){
    
    scorelist<-c(scorelist,outputlist[[i]]$score)   
    
  }
  
  theta1_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g1_shiftrot$theta * (180)/pi 
    theta1_list<-c(theta1_list,op)
  }
  
  
  theta2_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g2_shiftrot$theta * (180)/pi
    theta2_list<-c(theta2_list,op)
  }
  
  orderIndex<-order(scorelist, decreasing=FALSE)
  thetaDiff<-theta1_list-theta2_list
  thetaSatisfy<-which(abs(thetaDiff) < maxtheta )
  
  orderIndexThetaSatisfy<-which(orderIndex %in% thetaSatisfy) #the ordered indices of orderIndex, that satisfy theta contraint
  filteredOrderIndex<-orderIndex[orderIndexThetaSatisfy] #index of edgepairs that satisfy theta contsraint, ordered in increasing order of score.
  lowestScore<-outputlist[[filteredOrderIndex[1]]]$score
  
  minIndices<-which(scorelist[filteredOrderIndex]==lowestScore) #if there are more than 1 best match score, get all of the associated edges, graphs and ranks
  
  for(i in minIndices){
    bestIndex<- rbind(bestIndex,filteredOrderIndex[i])
    bestMatchScore<-rbind(bestMatchScore,outputlist[[filteredOrderIndex[i]]]$score)
    bestAlignEdge<-rbind(bestAlignEdge,outputlist[[filteredOrderIndex[i]]]$alignedEdgePair)
    g1_shiftrot<-rbind(g1_shiftrot,outputlist[[filteredOrderIndex[i]]]$g1_shiftrot)
    g2_shiftrot<-rbind(g2_shiftrot,outputlist[[filteredOrderIndex[i]]]$g2_shiftrot)
    bestRank<-rbind(bestRank,outputlist[[filteredOrderIndex[i]]]$rank)    
 #   bestAnchorScore<-rbind(bestAnchorScore, rapidscore_anchor(outputlist[[filteredOrderIndex[i]]]$g1_shiftrot, outputlist[[filteredOrderIndex[i]]]$g2_shiftrot, tol ))
  }
  
  
  output<-data.frame(graph1=g1_shiftrot, graph2=g2_shiftrot, bestEdgePair=bestAlignEdge, bestScore=bestMatchScore, bestRank=bestRank) 
  
  
  
  if(f_plot==1){
    for(i in 1:length(minIndices)){
      quartz()
      plotspatialgraph(output$graph1.g_aligned[[i]])
      quartz()
      plotspatialgraph(output$graph2.g_aligned[[i]])
      quartz()
      plotreggraphs(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]])
      #plotreggraphs_t(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]]) #for SNIR
      titlestring<-paste("Alignment on rank", i," best edge pair (# ",bestIndex[i], "), d=", round(bestMatchScore,2) , sep="")
      #title(titlestring)
      
    }#end of i loop
    
  }
  
  return(output)
  
}# end of function
# ########################################################################################################


registerGraphs_v1<-function(g1,g2,tol,maxtheta, f_plot, L){
  #L is the hard number of edgepairs to consider
  
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
      
      
      edge_matching[i,j] <- abs(g1_edgelength[i]-g2_edgelength[j])
      
      #  edge_matching[i,j] <- sqrt((g1_edgelength[i]-g2_edgelength[j])**2+(g1_edgeslope[i]-g2_edgeslope[j])**2)
      # edge_matching[i,j] <- abs( (g1_edgelength[i]-g2_edgelength[j]) ) + abs( tan( (g1_edgeslope[i]-g2_edgeslope[j])*pi/2  ) )
      
      
    }
  }
  
  ranking <- order(edge_matching, decreasing=FALSE)
  
  #########added code to test#############
  #r<-ranking[1: (0.02*length(ranking))]
  
  
  ##########################################
  
  #select all edge pairs with average length greater than 10. Put them in 1 structure
  #select all other edge pairs and put them in another structure
  edgePairs_long<-numeric() 
  edgePairs_short<-numeric()
  longTh<-mean(c(median(g1_edgelength),median(g2_edgelength) ) )
  
  for(i in 1:length(ranking)){
    #  for(i in r){  
    rank<-i
    
    edge2 <- ceiling(ranking[rank]/length(edge_matching[,1]))
    edge1 <- ranking[rank]%%length(edge_matching[,1])
    if (edge1==0) edge1 <- length(edge_matching[,1])
    tmp<-list(e1=edge1, e2=edge2, rank=rank)
    
    if( g1_edgelength[edge1]>= longTh && g2_edgelength[edge2] >= longTh ) edgePairs_long<-rbind(edgePairs_long, tmp)
    if( g1_edgelength[edge1] < longTh || g2_edgelength[edge2]< longTh ) edgePairs_short<-rbind(edgePairs_short, tmp)    
  }
  edgePairs_long<-as.data.frame(edgePairs_long)
  edgePairs_short<-as.data.frame(edgePairs_short)
  
  #   NEdgePairs_long<-max(1, round(0.02*length(edgePairs_long$e1)) )
  #   
  #   NEdgePairs_short<-max(1, round(0.02*length(edgePairs_short$e1)) )
  
  NEdgePairs_long<-min(length(edgePairs_long$e1), L/2 )
  NEdgePairs_short<-min(length(edgePairs_short$e1), L/2 )
  
  bestMatchScore<-numeric()
  bestAlignEdge<-numeric()
  #   g1_shiftrot<-g1
  #   g2_shiftrot<-g2
  g1_shiftrot<-numeric()
  g2_shiftrot<-numeric()
  scores<-numeric()
  alignedEdge<-numeric()
  bestRank<-numeric()
  bestIndex<-numeric()
  
  outputlist_long<-numeric()
  outputlist_short<-numeric()
  
  
  if(NEdgePairs_long>0){
    #registerDoMC(8)
    #outputlist_long<-foreach(i=1:NEdgePairs_long, .combine="cbind", .packages=c("network","sna", "clue")  )%dopar%{
    outputlist_long<-foreach(i=1:NEdgePairs_long, .combine="cbind", .packages=c("network","sna", "clue")  )%do%{
      source("basicFunctions.R")
      
      # for(i in 1:length(edgePairs_long$e1)){
      
      row<-edgePairs_long$e1[[i]]
      column<-edgePairs_long$e2[[i]]
      rank<-edgePairs_long$rank[[i]]
      edgepair<-c(row,column)
      
      g1_temp<-shiftrotate(g1,row)
      g2_temp<-shiftrotate(g2,column)
      tempscore<-rapidscore(g1_temp$g_aligned,g2_temp$g_aligned,tol)
      
      #     windows()
      #     plotreggraphs(g1_temp$g_aligned,g2_temp$g_aligned)
      #     title(paste("Long edge Rank#",rank," with rapidscore =", round(tempscore,2)))
      #     
      op<-list(score=tempscore, alignedEdgePair=edgepair, rank=rank, g1_shiftrot=g1_temp, g2_shiftrot=g2_temp)
      op
      
    } # end of foreach loop
    
    outputlist_long<-as.data.frame(outputlist_long)
    
  }
   
  
  
  #outputlist_short<-foreach(i=1:NEdgePairs_short, .combine="cbind", .packages=c("network","sna", "clue")  )%dopar%{
  outputlist_short<-foreach(i=1:NEdgePairs_short, .combine="cbind", .packages=c("network","sna", "clue")  )%do%{
    
    source("basicFunctions.R")
    
    row<-edgePairs_short$e1[[i]]
    column<-edgePairs_short$e2[[i]]
    rank<-edgePairs_short$rank[[i]]
    edgepair<-c(row, column)
    
    g1_temp<-shiftrotate(g1,row)
    g2_temp<-shiftrotate(g2,column)
    tempscore<-rapidscore(g1_temp$g_aligned,g2_temp$g_aligned,tol)
    #         windows()
    #         plotreggraphs(g1_temp$g_aligned,g2_temp$g_aligned)
    #         title(paste("Short edge Rank#",rank," with rapidscore =", round(tempscore,2)))
    
    op<-list(score=tempscore, alignedEdgePair=edgepair, rank=rank, g1_shiftrot=g1_temp, g2_shiftrot=g2_temp)
    op
    
  } # end of foreach loop
  
  outputlist_short<-as.data.frame(outputlist_short)
  
  if(NEdgePairs_long>0)
    outputlist<-cbind(outputlist_long,outputlist_short)
  if(NEdgePairs_long==0)
    outputlist<-outputlist_short
  
  #find lowest score here by converting output to dataframe, geting lowest score and corresponding features
  scorelist<-numeric()
  
  for(i in 1:dim(outputlist)[2]){
    
    scorelist<-c(scorelist,outputlist[[i]]$score)   
    
  }
  
  theta1_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g1_shiftrot$theta * (180)/pi 
    theta1_list<-c(theta1_list,op)
  }
  
  
  theta2_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g2_shiftrot$theta * (180)/pi
    theta2_list<-c(theta2_list,op)
  }
  
  orderIndex<-order(scorelist, decreasing=FALSE)
  thetaDiff<-theta1_list-theta2_list
  thetaSatisfy<-which(abs(thetaDiff) < maxtheta )
  
  orderIndexThetaSatisfy<-which(orderIndex %in% thetaSatisfy) #the ordered indices of orderIndex, that satisfy theta contraint
  filteredOrderIndex<-orderIndex[orderIndexThetaSatisfy] #index of edgepairs that satisfy theta contsraint, ordered in increasing order of score.
  lowestScore<-outputlist[[filteredOrderIndex[1]]]$score
  
  minIndices<-which(scorelist[filteredOrderIndex]==lowestScore) #if there are more than 1 best match score, get all of the associated edges, graphs and ranks
  
  for(i in minIndices){
    bestIndex<- rbind(bestIndex,filteredOrderIndex[i])
    bestMatchScore<-rbind(bestMatchScore,outputlist[[filteredOrderIndex[i]]]$score)
    bestAlignEdge<-rbind(bestAlignEdge,outputlist[[filteredOrderIndex[i]]]$alignedEdgePair)
    g1_shiftrot<-rbind(g1_shiftrot,outputlist[[filteredOrderIndex[i]]]$g1_shiftrot)
    g2_shiftrot<-rbind(g2_shiftrot,outputlist[[filteredOrderIndex[i]]]$g2_shiftrot)
    bestRank<-rbind(bestRank,outputlist[[filteredOrderIndex[i]]]$rank)    
    
  }
  
  
  output<-data.frame(graph1=g1_shiftrot, graph2=g2_shiftrot, bestEdgePair=bestAlignEdge, bestScore=bestMatchScore, bestRank=bestRank) 
  
  
  
  if(f_plot==1){
    for(i in 1:length(minIndices)){
      quartz()
      plotspatialgraph(output$graph1.g_aligned[[i]])
      quartz()
      plotspatialgraph(output$graph2.g_aligned[[i]])
      quartz()
      plotreggraphs(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]])
      #plotreggraphs_t(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]]) #for SNIR
      titlestring<-paste("Alignment on rank", i," best edge pair (# ",bestIndex[i], "), d=", round(bestMatchScore,2) , sep="")
      title(titlestring)
      
    }#end of i loop
    
  }
  
  output
  
}# end of function
# ########################################################################################################
# 











# 
# registerGraphs_v2<-function(g1,g2,tol,maxtheta, f_plot){
#   
#   g1_edgelength<-get.edge.value(g1,"edgelength")
#   g2_edgelength<-get.edge.value(g2,"edgelength")
#   g1_edgeslope<-get.edge.value(g1,"edgeslope")
#   g2_edgeslope<-get.edge.value(g2,"edgeslope")
#   
#   g1_edgelist<-as.matrix.network(g1, matrix.type="edgelist")
#   g2_edgelist<-as.matrix.network(g2, matrix.type="edgelist")
#   
#   
#   #match edges by comparing lengths and slopes
#   #===========================================
#   
#   edge_matching <- matrix(0, length(g1_edgelength), length(g2_edgelength) )
#   
#   for (i in 1:length(edge_matching[,1])){
#     for (j in 1:length(edge_matching[1,])){
#      i_n1<-get.neighborhood(g1,g1_edgelist[i,1], type="combined")
#      i_n2<-get.neighborhood(g1,g1_edgelist[i,2], type="combined")
#      j_n1<-get.neighborhood(g2,g2_edgelist[j,1], type="combined")
#      j_n2<-get.neighborhood(g2,g2_edgelist[j,2], type="combined")
#     
#      deg_i<-length(i_n1)+length(i_n2)
#      deg_j<-length(j_n1)+length(j_n2)
#      
#      n1<-numeric()
#      for(k in 1:length(i_n1)){
#        start<-g1_edgelist[i,1]
#        end<-i_n1[k]
#       tmp<-union( intersect( which (g1_edgelist[,1]==start) , which(g1_edgelist[,2]==end) ) , intersect( which(g1_edgelist[,1]==end), which(g1_edgelist[,2]==start) ) )
#        n1<-c(n1,tmp)
#      }
#      
#      n2<-numeric()
#      for(k in 1:length(i_n2)){
#        start<-g1_edgelist[i,2]
#        end<-i_n2[k]
#        tmp<-union( intersect( which (g1_edgelist[,1]==start) , which(g1_edgelist[,2]==end) ) , intersect( which(g1_edgelist[,1]==end), which(g1_edgelist[,2]==start) ) )
#        n2<-c(n2,tmp)
#      }
#       
#      i_sumedgelengths<-sum(g1_edgelength[n1])+sum(g1_edgelength[n2])-g1_edgelength[i] # the actual edge length is added twice, one from each vertex
#      
#      ############
#      
#      n1<-numeric()
#      for(k in 1:length(j_n1)){
#        start<-g2_edgelist[j,1]
#        end<-j_n1[k]
#        tmp<-union( intersect( which (g2_edgelist[,1]==start) , which(g2_edgelist[,2]==end) ) , intersect( which(g2_edgelist[,1]==end), which(g2_edgelist[,2]==start) ) )
#        n1<-c(n1,tmp)
#      }
#      
#      n2<-numeric()
#      for(k in 1:length(i_n2)){
#        start<-g2_edgelist[j,2]
#        end<-j_n2[k]
#        tmp<-union( intersect( which (g2_edgelist[,1]==start) , which(g2_edgelist[,2]==end) ) , intersect( which(g2_edgelist[,1]==end), which(g2_edgelist[,2]==start) ) )
#        n2<-c(n2,tmp)
#      }
#      
#      j_sumedgelengths<-sum(g2_edgelength[n1])+sum(g2_edgelength[n2])-g2_edgelength[j] # the actual edge length is added twice, one from each vertex
#      edge_matching[i,j] <- abs(i_sumedgelengths-j_sumedgelengths)
#       
#      
#     }
#   }
#   
#   ranking <- order(edge_matching, decreasing=FALSE)
#   
#   #select all edge pairs with average length greater than 10. Put them in 1 structure
#   #select all other edge pairs and put them in another structure
#   edgePairs_long<-numeric() 
#   edgePairs_short<-numeric()
#   for(i in 1:length(ranking)){
#     rank<-i
#     
#     edge2 <- ceiling(ranking[rank]/length(edge_matching[,1]))
#     edge1 <- ranking[rank]%%length(edge_matching[,1])
#     if (edge1==0) edge1 <- length(edge_matching[,1])
#     tmp<-list(e1=edge1, e2=edge2, rank=rank)
#     
#     if( g1_edgelength[edge1]>= 100 && g2_edgelength[edge2] >= 100 ) edgePairs_long<-rbind(edgePairs_long, tmp)
#     if( g1_edgelength[edge1] < 100 || g2_edgelength[edge2]< 100 ) edgePairs_short<-rbind(edgePairs_short, tmp)    
#   }
#   edgePairs_long<-as.data.frame(edgePairs_long)
#   edgePairs_short<-as.data.frame(edgePairs_short)
#   
#   NEdgePairs_long<-min(length(edgePairs_long$e1), round(0.005*length(ranking)) )
#   NEdgePairs_short<-min(length(edgePairs_short$e1), round(0.005*length(ranking)) )
#   
#   bestMatchScore<-numeric()
#   bestAlignEdge<-numeric()
#   #   g1_shiftrot<-g1
#   #   g2_shiftrot<-g2
#   g1_shiftrot<-numeric()
#   g2_shiftrot<-numeric()
#   scores<-numeric()
#   alignedEdge<-numeric()
#   bestRank<-numeric()
#   bestIndex<-numeric()
#   
#   
#   
#   #registerDoMC(8)
#   outputlist_long<-foreach(i=1:NEdgePairs_long, .combine="cbind", .packages=c("network","sna", "clue")  )%dopar%{
#     #      # outputlist<-foreach(i=1:NEdgePairs, .combine="cbind", .packages=c("network","sna", "clue")  )%do%{
#     source("basicGraphFunctions_vpac.R")
#     
#     # for(i in 1:length(edgePairs_long$e1)){
#     
#     row<-edgePairs_long$e1[[i]]
#     column<-edgePairs_long$e2[[i]]
#     rank<-edgePairs_long$rank[[i]]
#     edgepair<-c(row,column)
#     
#     g1_temp<-shiftrotate(g1,row)
#     g2_temp<-shiftrotate(g2,column)
#     tempscore<-rapidscore(g1_temp$g_aligned,g2_temp$g_aligned,tol)
#     
#     #     windows()
#     #     plotreggraphs(g1_temp$g_aligned,g2_temp$g_aligned)
#     #     title(paste("Long edge Rank#",rank," with rapidscore =", round(tempscore,2)))
#     #     
#     op<-list(score=tempscore, alignedEdgePair=edgepair, rank=rank, g1_shiftrot=g1_temp, g2_shiftrot=g2_temp)
#     op
#     
#   } # end of foreach loop
#   
#   outputlist_long<-as.data.frame(outputlist_long)
#   
#   
#   
#   outputlist_short<-foreach(i=1:NEdgePairs_short, .combine="cbind", .packages=c("network","sna", "clue")  )%dopar%{
#     source("basicGraphFunctions_vpac.R")
#     
#     row<-edgePairs_short$e1[[i]]
#     column<-edgePairs_short$e2[[i]]
#     rank<-edgePairs_short$rank[[i]]
#     edgepair<-c(row, column)
#     
#     g1_temp<-shiftrotate(g1,row)
#     g2_temp<-shiftrotate(g2,column)
#     tempscore<-rapidscore(g1_temp$g_aligned,g2_temp$g_aligned,tol)
#     #         windows()
#     #         plotreggraphs(g1_temp$g_aligned,g2_temp$g_aligned)
#     #         title(paste("Short edge Rank#",rank," with rapidscore =", round(tempscore,2)))
#     
#     op<-list(score=tempscore, alignedEdgePair=edgepair, rank=rank, g1_shiftrot=g1_temp, g2_shiftrot=g2_temp)
#     op
#     
#   } # end of foreach loop
#   
#   outputlist_short<-as.data.frame(outputlist_short)
#   
#   outputlist<-cbind(outputlist_long,outputlist_short)
#   
#   #find lowest score here by converting output to dataframe, geting lowest score and corresponding features
#   scorelist<-numeric()
#   
#   for(i in 1:dim(outputlist)[2]){
#     
#     scorelist<-c(scorelist,outputlist[[i]]$score)   
#     
#   }
#   
#   theta1_list<-numeric()
#   for(i in 1:dim(outputlist)[2]){
#     op<-outputlist[[i]]$g1_shiftrot$theta * (180)/pi 
#     theta1_list<-c(theta1_list,op)
#   }
#   
#   
#   theta2_list<-numeric()
#   for(i in 1:dim(outputlist)[2]){
#     op<-outputlist[[i]]$g2_shiftrot$theta * (180)/pi
#     theta2_list<-c(theta2_list,op)
#   }
#   
#   orderIndex<-order(scorelist, decreasing=FALSE)
#   thetaDiff<-theta1_list-theta2_list
#   thetaSatisfy<-which(abs(thetaDiff) < maxtheta )
#   
#   orderIndexThetaSatisfy<-which(orderIndex %in% thetaSatisfy) #the ordered indices of orderIndex, that satisfy theta contraint
#   filteredOrderIndex<-orderIndex[orderIndexThetaSatisfy] #index of edgepairs that satisfy theta contsraint, ordered in increasing order of score.
#   lowestScore<-outputlist[[filteredOrderIndex[1]]]$score
#   
#   minIndices<-which(scorelist[filteredOrderIndex]==lowestScore) #if there are more than 1 best match score, get all of the associated edges, graphs and ranks
#   
#   for(i in minIndices){
#     bestIndex<- rbind(bestIndex,filteredOrderIndex[i])
#     bestMatchScore<-rbind(bestMatchScore,outputlist[[filteredOrderIndex[i]]]$score)
#     bestAlignEdge<-rbind(bestAlignEdge,outputlist[[filteredOrderIndex[i]]]$alignedEdgePair)
#     g1_shiftrot<-rbind(g1_shiftrot,outputlist[[filteredOrderIndex[i]]]$g1_shiftrot)
#     g2_shiftrot<-rbind(g2_shiftrot,outputlist[[filteredOrderIndex[i]]]$g2_shiftrot)
#     bestRank<-rbind(bestRank,outputlist[[filteredOrderIndex[i]]]$rank)    
#     
#   }
#   
#   
#   output<-data.frame(graph1=g1_shiftrot, graph2=g2_shiftrot, bestEdgePair=bestAlignEdge, bestScore=bestMatchScore, bestRank=bestRank) 
#   
#   
#   
#   if(f_plot==1){
#     for(i in 1:length(minIndices)){
#       windows()
#       plotreggraphs(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]])
#       #plotreggraphs_t(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]]) #for SNIR
#       titlestring<-paste("Alignment on rank", i," best edge pair (# ",bestIndex[i], ")", sep="")
#       title(titlestring)
#       
#     }#end of i loop
#     
#   }
#   
#   output
#   
# }# end of function
# ########################################################################################################
# 
# 
# registerGraphs_v3<-function(g1,g2,tol,maxtheta, f_plot){
# 
#   g1_edgelength<-get.edge.value(g1,"edgelength")
#   g2_edgelength<-get.edge.value(g2,"edgelength")
#   g1_edgeslope<-get.edge.value(g1,"edgeslope")
#   g2_edgeslope<-get.edge.value(g2,"edgeslope")
# 
#   g1_edgelist<-as.matrix.network(g1, matrix.type="edgelist")
#   g2_edgelist<-as.matrix.network(g2, matrix.type="edgelist")
# 
# 
#   #match edges by comparing lengths and slopes
#   #===========================================
# 
#   edge_matching <- matrix(0, length(g1_edgelength), length(g2_edgelength) )
# 
#   for (i in 1:length(edge_matching[,1])){
#     for (j in 1:length(edge_matching[1,])){
# 
# 
#       #edge_matching[i,j] <- abs(g1_edgelength[i]-g2_edgelength[j])
# 
#       #edge_matching[i,j] <- sqrt((g1_edgelength[i]-g2_edgelength[j])**2+(g1_edgeslope[i]-g2_edgeslope[j])**2)
#        edge_matching[i,j] <- abs( (g1_edgelength[i]-g2_edgelength[j]) ) + abs( tan( (g1_edgeslope[i]-g2_edgeslope[j])*pi/2  ) )
# 
# 
#     }
#   }
# 
#   ranking <- order(edge_matching, decreasing=FALSE)
# 
#   #select all edge pairs with average length greater than 10. Put them in 1 structure
#   #select all other edge pairs and put them in another structure
#   edgePairs_long<-numeric()
#   edgePairs_short<-numeric()
#   longTh<-mean(c(median(g1_edgelength),median(g2_edgelength) ) )
# 
#   for(i in 1:length(ranking)){
#     rank<-i
# 
#     edge2 <- ceiling(ranking[rank]/length(edge_matching[,1]))
#     edge1 <- ranking[rank]%%length(edge_matching[,1])
#     if (edge1==0) edge1 <- length(edge_matching[,1])
#     tmp<-list(e1=edge1, e2=edge2, rank=rank)
# 
#     if( g1_edgelength[edge1]>= longTh && g2_edgelength[edge2] >= longTh ) edgePairs_long<-rbind(edgePairs_long, tmp)
#     if( g1_edgelength[edge1] < longTh || g2_edgelength[edge2]< longTh ) edgePairs_short<-rbind(edgePairs_short, tmp)
#   }
#   edgePairs_long<-as.data.frame(edgePairs_long)
#   edgePairs_short<-as.data.frame(edgePairs_short)
# 
#   NEdgePairs_long<-min(length(edgePairs_long$e1), round(0.01*length(ranking)) )
#   NEdgePairs_short<-min(length(edgePairs_short$e1), round(0.01*length(ranking)) )
# 
#   bestMatchScore<-numeric()
#   bestAlignEdge<-numeric()
#   #   g1_shiftrot<-g1
#   #   g2_shiftrot<-g2
#   g1_shiftrot<-numeric()
#   g2_shiftrot<-numeric()
#   scores<-numeric()
#   alignedEdge<-numeric()
#   bestRank<-numeric()
#   bestIndex<-numeric()
# 
# 
# 
#   #registerDoMC(8)
#   outputlist_long<-foreach(i=1:NEdgePairs_long, .combine="cbind", .packages=c("network","sna", "clue")  )%dopar%{
#     #      # outputlist<-foreach(i=1:NEdgePairs, .combine="cbind", .packages=c("network","sna", "clue")  )%do%{
#     source("basicGraphFunctions_vpac.R")
# 
#     # for(i in 1:length(edgePairs_long$e1)){
# 
#     row<-edgePairs_long$e1[[i]]
#     column<-edgePairs_long$e2[[i]]
#     rank<-edgePairs_long$rank[[i]]
#     edgepair<-c(row,column)
# 
#     g1_temp<-shiftrotate(g1,row)
#     g2_temp<-shiftrotate(g2,column)
#     tempscore<-rapidscore(g1_temp$g_aligned,g2_temp$g_aligned,tol)
# 
#     #     windows()
#     #     plotreggraphs(g1_temp$g_aligned,g2_temp$g_aligned)
#     #     title(paste("Long edge Rank#",rank," with rapidscore =", round(tempscore,2)))
#     #
#     op<-list(score=tempscore, alignedEdgePair=edgepair, rank=rank, g1_shiftrot=g1_temp, g2_shiftrot=g2_temp)
#     op
# 
#   } # end of foreach loop
# 
#   outputlist_long<-as.data.frame(outputlist_long)
# 
# 
# 
#   outputlist_short<-foreach(i=1:NEdgePairs_short, .combine="cbind", .packages=c("network","sna", "clue")  )%dopar%{
#     source("basicGraphFunctions_vpac.R")
# 
#     row<-edgePairs_short$e1[[i]]
#     column<-edgePairs_short$e2[[i]]
#     rank<-edgePairs_short$rank[[i]]
#     edgepair<-c(row, column)
# 
#     g1_temp<-shiftrotate(g1,row)
#     g2_temp<-shiftrotate(g2,column)
#     tempscore<-rapidscore(g1_temp$g_aligned,g2_temp$g_aligned,tol)
#     #         windows()
#     #         plotreggraphs(g1_temp$g_aligned,g2_temp$g_aligned)
#     #         title(paste("Short edge Rank#",rank," with rapidscore =", round(tempscore,2)))
# 
#     op<-list(score=tempscore, alignedEdgePair=edgepair, rank=rank, g1_shiftrot=g1_temp, g2_shiftrot=g2_temp)
#     op
# 
#   } # end of foreach loop
# 
#   outputlist_short<-as.data.frame(outputlist_short)
# 
#   outputlist<-cbind(outputlist_long,outputlist_short)
# 
#   #find lowest score here by converting output to dataframe, geting lowest score and corresponding features
#   scorelist<-numeric()
# 
#   for(i in 1:dim(outputlist)[2]){
# 
#     scorelist<-c(scorelist,outputlist[[i]]$score)
# 
#   }
# 
#   theta1_list<-numeric()
#   for(i in 1:dim(outputlist)[2]){
#     op<-outputlist[[i]]$g1_shiftrot$theta * (180)/pi
#     theta1_list<-c(theta1_list,op)
#   }
# 
# 
#   theta2_list<-numeric()
#   for(i in 1:dim(outputlist)[2]){
#     op<-outputlist[[i]]$g2_shiftrot$theta * (180)/pi
#     theta2_list<-c(theta2_list,op)
#   }
# 
#   orderIndex<-order(scorelist, decreasing=FALSE)
#   thetaDiff<-theta1_list-theta2_list
#   thetaSatisfy<-which(abs(thetaDiff) < maxtheta )
# 
#   orderIndexThetaSatisfy<-which(orderIndex %in% thetaSatisfy) #the ordered indices of orderIndex, that satisfy theta contraint
#   filteredOrderIndex<-orderIndex[orderIndexThetaSatisfy] #index of edgepairs that satisfy theta contsraint, ordered in increasing order of score.
#   lowestScore<-outputlist[[filteredOrderIndex[1]]]$score
# 
#   minIndices<-which(scorelist[filteredOrderIndex]==lowestScore) #if there are more than 1 best match score, get all of the associated edges, graphs and ranks
# 
#   for(i in minIndices){
#     bestIndex<- rbind(bestIndex,filteredOrderIndex[i])
#     bestMatchScore<-rbind(bestMatchScore,outputlist[[filteredOrderIndex[i]]]$score)
#     bestAlignEdge<-rbind(bestAlignEdge,outputlist[[filteredOrderIndex[i]]]$alignedEdgePair)
#     g1_shiftrot<-rbind(g1_shiftrot,outputlist[[filteredOrderIndex[i]]]$g1_shiftrot)
#     g2_shiftrot<-rbind(g2_shiftrot,outputlist[[filteredOrderIndex[i]]]$g2_shiftrot)
#     bestRank<-rbind(bestRank,outputlist[[filteredOrderIndex[i]]]$rank)
# 
#   }
# 
# 
#   output<-data.frame(graph1=g1_shiftrot, graph2=g2_shiftrot, bestEdgePair=bestAlignEdge, bestScore=bestMatchScore, bestRank=bestRank)
# 
# 
# 
#   if(f_plot==1){
#     for(i in 1:length(minIndices)){
#       windows()
#       plotreggraphs(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]])
#       #plotreggraphs_t(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]]) #for SNIR
#       titlestring<-paste("Alignment on rank", i," best edge pair (# ",bestIndex[i], ")", sep="")
#       title(titlestring)
# 
#     }#end of i loop
# 
#   }
# 
#   output
# 
# }# end of function
# ########################################################################################################
# 
registerGraphs_star<-function(g1,g2,tol,maxtheta, f_plot, L){



  g1_edgelength<-get.edge.value(g1,"edgelength")
  g2_edgelength<-get.edge.value(g2,"edgelength")
  g1_edgeslope<-get.edge.value(g1,"edgeslope")
  g2_edgeslope<-get.edge.value(g2,"edgeslope")
  g1_coords<-cbind(get.vertex.attribute(g1, "xcoord"), get.vertex.attribute(g1, "ycoord"))
  g2_coords<-cbind(get.vertex.attribute(g2, "xcoord"), get.vertex.attribute(g2, "ycoord"))

  g1_edgelist<-as.matrix.network(g1, matrix.type="edgelist")
  g2_edgelist<-as.matrix.network(g2, matrix.type="edgelist")

  g1_star_df<-get_star(g1)
  g2_star_df<-get_star(g2)

  #finding the best aligned pair of stars
  starmatching<-matrix(0, nrow=dim(g1_star_df)[1], ncol=dim(g2_star_df)[1])
  for(k1 in 1:dim(g1_star_df)[1]){
    star1<-g1_star_df[k1,]
    lengths1<-as.numeric( unlist(star1$star))
    angles1<-as.numeric( unlist( star1$angles))

      for(k2 in 1:dim(g2_star_df)[1]){
        star2<-g2_star_df[k2,]
        lengths2<-as.numeric( unlist(star2$star))
        angles2<-as.numeric( unlist( star2$angles))

        d_l<-sqrt( (lengths1[1]-lengths2[1])^2+(lengths1[2]-lengths2[2])^2+(lengths1[3]-lengths2[3])^2)
        a1<-min(abs(angles1[1]-angles2[1]), 360-abs(angles1[1]-angles2[1]))*pi/180
        a2<-min(abs(angles1[2]-angles2[2]), 360-abs(angles1[2]-angles2[2]))*pi/180
        #d_a<-sqrt( (mean(c(lengths1[2:3], lengths2[2:3]) )*a1)^2 + (mean(c(lengths1[2:3], lengths2[2:3]) )*a2 )^2 )
        d_a<-sqrt(a1^2+a2^2)

        #d<-d_a+d_l

        d<-sqrt( sum((lengths1-lengths2)^2)+a1^2+a2^2)
        starmatching[k1,k2]<-d

    }#end of k2 loop
  }#end of k1 loop

  ranking <- order(starmatching, decreasing=FALSE)
  #r<-round(p*length(ranking))
 r<-min(L, length(starmatching) )

  bestMatchScore<-numeric()
  bestAlignEdge<-numeric()
  #   g1_shiftrot<-g1
  #   g2_shiftrot<-g2
  g1_shiftrot<-numeric()
  g2_shiftrot<-numeric()
  scores<-numeric()
  alignedEdge<-numeric()
  bestRank<-numeric()
  bestIndex<-numeric()

  outputlist<-foreach(i=1:r, .combine="cbind", .packages=c("network","sna", "clue")  )%dopar%{
  #outputlist<-foreach(i=1:r, .combine="cbind", .packages=c("network","sna", "clue")  )%do%{
    source("basicGraphFunctions_vpac.R")

    rank<-i

    column<-ceiling(ranking[rank]/length(starmatching[,1]))
    row<-ranking[rank]%%length(starmatching[,1])
    if (row==0) row <- length(starmatching[,1])

    starpair<-c(row,column)

    star1<-g1_star_df[row,]
    star2<-g2_star_df[column,]
    #find the longest edge of the stars
    edge1<-c(unlist(star1$center), unlist(star1$neighb)[1])
    edge2<-c(unlist(star2$center), unlist(star2$neighb)[1])

    #get the indices of these edges in the edgelists
    edge1Index<-union(intersect( which(g1_edgelist[,1]==edge1[1]), which(g1_edgelist[,2]==edge1[2])), intersect( which(g1_edgelist[,1]==edge1[2]), which(g1_edgelist[,2]==edge1[1])) )[1]
    edge2Index<-union(intersect( which(g2_edgelist[,1]==edge2[1]), which(g2_edgelist[,2]==edge2[2])), intersect( which(g2_edgelist[,1]==edge2[2]), which(g2_edgelist[,2]==edge2[1])) )[1]

    # shift and rotate
    g1_temp<-shiftrotate(g1, edge1Index)
    g2_temp<-shiftrotate(g2,edge2Index)
    tempscore<-rapidscore(g1_temp$g_aligned,g2_temp$g_aligned,tol)

#         windows()
#         plotreggraphs(g1_temp$g_aligned,g2_temp$g_aligned)
#         title(paste("Long edge Rank#",rank," with rapidscore =", round(tempscore,2)))
#
    op<-list(score=tempscore, alignedEdgePair=starpair, rank=rank, g1_shiftrot=g1_temp, g2_shiftrot=g2_temp)
    op

  } # end of foreach loop
  outputlist<-as.data.frame(outputlist)

  #do same as other reg functions
  scorelist<-numeric()
  for(i in 1:dim(outputlist)[2]){

    scorelist<-c(scorelist,outputlist[[i]]$score)

  }

  theta1_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g1_shiftrot$theta * (180)/pi
    theta1_list<-c(theta1_list,op)
  }


  theta2_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g2_shiftrot$theta * (180)/pi
    theta2_list<-c(theta2_list,op)
  }

  orderIndex<-order(scorelist, decreasing=FALSE)
  thetaDiff<-theta1_list-theta2_list
  thetaSatisfy<-which(abs(thetaDiff) < maxtheta )

  orderIndexThetaSatisfy<-which(orderIndex %in% thetaSatisfy) #the ordered indices of orderIndex, that satisfy theta contraint
  filteredOrderIndex<-orderIndex[orderIndexThetaSatisfy] #index of edgepairs that satisfy theta contsraint, ordered in increasing order of score.
  lowestScore<-outputlist[[filteredOrderIndex[1]]]$score

  minIndices<-which(scorelist[filteredOrderIndex]==lowestScore) #if there are more than 1 best match score, get all of the associated edges, graphs and ranks

  for(i in minIndices){
    bestIndex<- rbind(bestIndex,filteredOrderIndex[i])
    bestMatchScore<-rbind(bestMatchScore,outputlist[[filteredOrderIndex[i]]]$score)
    bestAlignEdge<-rbind(bestAlignEdge,outputlist[[filteredOrderIndex[i]]]$alignedEdgePair)
    g1_shiftrot<-rbind(g1_shiftrot,outputlist[[filteredOrderIndex[i]]]$g1_shiftrot)
    g2_shiftrot<-rbind(g2_shiftrot,outputlist[[filteredOrderIndex[i]]]$g2_shiftrot)
    bestRank<-rbind(bestRank,outputlist[[filteredOrderIndex[i]]]$rank)

  }


  output<-data.frame(graph1=g1_shiftrot, graph2=g2_shiftrot, bestEdgePair=bestAlignEdge, bestScore=bestMatchScore, bestRank=bestRank, g1StarSize=dim(g1_star_df)[1], g2StarSize=dim(g2_star_df)[1])



  if(f_plot==1){
    windows()
    plotspatialgraph(g1)
    title("g1")
    windows()
    plotspatialgraph(g2)
    title("g2")

    for(i in 1:length(minIndices)){
      windows()
      plotreggraphs(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]])
      #plotreggraphs_t(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]]) #for SNIR
      titlestring<-paste("Alignment on rank", i," best star pair (# ",bestIndex[i], ") d=", round(bestMatchScore[i],2) , sep="")
      title(titlestring)

    }#end of i loop

  }

  output

}#end of function
# #################################################################################################################################################

registerGraphs_twostars<-function(g1,g2,tol,maxtheta, f_plot, L){

  g1_edgelength<-get.edge.value(g1,"edgelength")
  g2_edgelength<-get.edge.value(g2,"edgelength")
  g1_edgeslope<-get.edge.value(g1,"edgeslope")
  g2_edgeslope<-get.edge.value(g2,"edgeslope")
  g1_coords<-cbind(get.vertex.attribute(g1, "xcoord"), get.vertex.attribute(g1, "ycoord"))
  g2_coords<-cbind(get.vertex.attribute(g2, "xcoord"), get.vertex.attribute(g2, "ycoord"))

  g1_edgelist<-as.matrix.network(g1, matrix.type="edgelist")
  g2_edgelist<-as.matrix.network(g2, matrix.type="edgelist")

  g1_twostar_df<-get_twostar(g1)
  g2_twostar_df<-get_twostar(g2)


  #finding the best aligned pair of stars
  twostarmatching<-matrix(0, nrow=dim(g1_twostar_df)[1], ncol=dim(g2_twostar_df)[1])
  for(k1 in 1:dim(g1_twostar_df)[1]){
    g1_star1<-g1_twostar_df$star1[[k1]]
    g1_star2<-g1_twostar_df$star2[[k1]]

    g1_lengths1<-as.numeric( unlist(g1_star1$star))
    g1_angles1<-as.numeric( unlist( g1_star1$angles))

    g1_lengths2<-as.numeric( unlist(g1_star2$star))
    g1_angles2<-as.numeric( unlist( g1_star2$angles))

    g1_twostar_edgeLength<-as.numeric( g1_edgelength[g1_twostar_df$edgeIndex[[k1]]] )

    for(k2 in 1:dim(g2_twostar_df)[1]){
      g2_star1<-g2_twostar_df$star1[[k2]]
      g2_star2<-g2_twostar_df$star2[[k2]]

      g2_lengths1<-as.numeric( unlist(g2_star1$star))
      g2_angles1<-as.numeric( unlist( g2_star1$angles))

      g2_lengths2<-as.numeric( unlist(g2_star2$star))
      g2_angles2<-as.numeric( unlist( g2_star2$angles))

      g2_twostar_edgeLength<-as.numeric( g2_edgelength[g2_twostar_df$edgeIndex[[k2]]] )

      l1_diff<-g1_lengths1-g2_lengths1
      l2_diff<-g1_lengths2-g2_lengths2
      a1_diff<-g1_angles1-g2_angles1
      a2_diff<-g1_angles2-g2_angles2
      edge_diff<-g1_twostar_edgeLength-g2_twostar_edgeLength

      d<-sqrt(sum(l1_diff^2)) + sqrt(sum(l2_diff^2)) + sqrt(sum(a1_diff^2)) + sqrt(sum(a2_diff^2)) + sqrt(edge_diff^2)
      d<-d/5
      twostarmatching[k1,k2]<-d

    }#end of k2 loop
  }#end of k1 loop



  #####
  ranking <- order(twostarmatching, decreasing=FALSE)
  #r<-round(p*length(ranking))
  r<-min(L, length(twostarmatching) )

  bestMatchScore<-numeric()
  bestAlignEdge<-numeric()
  #   g1_shiftrot<-g1
  #   g2_shiftrot<-g2
  g1_shiftrot<-numeric()
  g2_shiftrot<-numeric()
  scores<-numeric()
  alignedEdge<-numeric()
  bestRank<-numeric()
  bestIndex<-numeric()

  outputlist<-foreach(i=1:r, .combine="cbind", .packages=c("network","sna", "clue")  )%dopar%{
    #outputlist<-foreach(i=1:r, .combine="cbind", .packages=c("network","sna", "clue")  )%do%{
    source("basicGraphFunctions_vpac.R")

    rank<-i

    column<-ceiling(ranking[rank]/length(twostarmatching[,1]))
    row<-ranking[rank]%%length(twostarmatching[,1])
    if (row==0) row <- length(twostarmatching[,1])

    twostarpair<-c(row,column)

    twostar1<-g1_twostar_df[row,]
    twostar2<-g2_twostar_df[column,]
    #find the joining edge of the stars
    edge1Index<-as.numeric(twostar1$edgeIndex)
    edge2Index<-as.numeric(twostar2$edgeIndex)

    # shift and rotate
    g1_temp<-shiftrotate(g1, edge1Index)
    g2_temp<-shiftrotate(g2,edge2Index)
    tempscore<-rapidscore(g1_temp$g_aligned,g2_temp$g_aligned,tol)

    #         windows()
    #         plotreggraphs(g1_temp$g_aligned,g2_temp$g_aligned)
    #         title(paste("Long edge Rank#",rank," with rapidscore =", round(tempscore,2)))
    #
    op<-list(score=tempscore, alignedEdgePair=twostarpair, rank=rank, g1_shiftrot=g1_temp, g2_shiftrot=g2_temp)
    op

  } # end of foreach loop
  outputlist<-as.data.frame(outputlist)

  #do same as other reg functions
  scorelist<-numeric()
  for(i in 1:dim(outputlist)[2]){

    scorelist<-c(scorelist,outputlist[[i]]$score)

  }

  theta1_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g1_shiftrot$theta * (180)/pi
    theta1_list<-c(theta1_list,op)
  }


  theta2_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g2_shiftrot$theta * (180)/pi
    theta2_list<-c(theta2_list,op)
  }

  orderIndex<-order(scorelist, decreasing=FALSE)
  thetaDiff<-theta1_list-theta2_list
  thetaSatisfy<-which(abs(thetaDiff) < maxtheta )

  orderIndexThetaSatisfy<-which(orderIndex %in% thetaSatisfy) #the ordered indices of orderIndex, that satisfy theta contraint
  filteredOrderIndex<-orderIndex[orderIndexThetaSatisfy] #index of edgepairs that satisfy theta contsraint, ordered in increasing order of score.
  lowestScore<-outputlist[[filteredOrderIndex[1]]]$score

  minIndices<-which(scorelist[filteredOrderIndex]==lowestScore) #if there are more than 1 best match score, get all of the associated edges, graphs and ranks

  for(i in minIndices){
    bestIndex<- rbind(bestIndex,filteredOrderIndex[i])
    bestMatchScore<-rbind(bestMatchScore,outputlist[[filteredOrderIndex[i]]]$score)
    bestAlignEdge<-rbind(bestAlignEdge,outputlist[[filteredOrderIndex[i]]]$alignedEdgePair)
    g1_shiftrot<-rbind(g1_shiftrot,outputlist[[filteredOrderIndex[i]]]$g1_shiftrot)
    g2_shiftrot<-rbind(g2_shiftrot,outputlist[[filteredOrderIndex[i]]]$g2_shiftrot)
    bestRank<-rbind(bestRank,outputlist[[filteredOrderIndex[i]]]$rank)

  }


  output<-data.frame(graph1=g1_shiftrot, graph2=g2_shiftrot, bestEdgePair=bestAlignEdge, bestScore=bestMatchScore, bestRank=bestRank, g1StarSize=dim(g1_twostar_df)[1], g2StarSize=dim(g2_twostar_df)[1])



  if(f_plot==1){
    windows()
    plotspatialgraph(g1)
    title("g1")
    windows()
    plotspatialgraph(g2)
    title("g2")

    for(i in 1:length(minIndices)){
      windows()
      plotreggraphs(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]])
      #plotreggraphs_t(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]]) #for SNIR
      titlestring<-paste("Alignment on rank", i," best twostar pair (# ",bestIndex[i], ") d=", round(bestMatchScore[i],2) , sep="")
      title(titlestring)

    }#end of i loop

  }

  output

  #######

}#end of function



#returns a star data structure
get_star<-function(g){

  g_edgelength<-get.edge.value(g,"edgelength")
  g_edgeslope<-get.edge.value(g,"edgeslope")
  g_coords<-cbind(get.vertex.attribute(g, "xcoord"), get.vertex.attribute(g, "ycoord"))
  g_edgelist<-as.matrix.network(g, matrix.type="edgelist")

  g_star<-numeric()

  for(y in 1:network.size(g)){
    neighbors<-get.neighborhood(g,y, type="combined")
    if(length(neighbors)!=3) next

    starEdges<-numeric()
    for(z in 1:length(neighbors)){
      edge<-numeric()
      if(y<neighbors[z]) edge<-c(y,neighbors[z])
      if(y>neighbors[z]) edge<- c(neighbors[z],y)
      edgeIndex<-union( intersect(which(g_edgelist[,1]==edge[1]),which(g_edgelist[,2]==edge[2]) ), intersect(which(g_edgelist[,1]==edge[2]),which(g_edgelist[,2]==edge[1]) ) )[1]
      length<-g_edgelength[edgeIndex]
      starEdges<-c(starEdges, length )
    }#end of z loop

    #order edge in increasing order of length
    orderInd<-order(starEdges, decreasing=TRUE)
    starEdges<-starEdges[orderInd]
    neighbOrd<-neighbors[orderInd]

    #find internal angles
    thetas<-numeric()
    for(z in 1:length(starEdges)){
      if(z==3) next

      l1<-starEdges[z]
      v1<-neighbOrd[z]

      l2<-starEdges[z+1]
      v2<-neighbOrd[z+1]

      l3<-as.numeric(dist(g_coords[c(v1,v2),], method="euclidean"))
      tmp<-round((l1^2+l2^2-l3^2)/(2*l1*l2),2)

      theta<-min(acos(tmp)*180/pi, 360-acos(tmp)*180/pi)
      thetas<-c(thetas, theta)
    }#end of z loop
    tmp<-list(center=y, neighb=neighbOrd, star=starEdges, angles=thetas)
    g_star<-rbind(g_star,tmp)
  }#end of y loop
  g_star_df<-as.data.frame( g_star)


  g_star_df


}# end of get_star function

get_twostar<-function(g){

  g_edgelist<-as.matrix.network(g, matrix.type="edgelist")
  g_adjacency<-as.matrix.network(g, matrix.type="adjacency")

  g_star_df<-get_star(g)
  #find which of the star centers are connected
  centerlist<-as.numeric(unlist(g_star_df$center))
  g_twostar<-numeric()
  for(i in centerlist){
    for(j in centerlist){
      if(i==j) next
      if(i>j) next
      if(g_adjacency[i,j]==1){
        index_i<-which(centerlist==i)
        index_j<-which(centerlist==j)
        index_edge<-intersect(which(g_edgelist[,1]==i),which(g_edgelist[,2]==j ) )
        tmp<-list(star1=g_star_df[index_i,], star2=g_star_df[index_j,], edgeIndex=index_edge)
        g_twostar<-rbind(g_twostar, tmp)
      }  #end of if
    }# end of j loop
  }# end of i loop

  g_twostar_df<-as.data.frame(g_twostar)

  g_twostar_df

}#end of function
#




#################################################################################################################################################
# require(pdist)
# registerGraphs_cycle<-function(g1,g2,tol,maxtheta, f_plot){
#   coords_1<-cbind(get.vertex.attribute(g1, "xcoord"), get.vertex.attribute(g1, "ycoord") )
#   coords_2<-cbind(get.vertex.attribute(g2, "xcoord"), get.vertex.attribute(g2, "ycoord") )
#
#   MaxLength<-10
#   cycList1<-get_cycle(g1, MaxLength)
#   cycList2<-get_cycle(g2, MaxLength)
#   l1<-as.numeric(cycList1$l)
#   l2<-as.numeric(cycList2$l)
#   d_l1l2<-as.matrix( pdist(as.matrix(l1),as.matrix(l2) ) )
#   matchpair<-which(d_l1l2==0, arr.ind=TRUE)
#
#   matchedCycles<-numeric()
#   for(i in 1:dim(matchpair)[1]){
#     g1_ind<-matchpair[i,1]
#     g2_ind<-matchpair[i,2]
#     g1_cyc<-cycList1$gCyc[[g1_ind]]
#     g2_cyc<-cycList2$gCyc[[g2_ind]]
#     windows()
#     par(mfrow=c(2,1))
#     plotspatialgraph(g1_cyc)
#     plotspatialgraph(g2_cyc)
#
#     path1<-cycList1$path[[g1_ind]]
#     path2<-cycList2$path[[g2_ind]]
#
#     #create the cycle's traversal path as an ordered list of vertices
#     p1<-get_path(g1_cyc,1)
#     path1<-path1[p1]
#     coords1<-coords_1[path1,]
#     p2<-get_path(g2_cyc,1)
#     path2<-path2[p2]
#     coords2<-coords_2[p2,]
#
#     #get the included angles
#     angles1<-rep(0, times=length(path1))
#     angles2<-rep(0, times=length(path2))
#     for(j in 1:length(path1)){
#       #vertexids
#       ind1<-j
#       ind2<-(j+1)%%length(path1)
#       if(ind2==0)ind2<-length(path1)
#       ind3<-(j+2)%%length(path1)
#       if(ind3==0) ind3<-length(path1)
#       A<-path1[ind1]
#       B<-path1[ind2]
#       C<-path1[ind3]
#       #lengths
#       AB<-as.numeric(dist(rbind(coords1[ind1,], coords1[ind2,]), method="euclidean" ) )
#       BC<-as.numeric(dist(rbind(coords1[ind2,], coords1[ind3,]), method="euclidean" ) )
#       CA<-as.numeric(dist(rbind(coords1[ind3,], coords1[ind1,]), method="euclidean" ) )
#
#       angles1[ind2]<-get_angle(AB, BC, CA)
#     }
#
#     for(j in 1:length(path2)){
#       #vertexids
#       ind1<-j
#       ind2<-(j+1)%%length(path2)
#       if(ind2==0)ind2<-length(path2)
#       ind3<-(j+2)%%length(path2)
#       if(ind3==0) ind3<-length(path2)
#       A<-path2[ind1]
#       B<-path2[ind2]
#       C<-path2[ind3]
#       #lengths
#       AB<-as.numeric(dist(rbind(coords2[ind1,], coords2[ind2,]), method="euclidean" ) )
#       BC<-as.numeric(dist(rbind(coords2[ind2,], coords2[ind3,]), method="euclidean" ) )
#       CA<-as.numeric(dist(rbind(coords2[ind3,], coords2[ind1,]), method="euclidean" ) )
#
#       angles2[ind2]<-get_angle(AB, BC, CA)
#     }#end of j loop
#
#     #find a distance measure with length(angles2) rotations of the angles2 vector
#     d_rot<-rep(0, times=length(angles2))
#     tmp<-angles2
#     for(j in 1:length(angles2)){
#
#       a<-angles1-tmp
#       d<-sqrt(sum(a*a))/length(angles2)
#       d_rot[j]<-d
#
#       tmp<-tmp[rot(length(angles2))]
#
#     }#end of j loop
#     minIndex<-which(d_rot==min(d_rot))
#     #shift minIndex-1 times, if minIndex=1, dont rotate.
#
#     if(minIndex>1){
#       for(j in 1:(minIndex-1)){
#         path2<-path2[rot(length(path2))]
#         coords2[,1]<-coords2[,1][rot(length(path2))]
#         coords2[,2]<-coords2[,2][rot(length(path2))]
#         angles2<-angles2[rot(length(path2))]
#       }#end of j loop
#     }
#    #Align on all possible edge pairs going along the cycle
#
#     tmp<-list(d=min(d_rot), path1=path1, path2=path2)
#     matchedCycles<-rbind(matchedCycles, tmp)
#
#
#   }#end of i loop
#   matchedCycles<-as.data.frame(matchedCycles)
#
#   ordInd<-order(as.numeric(matchedCycles$d))
#   ordMatchedCycles<-matchedCycles[ordInd,]
#
#
#
#
# }#end of function
#
#
# ##############################################################################################################################################
# #this function is called by registerGraphs_cycle
#
# get_cycle<-function(g1, MaxLength){
#
#   cycList<-numeric() #cycles for this graph
#
#   cyc<-kcycle.census(g1, maxlen = MaxLength, mode = "graph", tabulate.by.vertex = TRUE, cycle.comembership = "bylength")
#   # a MaxLength-1 x |V| matrix, first column gives number of cycles of rowIndex+1 length, column gives numbr of times the cycle of 'row+1' length pass through a vertex
#   cycCountMx<-cyc$cycle.count
#   cycMxRowInd<-sort(as.numeric(which(cycCountMx[,1]>0)), decreasing=FALSE) #get all cycles of length greater than 2
#
#   for(j in 1:length(cycMxRowInd)){
#     length<-cycMxRowInd[j]+1
#     Row<-cycMxRowInd[j]
#
#     #1. there is only one such cycle of this 'length'
#     if(cycCountMx[Row,1]==1){
#       path<-as.numeric( which(cycCountMx[Row,-1]==1))
#       g_path<-g1
#       indToDel<-setdiff(1:network.size(g1), path)
#       g_path<-delete.vertices(g_path, indToDel)
#
#
#
#       #         windows()
#       #         plotspatialgraph(g_path)
#       tmp<-list(gCyc=g_path, path=path, l=length, person=i)
#       cycList<-rbind(cycList,tmp)
#
#
#
#     }#end of if
#
#     #2. there are more than one cycles of this length in the graph
#     if(cycCountMx[Row,1]>1){
#       ncycles<-cycCountMx[Row,1]
#       oneTimePassers<-as.numeric( which(cycCountMx[Row,-1]==1))
#       if(length(oneTimePassers)==0) next
#       for(k in 1:length(oneTimePassers)){
#         v<-oneTimePassers[k]
#         path<-as.numeric(which(cyc$cycle.comemb[Row,,v]>0))
#         g_path<-g1
#         indToDel<-setdiff(1:network.size(g1), path)
#         g_path<-delete.vertices(g_path, indToDel)
#
#
#         #           windows()
#         #           plotspatialgraph(g_path)
#         tmp<-list(gCyc=g_path, path=path, l=length,person=i)
#         cycList<-rbind(cycList,tmp)
#
#       }#end of k loop
#
#     }#end of if
#   }#endf of j loop
#
#   cycList<-as.data.frame(cycList)
#   newList<-cycList
#
# #   j<-1
# #   repeat{
# #     path1<-newList$path[[j]]
# #
# #     colToDel<-numeric()
# #     for(k in (j+1): dim(newList)[1]){
# #       path2<-newList$path[[k]]
# #       overlap<-length(intersect(path1,path2))/sqrt((length(path1)*length(path2)))
# #       if(overlap>0.9) colToDel<-c(colToDel,k)
# #     }#end of k loop
# #     if(length(colToDel)>0) newList<-newList[-colToDel,]
# #     j<-j+1
# #     if(j>=dim(newList)[1]) break
# #   }#end of repeat loop
#
#   #if one cycle is a complete subset of a bigger cycle, delete the bigger cycle
#   j<-1
#   repeat{
#     path1<-newList$path[[j]]
#
#     colToDel<-numeric()
#     for(k in (j+1): dim(newList)[1]){
#       path2<-newList$path[[k]]
#       l<-length( which(is.element(path1,path2)==TRUE)) #how many of path1 is present in path2
#       overlap<-l/length(path1)
#       if(overlap>0.9) colToDel<-c(colToDel,k)
#     }#end of k loop
#     if(length(colToDel)>0) newList<-newList[-colToDel,]
#     j<-j+1
#     if(j>=dim(newList)[1]) break
#   }#end of repeat loop
#
#   #return the dataframe
#   newList
# }#end of function
#
# ##################################################################################################################################################
# #function to get included angle for sides1 and 2, given the length of three sides.
# get_angle<-function(s1,s2,s3){
#   tmp<-(s1^2+s2^2-s3^2)/(2*s1*s2)
#   theta<-min(acos(tmp)*180/pi, 360-acos(tmp)*180/pi)
#
#   round(theta,2)
# }
#
# ####################################################################################################
#
# get_path<-function(g, startnode){
#   g_temp<-g
#
#
# #   windows()
# #   plotspatialgraph(g_temp)
#   start<-startnode
#   neighb<-get.edges(g_temp, start, neighborhood="combined")[[1]]$inl
#   g_temp[start, neighb]<-0
#
#   gdist<-geodist(g_temp, predecessors=TRUE)
#   distMx<-gdist$gdist
#   nList<-gdist$predecessor
#
#   d<-degree(g_temp, gmode="graph")
#   endNodes<-which(d==1)
#   end<-endNodes[which(endNodes!=start)]
#   #list cycle in order of nodes
#
#   #get the path
#   path<-end
#   hopNode<-unlist(nList[[start]][end])[1]
#            repeat{
#              path<-c(path, hopNode)
#              hopNode<-unlist(nList[[start]][hopNode])[1]
#              if(hopNode==start) break
#            }
#   path<-c(path,start)
#
#   path
# }
# ###########################
# #vector rotation functions
# rot <- function(x) (1:x %% x) +1
# #rot(9)
# rotvec <- function(vec) vec[rot(length(vec))]