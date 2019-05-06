### Functions around getting the anchor
library(ggplot2)

getAnchor<-function(biometric, person, LR){
  anchor<-network(as.matrix(1) )
  gList<-numeric()
  
  session_range<-1:3
  number_range<-1:4
  for(si in session_range){
    for(ni in number_range){
      session<-session_range[si]
      number<-number_range[ni]
      g<-getGraph(biometric,person,session,number, LR)
      gtemp<-list(graph=g, person=person, session=session, number=number, LR=LR)
      gList<-rbind(gList,gtemp)
    }
  }
  gList<-as.data.frame(gList)
  
  # get list of samples for creating anchor 
  sample_index<-setdiff(1:12, testSampleList)
 
  parameters_anchor<-parameters
  edgeCost_list<-seq(from=4, to=16, by=2)
  for(e in 1:length(edgeCost_list)){
    parameters_anchor$cid_e<- edgeCost_list[e]
    
    g_input_1<-gList$graph[[sample_index[1]]]
    g_input_2<-NULL
    
      for(i in 2:length(sample_index) ){
      g_input_2<-gList$graph[[sample_index[i]]]
      g_common<-getCommonGraph(g_input_1, g_input_2, parameters_anchor)
      
       #print(i)
      g_input_1<-g_common$mcs

    
      # print(network.size(g_common$mcs))
      
    }# end of i loop
    tmp_anchor<-getComponent_max10(g_common$mcs)
    # two tests for the anchor
    # 1. Is it of size >5 
    size_test<-(network.size(tmp_anchor)>5)
    #2. does it have at least one claw
    claw_test<-max(degree(tmp_anchor, gmode="graph"))>2
   
     if(size_test==TRUE && claw_test==TRUE){
      anchor<-tmp_anchor
      break
     }
  }#end of e loop
  

return(anchor)
}


getCommonGraph<-function(g1, g2, parameters){
  tol<-parameters$tol
  maxtheta<-parameters$maxtheta
  f_plot<-parameters$f_plot
  L<-parameters$L
  cid_n<-parameters$cid_n
  cid_e<-parameters$cid_e
  f_edgeinclude<-parameters$f_edgeinclude
  
  reg_graphs<-registerGraphs_v1(g1,g2,tol,maxtheta, f_plot, L)
  g1_reg<-reg_graphs$graph1.g_aligned[[1]]
  g2_reg<-reg_graphs$graph2.g_aligned[[1]]
  path_edge<-getGraphEditPath_edgeBased(g1_reg,g2_reg,cid_n, cid_e)
  g_mcs<-getMCS_edgeBased(g1_reg, g2_reg, path_edge$editPath) # mcs between g1 and g2
  if(network.size(g_mcs)<2)
    g_mcs<-network(as.matrix(1) )
  
  list(mcs=g_mcs, g1_reg=g1_reg, g2_reg=g2_reg)
}

findAnchor<-function(g, g_anchor, parameters){
  tol<-parameters$tol
  maxtheta<-parameters$maxtheta
  f_plot<-1
  L<-parameters$L
  cid_n<-parameters$cid_n
  cid_e<-parameters$cid_e
  f_edgeinclude<-parameters$f_edgeinclude
  
  tmp<-registerAnchor(g, g_anchor,tol,maxtheta, f_plot, L)
  return(tmp)
}

# g_anchor<-anchor
# g<-gList$graph[[3]]
# findAnchor(g,g_anchor,  parameters)

getLargestComponent<-function(g){
  lgc<-component.largest(g,connected="weak") # list of node associations
  largest_component<-g
  v_id<-which(lgc==FALSE)
  delete.vertices(largest_component,v_id)
  # quartz()
  # plotspatialgraph(largest_component)
  if(network.size(largest_component)<2)
    largest_component<-network.size(1)
  return(largest_component)
}

# Function to return a connected component of the graph whose size is no more than 10
getComponent_max10<-function(g){
  component<-g
  
  cd<-component.dist(g,connected="weak")
  size_sort<-order(cd$csize, decreasing =TRUE)
  for(i in 1:length(size_sort)){
    if(cd$csize[size_sort[i]]>10) next
   comp<-which(cd$csize==cd$csize[size_sort[i]])[1]
   v_id<-which(cd$membership!=comp)
   delete.vertices(component, v_id)
   break
  }
  # quartz()
  # plotspatialgraph(component)
  if(network.size(component)<2)
    component<-network.size(1)
  return(component)
  
}

plotAnchorAlignment<-function(person, anchor, testSampleList){
  gList<-numeric()
  
  session_range<-1:3
  number_range<-1:4
  for(si in session_range){
    for(ni in number_range){
      session<-session_range[si]
      number<-number_range[ni]
      g<-getGraph(biometric,person,session,number, LR)
      gtemp<-list(graph=g, person=person, session=session, number=number, LR=LR)
      gList<-rbind(gList,gtemp)
    }
  }
  gList<-as.data.frame(gList)
  
  quartz()
  par(mfrow=c(3,2), oma=c(0,0,2,0))
  for(t in 1:length(testSampleList)){
    t_index<-testSampleList[t]
    tmp<-findAnchor(gList$graph[[t_index]], anchor, parameters)
    plotreggraphs(tmp$graph1.g_aligned[[1]],tmp$graph2.g_aligned[[1]])
    title(paste0("Anchor_overlap with sample ",t_index," = ",round(1-tmp$bestScore[[1]],4)*100 ) )
  }
}

# This function returns a list of overlaps of the anchor with the test samples. 
# Overlap is expressed in percentage.
getAnchorOverlap<-function(person, anchor, testSampleList){
  
  overlap_list<-numeric()
  
  gList<-numeric()
  
  session_range<-1:3
  number_range<-1:4
  for(si in session_range){
    for(ni in number_range){
      session<-session_range[si]
      number<-number_range[ni]
      g<-getGraph(biometric,person,session,number, LR)
      gtemp<-list(graph=g, person=person, session=session, number=number, LR=LR)
      gList<-rbind(gList,gtemp)
    }
  }
  gList<-as.data.frame(gList)
  
  for(t in 1:length(testSampleList)){
    t_index<-testSampleList[t]
    tmp<-findAnchor(gList$graph[[t_index]], anchor, parameters)
    overlap_list<-c(overlap_list, round(1-tmp$bestScore[[1]], 4)*100  )
    }
  return(overlap_list)
}

# function to analyse the results of finding an anchor
getAnchorAnalysis<-function(biometric, LR){
  anchorList_all<-numeric()
  overlap_mx_all<-numeric()
  
 
    outputfile<-paste0(getwd(),"/Output/AnchorList_",biometric,"_", LR,".Rdata")
    load(outputfile)
    outputfile<-paste0(getwd(),"/Output/OverlapMatrix_",biometric,"_", LR,".Rdata")
    load(outputfile)
    anchorList_all<-rbind(anchorList_all, anchorList)
    overlap_mx_all<-rbind(overlap_mx_all, overlap_mx)
    
    outputfile<-paste0(getwd(),"/Output/AnchorList2_",biometric,"_", LR,".Rdata")
    load(outputfile)
    outputfile<-paste0(getwd(),"/Output/OverlapMatrix2_",biometric,"_", LR,".Rdata")
    load(outputfile)
    anchorList_all<-rbind(anchorList_all, anchorList)
    overlap_mx_all<-rbind(overlap_mx_all, overlap_mx)
    
    
    outputfile<-paste0(getwd(),"/Output/AnchorList3_",biometric,"_", LR,".Rdata")
    load(outputfile)
    outputfile<-paste0(getwd(),"/Output/OverlapMatrix3_",biometric,"_", LR,".Rdata")
    load(outputfile)
    anchorList_all<-rbind(anchorList_all, anchorList)
    overlap_mx_all<-rbind(overlap_mx_all, overlap_mx)
    
    
    outputfile<-paste0(getwd(),"/Output/AnchorList4_",biometric,"_", LR,".Rdata")
    load(outputfile)
    outputfile<-paste0(getwd(),"/Output/OverlapMatrix4_",biometric,"_", LR,".Rdata")
    load(outputfile)
    anchorList_all<-rbind(anchorList_all, anchorList)
    overlap_mx_all<-rbind(overlap_mx_all, overlap_mx)
    
 
  #main questions we ask in the experiment are:
  # 1. what proportion fails to find an anchor?
  anchor_sizes<-numeric()
  for(i in 1:dim(anchorList_all)[1])
    anchor_sizes<-c(anchor_sizes, network.size(anchorList_all$anchor[[i]]) )
  failed_anchor<-length(which(anchor_sizes==1))/length(anchor_sizes)
  failed_anchor_indices<-which(anchor_sizes==1)
  
  #2. When we find an anchor what is the anchor overlap overall for the database?
  hist_overlap<-hist(c(overlap_mx_all[-failed_anchor_indices,]), breaks=seq(from=0, to=100, by=10), plot=FALSE )
  mean(overlap_mx_all[-failed_anchor_indices,])
  sqrt(var(as.numeric(overlap_mx_all[-failed_anchor_indices,]) ))
  # 3. When we did find an anchor, how many times was it found in the test samples?
  # setting threshold PV,WV=70% as the min overlap needed for successful anchor alignment in test sample
  success_counts<-rep(0, times=dim(overlap_mx_all)[1])
  for(i in 1:dim(overlap_mx_all)[1]){
    success_counts[i]<-length(which(overlap_mx_all[i,]>=70))
  }
anchorFindingLikelihood_mean<-mean(success_counts[-failed_anchor_indices])
anchorFindingLikelihood_sd<-sqrt(var(success_counts[-failed_anchor_indices]))
hist_likelihood<-hist(success_counts[-failed_anchor_indices], breaks=seq(from=0, to=6, by=1), plot = FALSE)
p1<-plotHistogram(hist_overlap, "Anchor overlap (measure of anchor reliability) ", "Anchor overlap (%)", "Frequency")
p2<-plotHistogram(hist_likelihood, "Number of correct registrations (overlap>=70%)", "Number of BGs that were successfully registered", "Frequency")
#quartz()
p1
#quartz()
p2

list(plot1=p1, plot2=p2, failureRate=failed_anchor)
}

plotHistogram<-function(hist_data, titlestring, xlab, ylab){
  plotData<-data.frame(overlap=hist_data$mids, counts=hist_data$counts)
  g<-ggplot(plotData, aes(x=overlap, y=counts))+geom_bar(stat="identity", color="black", fill="steelblue")+
  ggtitle(titlestring)+
    xlab(xlab)+
    ylab(ylab)+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"))+
    theme(plot.title = element_text(size = 20))
  g
  
}
