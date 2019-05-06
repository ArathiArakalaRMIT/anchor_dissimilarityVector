

## visualise the data

source("basicFunctions.R")
source("graphRegFunctions.R")
source("graphRegMatchFunc.R")
source("graphMatchFunc_lib.R")
source("anchorFunc.R")

biometric<-"W"
person_range<-1:50
session_range<-1:3
number_range<-1:4
LR<-"L"
anchorList<-numeric()
#parameter values for the chosen biometric
# for "P", the L depth is 220 (based on IET paper), tol=8
#parameters<-list(tol=8, maxtheta=360, f_plot=0, L=220, cid_n=1, cid_e=8, f_edgeinclude=1 )

# for "W", tol=15 (bit more noise), L=60
parameters<-list(tol=15, maxtheta=360, f_plot=0, L=60, cid_n=1, cid_e=8, f_edgeinclude=1 )

#1.
testSampleList<-c(3,4,7,8,11,12)
#2. testSampleList<-c(1,2,5,6,9,10)
#3. testSampleList<-c(1,3,5,7,9,11)
#4. testSampleList<-c(2,4,6,8,10,12)

overlap_mx<-matrix(0, nrow=length(person_range), ncol=length(testSampleList))
colnames(overlap_mx)<-testSampleList
rownames(overlap_mx)<-person_range
for(pi in 1:length(person_range) ){
  
     person<-person_range[pi]
     anchor<-getAnchor(biometric, person, LR) 
    # plotAnchorAlignment(person, anchor, testSampleList)
     
     tmp_anchor<-list(anchor=anchor,biometric=biometric, person=person, LR=LR)
     anchorList<-rbind(anchorList, tmp_anchor)
     
     if(network.size(anchor)==1)
        tmp_list<-rep(0, times=length(testSampleList))
     if(network.size(anchor)>1)
        tmp_list<-getAnchorOverlap(person, anchor, testSampleList)
       
     overlap_mx[pi,]<-tmp_list
     
 
 print(person)
 print(tmp_list)
}
# with file naming, put 2,3,4 to indicate diff selection of test samples
anchorList<-as.data.frame(anchorList)
outputfile<-paste0(getwd(),"/Output/AnchorList4_",biometric,"_", LR,".Rdata")
save(anchorList, file=outputfile)
load(outputfile)

outputfile<-paste0(getwd(),"/Output/OverlapMatrix4_",biometric,"_", LR,".Rdata")
save(overlap_mx, file=outputfile)
load(outputfile)

result<-getAnchorAnalysis(biometric, LR)
result$failureRate
quartz()
result$plot1
quartz()
result$plot2

filename<-paste0(getwd(),"/Output/hist_overlap_",biometric,"_", LR,".eps")
ggsave(filename, plot = result$plot1)

filename<-paste0(getwd(),"/Output/hist_likelihood_",biometric,"_", LR,".eps")
ggsave(filename, plot = result$plot2)

# #testing
# quartz()
# plotspatialgraph(anchorList$anchor[[1]])


# code to create the 4 types of induced subgraphs
person<-10
session<-1
number1<-1
number2<-2
g1<-getGraph(biometric,person,session,number1, LR)
g2<-getGraph(biometric,person,session,number2, LR)
reg_graphs<-registerGraphs_v1(g1,g2,parameters$tol,parameters$maxtheta, parameters$f_plot, parameters$L)

quartz()
plotspatialgraph(reg_graphs$graph1.g_aligned[[1]])
quartz()
plotspatialgraph(reg_graphs$graph2.g_aligned[[1]])

g1_reg<-reg_graphs$graph1.g_aligned[[1]]
g2_reg<-reg_graphs$graph2.g_aligned[[1]]
cid_n<-parameters$cid_n
cid_e<-parameters$cid_e
f_edgeinclude<-1
path_node<-getGraphEditPath_nodeBased(g1_reg,g2_reg,cid_n, cid_e, f_edgeinclude)
g_mcs<-getMCS_nodeBased(g1_reg, g2_reg, path_node$editPath) # mcs between g1 and g2
quartz()
plotspatialgraph(g_mcs)

path_edge<-getGraphEditPath_edgeBased(g1_reg,g2_reg,cid_n, cid_e)
g_mcs<-getMCS_edgeBased(g1_reg, g2_reg, path_edge$editPath) # mcs between g1 and g2
quartz()
plotspatialgraph(g_mcs)

path_star<-getGraphEditPath_starBased(g1_reg,g2_reg,cid_n, cid_e)
g_mcs<-getMCS_starBased(g1_reg, g2_reg, path_star$editPath) # mcs between g1 and g2
quartz()
plotspatialgraph(g_mcs)

path_twostar<-getGraphEditPath_twostarBased(g1_reg,g2_reg,cid_n, cid_e)
g_mcs<-getMCS_twostarBased(g1_reg, g2_reg, path_twostar$editPath) # mcs between g1 and g2
quartz()
plotspatialgraph(g_mcs)

# code to get the statistics of ESRID database
biometric<-"RE" 
sourceid <- "ECGHao"
wd<-"/Users/Arathi/Documents/2018/RMIT/Biometrics/Dissimilarity Vector/Outputs/RetinaGraphs/"

outputfile<-paste(wd, sourceid, "/graphList_",sourceid,".Rdata", sep="")
#save(graphList, file=outputfile)
load(outputfile)

outputfile<-paste(wd, sourceid, "/graphList_rescaled_",sourceid,".Rdata", sep="")
#save(graphList_rescaled, file=outputfile)
load(outputfile)

outputfile<-paste(wd, sourceid, "/graphList_od_",sourceid,".Rdata", sep="")
#save(graphList_od, file=outputfile)
load(outputfile)
graphList_ECG<-graphList_od

outputfile<-paste(wd, sourceid, "/graphList_rescaled_od_",sourceid,".Rdata", sep="")
#save(graphList_rescaled_od, file=outputfile)
load(outputfile)
graphList_rescaled_ECG<-graphList_rescaled_od

topology_results<-numeric()
for(i in 1:dim(graphList_rescaled_ECG)[1]){
  g<-graphList_rescaled_ECG$graph[[i]]
  tmp<-getTopology(g)
  topology_results<-rbind(topology_results, tmp)
  print(i)
}

topology_results<-as.data.frame(topology_results)

#results for Springer book chapter
#|V|
print("|V|")
mean(as.numeric(topology_results$Size), na.rm=TRUE)
sqrt(var(as.numeric(topology_results$Size), na.rm=TRUE))

#|E|
print("|E|")

mean(as.numeric(topology_results$Edgecount), na.rm=TRUE)
sqrt(var(as.numeric(topology_results$Edgecount), na.rm=TRUE))

#|E|/|V|
print("#|E|/|V|")

mean(as.numeric(topology_results$Edgenoderatio), na.rm=TRUE)
sqrt(var(as.numeric(topology_results$Edgenoderatio), na.rm=TRUE))

#|C_1|
print("|C_1|")

mean(as.numeric(topology_results$C1), na.rm=TRUE)
sqrt(var(as.numeric(topology_results$C1), na.rm=TRUE))

#|claws|
print("|claws|")

mean(as.numeric(topology_results$starSize), na.rm=TRUE)
sqrt(var(as.numeric(topology_results$starSize), na.rm=TRUE))

#|twoclaws|
print("|twoclaws|")

mean(as.numeric(topology_results$twostarSize), na.rm=TRUE)
sqrt(var(as.numeric(topology_results$twostarSize), na.rm=TRUE))
