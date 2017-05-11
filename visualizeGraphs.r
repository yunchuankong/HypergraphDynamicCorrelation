# setwd("C:\\Users\\DELL\\Dropbox\\Research_Yu\\fortress")
# setwd("C:\\Users\\kyccw\\Dropbox\\Research_Yu\\fortress")
# setwd("H:\\Research_Tianwei")

##############################################################################################
## Supervised
##############################################################################################

load('net_hum_GO.bin')
load('GO_select_human_ 0.7 0.4 1000 100 .bin')
load('label_hum_GO.bin')

load('net_yea_GO.bin')
load('GO_select_yeast_ 0.7 0.4 1000 20 .bin')
load('label_yea_GO.bin')

library(igraph)
library(GO.db)

# Match GO terms with GO IDs
# library(GO.db)
# GO.names<-mget(names(GO.select), GOTERM, ifnotfound=NA)
# findTerm=function(obj){
#     if (is.na(obj)!=T)
#         return(Term(obj))
#     else
#         return(NA)
# }
# GOdict<-as.vector(unlist(lapply(GO.names, findTerm)))

# pick up those larger than f fold change
k=dim(net)[1] # net is fully normalized 
f=1 # fold change. KEY parameter (should not be higher than ?)
length(which(net>f))
length(which(net[net>0]<=f))
hist(net[net>0],1000,freq=F,xlim=c(0,25))
abline(v=f,col='red')
net_vec_sel=as.vector(net)
net_vec_sel[net_vec_sel<f]=0
net_sel=array(net_vec_sel,c(k,k,k))
elist=which(net_sel>0,arr.ind=TRUE)
ewgt=net_sel[elist]
elist=cbind(elist,ewgt)

V=length(unique(as.vector(elist[,1:3]))) # num of modules that contain selected edges
# save(elist,file='elist_GO.bin')

# use {igraph}
# fold as weight, allow multi-edge
el=rbind(elist[,c(1,2,4)],elist[,c(1,3,4)],elist[,c(2,3,4)])
g=graph.edgelist(apply(el[,1:2],2,as.character), directed = F)
length(V(g))==V
E(g)$weight=el[,3]
# hist(unique(el[,1:2]) # degree distribution for the collapsed graph
     # # ,ylim=c(0,400)
     # )
degrees=table(el[,1:2]) # number of connections for each node
dlevels=sort(unique(degrees),decreasing=T) # degree levels
nocol=length(dlevels) # number of colors
rgb.palette <- colorRampPalette(c("red","orange","yellow" 
                                  # ,"green","cyan","blue"
                                  ), space = "rgb",bias=0.5)
col<-rgb.palette(nocol)
color=rep(NA,V)
ind=match(names(degrees),names(V(g)))
for (i in 1:nocol){
    color[ind[which(degrees==dlevels[i])]]=col[i]
}
V(g)$color=color
vlabels=names(GO.select)[as.numeric(names(V(g)))]
ewidth=(E(g)$weight-f)/(mean(E(g)$weight)-f)

# pdf(paste("full_hum_GO (f=",f,").pdf",sep=''))
# pdf(paste("full_yea_GO (f=",f,").pdf",sep=''))
# plot(g, # http://igraph.org/r/doc/plot.common.html
#      edge.width=ewidth,
#      edge.color='grey75',
#      vertex.label=vlabels,
#      vertex.label.cex=0.3,
#      vertex.label.color='black',
#      vertex.frame.color=NA,
#      vertex.size=5)
# dev.off()

#############################################################################################
## Plot 3-uniform hypergraph (specified)

pdf(paste("All_hum_GO (f=",f,").pdf",sep='')
    ,paper='USr'
    # ,width=10,height = 10
    )
for (group in 1:length(GO.select)){
    # for (group in 1:30){
    GOID = names(GO.select)[group]
    # GOID="GO:0008361"
    term = Term(GOID)
    groupID = which(names(GO.select)==GOID)
    find_keep=function(x)return(all(groupID %in% x))
    found_keep=apply(elist[,1:3],1,find_keep)
    if (length(which(found_keep==T))==0){
        cat("WARNING: the group is not involved under the current fold change threshold!\n")
        next
    }
    helist=elist[found_keep,]
    if (length(which(found_keep==T))==1){
        cat("WARNING: the group has only one edge. No hypergraph will be drawn.\n")
        cat("The edge connects:\n",Term(names(GO.select)[helist[1]]),"\n",
            Term(names(GO.select)[helist[2]]),"\n",
            Term(names(GO.select)[helist[3]]),"\n")
        next
    # } 
    # # If want a limited number of connections per plot: (15)
    # else if (nrow(helist)>15){
    #     ord = order(helist[,4],decreasing = T)
    #     helist = helist[ord[1:15],]
    #     noe=dim(helist)[1]
    } else {
        noe=dim(helist)[1]
    }
    
    hel=NULL
    hedge_type=NULL
    for (i in 1:noe){
        he1=c(as.character(helist[i,1]),paste("E",i,sep=''),helist[i,4])
        he2=c(as.character(helist[i,2]),paste("E",i,sep=''),helist[i,4])
        he3=c(as.character(helist[i,3]),paste("E",i,sep=''),helist[i,4])
        if (length(unique(helist[i,1:3]))==1)
            hedge_type=rbind(hedge_type,c(paste("E",i,sep=''),"blue"))
        else if (length(unique(helist[i,1:3]))==2)
            hedge_type=rbind(hedge_type,c(paste("E",i,sep=''),"cyan"))
        else
            hedge_type=rbind(hedge_type,c(paste("E",i,sep=''),"green"))
        hel=rbind(hel,he1,he2,he3)
    }
    rownames(hel)=NULL
    hg=graph.edgelist(hel[,1:2], directed = F)
    hind=match(names(V(hg)),names(V(g)))
    V(hg)$color=V(g)$color[hind]
    for (i in 1:length(hind)){
        if (is.na(hind[i]))
            V(hg)$color[i]=hedge_type[which(hedge_type==names(V(hg)[i]),arr.ind = T)[1],2]
    }
    V(hg)$shape='circle'
    V(hg)$shape[is.na(hind)]='csquare'
    V(hg)$size[is.na(hind)]=2
    hvlabels=rep(NA,length(hind))
    hvlabels[!is.na(hind)]=names(GO.select)[as.numeric(names(V(hg))[!is.na(hind)])]
    hvlabels=as.vector(Term(hvlabels))
    hewgt=as.numeric(hel[,3])
    hewidth=(hewgt-f)/(mean(hewgt)-f)
    
    noc=ncol(label)-1 # num of clusters
    labelSum=apply(label[,1:noc],2,sum)
    hvsizes=labelSum[as.numeric(names(V(hg))[!is.na(hind)])]
    hvsizes=(hvsizes-min(hvsizes))/(mean(hvsizes)-min(hvsizes))
    hvsizes=4 * (hvsizes+2)
    V(hg)$size[!is.na(hind)]=hvsizes
    
    # pdf(paste(term,"_hum_GO (f=",f,").pdf",sep=''))
    # pdf(paste(term,"_yea_GO (f=",f,").pdf",sep=''))
    # set.seed(2)
    plot(hg, # http://igraph.org/r/doc/plot.common.html
         edge.width=hewidth,
         edge.color='grey80',
         vertex.label=hvlabels,
         vertex.label.cex=1,
         vertex.label.color='black',
         vertex.frame.color=NA,
         margin=0)
    
    legend("bottomleft"
           # 0.85,-0.75
           ,c("type1","type2","type3") 
           ,col=c("blue","cyan","green")
           ,pch=15 # square
           # ,lwd=4
           # ,title="legend"
    ) 
    cat("Group ",group,": ",term,"plotted.\n")
}
dev.off()

###############################################################################################
## Plot 3-uniform hypergraph (most connected)

hdegrees=table(elist[,1:3])
nokeep=15 # number of nodes to keep
# nokeep=k # fully kept - CAREFULLY use for supervised
keep=as.numeric(names(sort(hdegrees,decreasing = T)[1:nokeep]))
find_keep=function(x)return(all(x %in% keep))
found_keep=apply(elist[,1:3],1,find_keep)
length(which(found_keep==T))
helist=elist[found_keep,]
noe=dim(helist)[1]
hel=NULL
hedge_type=NULL
for (i in 1:noe){
    he1=c(as.character(helist[i,1]),paste("E",i,sep=''),helist[i,4])
    he2=c(as.character(helist[i,2]),paste("E",i,sep=''),helist[i,4])
    he3=c(as.character(helist[i,3]),paste("E",i,sep=''),helist[i,4])
    if (length(unique(helist[i,1:3]))==1)
        hedge_type=rbind(hedge_type,c(paste("E",i,sep=''),"blue"))
    else if (length(unique(helist[i,1:3]))==2)
        hedge_type=rbind(hedge_type,c(paste("E",i,sep=''),"cyan"))
    else
        hedge_type=rbind(hedge_type,c(paste("E",i,sep=''),"green"))
    hel=rbind(hel,he1,he2,he3)
}
rownames(hel)=NULL
hg=graph.edgelist(hel[,1:2], directed = F)
hind=match(names(V(hg)),names(V(g)))
V(hg)$color=V(g)$color[hind]
for (i in 1:length(hind)){
    if (is.na(hind[i]))
        V(hg)$color[i]=hedge_type[which(hedge_type==names(V(hg)[i]),arr.ind = T)[1],2]
}
V(hg)$shape='circle'
V(hg)$shape[is.na(hind)]='csquare'
V(hg)$size[is.na(hind)]=2
hvlabels=rep(NA,length(hind))
hvlabels[!is.na(hind)]=names(GO.select)[as.numeric(names(V(hg))[!is.na(hind)])]
hvlabels=as.vector(Term(hvlabels))
hewgt=as.numeric(hel[,3])
hewidth=(hewgt-f)/(mean(hewgt)-f)

noc=ncol(label)-1 # num of clusters
labelSum=apply(label[,1:noc],2,sum)
hvsizes=labelSum[as.numeric(names(V(hg))[!is.na(hind)])]
hvsizes=(hvsizes-min(hvsizes))/(mean(hvsizes)-min(hvsizes))
hvsizes=4 * (hvsizes+2)
V(hg)$size[!is.na(hind)]=hvsizes

pdf(paste("core_hum_GO_",nokeep," (f=",f,").pdf",sep=''))
# pdf(paste("core_yea_GO_",nokeep," (f=",f,").pdf",sep=''))
set.seed(1)
plot(hg, # http://igraph.org/r/doc/plot.common.html
     edge.width=hewidth,
     edge.color='grey80',
     vertex.label=hvlabels,
     vertex.label.cex=1,
     vertex.label.color='black',
     vertex.frame.color=NA)

legend("bottomleft"
        # 0.85,-0.75
       ,c("type1","type2","type3")
       ,col=c("blue","cyan","green")
       ,pch=15 # square
       # ,lwd=4
       # ,title="legend"
       ) 
dev.off()


##################################################################################################
##################################################################################################
## Unupervised ##
##################################################################################################

load('net_hum_hclus.bin')
load('label_hum_hclus.bin')

load('net_yea_hclus.bin')
load('label_yea_hclus.bin')

library(igraph)

# pick up those larger than f fold change
k=dim(net)[1] # net is fully normalized 
f=10 # fold change. KEY parameter (should not be higher than 19)
length(which(net>f))
length(which(net[net>0]<=f))
hist(net[net>0],1000,freq=F,xlim=c(0,25))
abline(v=f,col='red')
net_vec_sel=as.vector(net)
net_vec_sel[net_vec_sel<f]=0
net_sel=array(net_vec_sel,c(k,k,k))
elist=which(net_sel>0,arr.ind=TRUE)
ewgt=net_sel[elist]
elist=cbind(elist,ewgt)
V=length(unique(as.vector(elist[,1:3]))) # num of modules that contain selected edges
# save(elist,file='elist_GO.bin')

# use {igraph}
# fold as weight, allow multi-edge
el=rbind(elist[,c(1,2,4)],elist[,c(1,3,4)],elist[,c(2,3,4)])
g=graph.edgelist(apply(el[,1:2],2,as.character), directed = F)
length(V(g))==V
E(g)$weight=el[,3]
# hist(unique(el[,1:2]) # degree distribution for the collapsed graph
     # # ,ylim=c(0,400)
     # )
degrees=table(el[,1:2]) # number of connections for each node
dlevels=sort(unique(degrees),decreasing=T) # degree levels
nocol=length(dlevels) # number of colors
rgb.palette <- colorRampPalette(c("red","orange","yellow" 
                                  # ,"green","cyan","blue"
                                  ), space = "rgb",bias=0.5)
col<-rgb.palette(nocol)
color=rep(NA,V)
ind=match(names(degrees),names(V(g)))
for (i in 1:nocol){
    color[ind[which(degrees==dlevels[i])]]=col[i]
}
V(g)$color=color
# vlabels=names(GO.select)[as.numeric(names(V(g)))]
ewidth=(E(g)$weight-f)/(mean(E(g)$weight)-f)

pdf(paste("full_hum_hclus (f=",f,").pdf",sep=''))
# pdf(paste("full_yea_hclus (f=",f,").pdf",sep=''))
set.seed(03302017)
plot(g, # http://igraph.org/r/doc/plot.common.html
     edge.width=ewidth,
     edge.color='grey75',
     # vertex.label=vlabels,
     vertex.label.cex=0.5,
     vertex.label.color='black',
     vertex.frame.color=NA,
     vertex.size=5)
dev.off()

#####################################################################################
# Plot 3-uniform hypergraph (most connected)


hdegrees=table(elist[,1:3])
nokeep=15 # number of nodes to keep
# nokeep=k # fully kept
keep=as.numeric(names(sort(hdegrees,decreasing = T)[1:nokeep]))
find_keep=function(x)return(all(x %in% keep))
found_keep=apply(elist[,1:3],1,find_keep)
length(which(found_keep==T))
helist=elist[found_keep,]
noe=dim(helist)[1]
hel=NULL
hedge_type=NULL
for (i in 1:noe){
    he1=c(as.character(helist[i,1]),paste("E",i,sep=''),helist[i,4])
    he2=c(as.character(helist[i,2]),paste("E",i,sep=''),helist[i,4])
    he3=c(as.character(helist[i,3]),paste("E",i,sep=''),helist[i,4])
    if (length(unique(helist[i,1:3]))==1)
        hedge_type=rbind(hedge_type,c(paste("E",i,sep=''),"blue"))
    else if (length(unique(helist[i,1:3]))==2)
        hedge_type=rbind(hedge_type,c(paste("E",i,sep=''),"cyan"))
    else
        hedge_type=rbind(hedge_type,c(paste("E",i,sep=''),"green"))
    hel=rbind(hel,he1,he2,he3)
}
rownames(hel)=NULL
hg=graph.edgelist(hel[,1:2], directed = F)
hind=match(names(V(hg)),names(V(g)))
V(hg)$color=V(g)$color[hind]
for (i in 1:length(hind)){
    if (is.na(hind[i]))
        V(hg)$color[i]=hedge_type[which(hedge_type==names(V(hg)[i]),arr.ind = T)[1],2]
}
V(hg)$shape='circle'
V(hg)$shape[is.na(hind)]='csquare'
V(hg)$size[is.na(hind)]=2
hvlabels=rep(NA,length(hind))
hvlabels[!is.na(hind)]=names(V(hg))[!is.na(hind)]
# hvlabels=as.vector(Term(hvlabels))
hewgt=as.numeric(hel[,3])
hewidth=(hewgt-f)/(mean(hewgt)-f)

noc=ncol(label)-1 # num of clusters
labelSum=apply(label[,1:noc],2,sum)
hvsizes=labelSum[as.numeric(names(V(hg))[!is.na(hind)])]
hvsizes=(hvsizes-min(hvsizes))/(mean(hvsizes)-min(hvsizes))
hvsizes=3 * (hvsizes+2)
V(hg)$size[!is.na(hind)]=hvsizes

pdf(paste("core_hum_hclus_",nokeep," (f=",f,").pdf",sep=''))
# pdf(paste("core_yea_hclus_",nokeep," (f=",f,").pdf",sep=''))
set.seed(03292017)
plot(hg, # http://igraph.org/r/doc/plot.common.html
     edge.width=hewidth,
     edge.color='grey80',
     vertex.label=hvlabels,
     vertex.label.cex=1,
     vertex.label.color='black',
     vertex.frame.color=NA)

legend("bottomleft"
        # 0.85,-0.75
       ,c("type1","type2","type3")
       ,col=c("blue","cyan","green")
       ,pch=15 # square
       # ,lwd=4
       # ,title="legend"
       ) 
dev.off()


