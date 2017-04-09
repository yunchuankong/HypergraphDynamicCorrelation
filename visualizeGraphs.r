# setwd("C:\\Users\\DELL\\Dropbox\\Research_Yu\\darkness")
# setwd("C:\\Users\\kyccw\\Dropbox\\Research_Yu\\darkness")

####################################################################################
# Supervised

load('net_GO.bin')
load('GO_select_yeast_ 0.7 0.4 1000 20 .bin')

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
f=10 # fold change. KEY parameter (should not be higher than 19)
length(which(net>f))
# length(which(net[net>0]<=f))
# hist(net[net>0],1000,freq=F,xlim=c(0,25))
# abline(v=f,col='red')
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
# pdf(paste("full_GO (f=",f,").pdf",sep=''))
plot(g, # http://igraph.org/r/doc/plot.common.html
     edge.width=ewidth,
     edge.color='grey75',
     vertex.label=vlabels,
     vertex.label.cex=0.3,
     vertex.label.color='black',
     vertex.frame.color=NA,
     vertex.size=5)
# dev.off()

# Plot 3-uniform hypergraph (partial)
hdegrees=table(elist[,1:3])
nokeep=4 # number of nodes to keep
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
V(hg)$size=10
V(hg)$size[is.na(hind)]=2
hvlabels=rep(NA,length(hind))
hvlabels[!is.na(hind)]=names(GO.select)[as.numeric(names(V(hg))[!is.na(hind)])]
hvlabels=as.vector(Term(hvlabels))
hewgt=as.numeric(hel[,3])
hewidth=(hewgt-f)/(mean(hewgt)-f)

pdf(paste("core_GO_",nokeep,".pdf",sep=''))
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


#####################################################################################
# Unupervised

load('net_hclus.bin')

library(igraph)

# pick up those larger than f fold change
k=dim(net)[1] # net is fully normalized 
f=10 # fold change. KEY parameter (should not be higher than 19)
length(which(net>f))
# length(which(net[net>0]<=f))
# hist(net[net>0],1000,freq=F,xlim=c(0,25))
# abline(v=f,col='red')
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

# pdf(paste("full_hclus (f=",f,").pdf",sep=''))
set.seed(03302017)
plot(g, # http://igraph.org/r/doc/plot.common.html
     edge.width=ewidth,
     edge.color='grey75',
     # vertex.label=vlabels,
     vertex.label.cex=0.5,
     vertex.label.color='black',
     vertex.frame.color=NA,
     vertex.size=5)
# dev.off()

# Plot 3-uniform hypergraph (partial)
hdegrees=table(elist[,1:3])
nokeep=33 # number of nodes to keep
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
V(hg)$size=10
V(hg)$size[is.na(hind)]=2
hvlabels=rep(NA,length(hind))
hvlabels[!is.na(hind)]=names(V(hg))[!is.na(hind)]
# hvlabels=as.vector(Term(hvlabels))
hewgt=as.numeric(hel[,3])
hewidth=(hewgt-f)/(mean(hewgt)-f)

pdf(paste("core_hclus_",nokeep,".pdf",sep=''))
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


