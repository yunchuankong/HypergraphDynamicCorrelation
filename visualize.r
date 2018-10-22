## Author: Yunchuan Kong
## Copyright Reserved 2018

library(igraph)

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

##############################################################################################
## Supervised section

## INPUT: net, label, GO.select
plot_entrie_sup <- function(net, label, GO.select, folds=10){

  # pick up those larger than f fold change
  k=dim(net)[1] # net is fully normalized 
  f=folds # fold change
  net_vec_sel=as.vector(net)
  net_vec_sel[net_vec_sel<f]=0
  net_sel=array(net_vec_sel,c(k,k,k))
  elist=which(net_sel>0,arr.ind=TRUE)
  ewgt=net_sel[elist]
  elist=cbind(elist,ewgt)
  V=length(unique(as.vector(elist[,1:3]))) # num of modules that contain selected edges
  # print(V)
  
  # fold as weight, allow multi-edge
  el=rbind(elist[,c(1,2,4)],elist[,c(1,3,4)],elist[,c(2,3,4)])
  g=graph.edgelist(apply(el[,1:2],2,as.character), directed = F)
  E(g)$weight=el[,3]
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

  # ewidth=(E(g)$weight-f)/(mean(E(g)$weight)-f)
  ewidth = log(E(g)$weight-f+1)

  # pdf(paste("full_sup (f=",f,").pdf",sep=''))
  plot(g, # http://igraph.org/r/doc/plot.common.html
      edge.width=ewidth,
      edge.color='grey75',
      vertex.label=vlabels,
      vertex.label.cex=0.3,
      vertex.label.color='black',
      vertex.frame.color=NA,
      vertex.size=5)
  # dev.off()
  # cat("A full hypergraph plot (2-D version) for the supervised approach has been saved.\n")
  res <- NULL
  res$elist <- elist
  res$g <- g
  res$folds <- f
  return(res)
}


## INPUT: g, elist from plot_entire_sup(), and folds should be consistent
plot_top_sup <- function(g, elist, GO.select, folds, n_keep=15){
  library(GO.db)
  f = folds
  hdegrees=table(elist[,1:3])
  nokeep=n_keep # number of nodes to keep
  keep=as.numeric(names(sort(hdegrees,decreasing = T)[1:nokeep]))
  find_keep=function(x)return(all(x %in% keep))
  found_keep=apply(elist[,1:3],1,find_keep)
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
  V(hg)$size[is.na(hind)]=4
  hvlabels=rep(NA,length(hind))
  hvlabels[!is.na(hind)]=names(GO.select)[as.numeric(names(V(hg))[!is.na(hind)])]
  hvlabels=as.vector(Term(hvlabels))
  hewgt=as.numeric(hel[,3])
  hewidth=log(hewgt-f+1)
  noc=ncol(label)-1 # num of clusters
  labelSum=apply(label[,1:noc],2,sum)
  hvsizes=labelSum[as.numeric(names(V(hg))[!is.na(hind)])]
  hvsizes=(hvsizes-min(hvsizes))/(mean(hvsizes)-min(hvsizes))
  hvsizes=4 * (hvsizes+2)
  V(hg)$size[!is.na(hind)]=hvsizes
  
  # pdf(paste("Top",nokeep,"Connected_sup (f=",f,").pdf",sep=''))
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
  # dev.off()
  # cat("A hypergraph plot of the",nokeep,"most connected vertices for the supervised approach has been saved.\n")
}

## INPUT: g, elist from plot_entire_sup(), and folds should be consistent
plot_one_sup <- function(g, elist, GO.select, folds, GOID){
  library(GO.db)
  f = folds
  term = Term(GOID)
  groupID = which(names(GO.select)==GOID)
  find_keep=function(x)return(all(groupID %in% x))
  found_keep=apply(elist[,1:3],1,find_keep)
  if (length(which(found_keep==T))==0){
    cat("WARNING: the group is not involved under the current fold change threshold!\n")
    return(NULL)
  }
  helist=elist[found_keep,]
  if (length(which(found_keep==T))==1){
    cat("WARNING: the group has only one edge. No hypergraph will be drawn.\n")
    cat("The edge connects:\n",Term(names(GO.select)[helist[1]]),"\n",
        Term(names(GO.select)[helist[2]]),"\n",
        Term(names(GO.select)[helist[3]]),"\n")
    return(NULL)
  # } 
  # ## If want a limited number of connections per plot: (15)
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
  V(hg)$size[is.na(hind)]=6
  hvlabels=rep(NA,length(hind))
  hvlabels[!is.na(hind)]=names(GO.select)[as.numeric(names(V(hg))[!is.na(hind)])]
  hvlabels=as.vector(Term(hvlabels))
  hewgt=as.numeric(hel[,3])
  # hewidth=(hewgt-f)/(mean(hewgt)-f)
  hewidth=log(hewgt-f+1)
  
  noc=ncol(label)-1 # num of clusters
  labelSum=apply(label[,1:noc],2,sum)
  hvsizes=labelSum[as.numeric(names(V(hg))[!is.na(hind)])]
  hvsizes=(hvsizes-min(hvsizes))/(mean(hvsizes)-min(hvsizes))
  hvsizes=4 * (hvsizes+2)
  V(hg)$size[!is.na(hind)]=hvsizes

  # pdf(paste(term,"_sup (f=",f,").pdf",sep=''))
  plot(hg, # http://igraph.org/r/doc/plot.common.html
       edge.width=hewidth,
       edge.color='grey80',
       vertex.label=hvlabels,
       vertex.label.cex=1,
       vertex.label.color='black',
       vertex.frame.color=NA,
       margin=c(0,0,0,0)
  )
  
  legend("bottomleft"
         # 0.85,-0.75
         ,c("type1","type2","type3") 
         ,col=c("blue","cyan","green")
         ,pch=15 # square
         # ,lwd=4
         # ,title="legend"
  ) 
  # dev.off()
  # cat("A hypergraph plot of",GOID,"and its connected vertices has been saved.\n")
}

##############################################################################################
## plot gene level hypergraph for a given module level hyper-edge

query <- function(triplet, c1, c2, c3){
  n=dim(triplet)[1]    
  k1 <- length(c1)
  k2 <- length(c2)
  k3 <- length(c3)
  res <- rep(c(0,0,0), k1*k2*k3)
  .C('query',tri=as.integer(as.vector(t(triplet))),
     n=as.integer(n),
     c1=as.integer(c1),
     c2=as.integer(c2),
     c3=as.integer(c3),
     k1=as.integer(k1),
     k2=as.integer(k2),
     k3=as.integer(k3),
     res=as.integer(res))$res
}

find_Vhg <- function(str){
  if (unlist(strsplit(str,""))[1] == "E" ){
    return(NA)
  } else {
    return(1)
  }
}

plot_gene_level <- function(hyperedge, module_names, 
                            net, label, triplets, glist, 
                            folds){
  net_sel <- net
  net_sel[net<folds] <- 0
  hyperedge <- sort(hyperedge, decreasing=T)
  if (net_sel[hyperedge[1],hyperedge[2],hyperedge[3]] == 0){
    cat("The hyperedge does not exist under the current fold change.\n")
  }
  
  c1 <- label[,hyperedge[1]]
  c1 <- which(c1==1)
  c2 <- label[,hyperedge[2]]
  c2 <- which(c2==1)
  c3 <- label[,hyperedge[3]]
  c3 <- which(c3==1)
  
  sel_list <- query(triplets, c1, c2, c3) 
  sel_list <- sel_list[sel_list!=0]
  sel_list <- matrix(sel_list, ncol=3, byrow=T)

  helist <- sel_list
  noe=dim(helist)[1]
  hel=NULL
  hedge_type=NULL
  for (j in 1:noe){
    he1=c(as.character(helist[j,1]),paste("E",j,sep=''))
    he2=c(as.character(helist[j,2]),paste("E",j,sep=''))
    he3=c(as.character(helist[j,3]),paste("E",j,sep=''))
    if (length(unique(hyperedge))==1)
      hedge_type=rbind(hedge_type,c(paste("E",j,sep=''),adjustcolor("blue", alpha.f = .5)))
    else if (length(unique(hyperedge))==2)
      hedge_type=rbind(hedge_type,c(paste("E",j,sep=''),adjustcolor("cyan", alpha.f = .5)))
    else
      hedge_type=rbind(hedge_type,c(paste("E",j,sep=''),adjustcolor("green", alpha.f = .5)))
    hel=rbind(hel,he1,he2,he3)
  }
  rownames(hel)=NULL
  hg=graph.edgelist(hel[,1:2], directed = F)
  hind <- apply(as.matrix(names(V(hg))),1,find_Vhg)
  
  for (k in 1:length(hind)){
    if (is.na(hind[k])){
      V(hg)$color[k]=hedge_type[which(hedge_type==names(V(hg)[k]),arr.ind = T)[1],2]
    } else if (hyperedge[1] %in% which(label[as.numeric(names(V(hg)[k])),]==1)){
      V(hg)$color[k] <- adjustcolor("purple1", alpha.f = .5)
    } else if (hyperedge[2] %in% which(label[as.numeric(names(V(hg)[k])),]==1)){
      V(hg)$color[k] <- adjustcolor("orchid1", alpha.f = .5)
    } else {
      V(hg)$color[k] <- adjustcolor("navajowhite", alpha.f = .5) 
    }
  } ## for k
  V(hg)$shape='circle'
  V(hg)$shape[is.na(hind)]='csquare'
  V(hg)$size[is.na(hind)]=1
  # V(hg)$size[!is.na(hind)]=6
  hvlabels=rep(NA,length(hind))
  hvlabels[!is.na(hind)]=glist[as.numeric(names(V(hg))[!is.na(hind)])]
  
  noc=ncol(label)-1 # num of clusters
  freq <- table(helist)
  hvsizes=freq[match(names(V(hg))[!is.na(hind)], names(freq))]
  hvsizes=(hvsizes-min(hvsizes))/(mean(hvsizes)-min(hvsizes))
  hvsizes=3 * (hvsizes+1)
  V(hg)$size[!is.na(hind)]=hvsizes
  
  # pdf(paste0("graph_gene_level"))
  # layout <- layout_with_kk(hg)
  # set.seed(2)
  plot(hg, # http://igraph.org/r/doc/plot.common.html
       layout=layout_with_lgl(hg), # layout_on_sphere(hg),
       edge.width=0.4,
       edge.color='grey80',
       vertex.label=hvlabels,
       vertex.label.cex=0.8,
       vertex.label.color=adjustcolor("black", alpha.f = .5),
       vertex.frame.color=NA,
       margin=c(0,0,0,0)
  )
  legend("bottomleft"
         # 0.85,-0.75
         ,module_names[hyperedge]
         ,col=c("purple1","orchid1","navajowhite")
         ,pch=15 # square
         # ,lwd=4
         # ,title="legend"
  ) 
  # dev.off()
}



##############################################################################################
## Unsupervised section

## INPUT: label, net
plot_entrie_unsup <- function(net, label, folds=10){
  # pick up those larger than f fold change
  k=dim(net)[1] 
  f=folds 
  net_vec_sel=as.vector(net)
  net_vec_sel[net_vec_sel<f]=0
  net_sel=array(net_vec_sel,c(k,k,k))
  elist=which(net_sel>0,arr.ind=TRUE)
  ewgt=net_sel[elist]
  elist=cbind(elist,ewgt)
  V=length(unique(as.vector(elist[,1:3]))) # num of modules that contain selected edges
  # print(V)
  
  # fold as weight, allow multi-edge
  el=rbind(elist[,c(1,2,4)],elist[,c(1,3,4)],elist[,c(2,3,4)])
  g=graph.edgelist(apply(el[,1:2],2,as.character), directed = F)
  E(g)$weight=el[,3]
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
  # ewidth=(E(g)$weight-f)/(mean(E(g)$weight)-f)
  ewidth = log(E(g)$weight-f+1)
  
  # pdf(paste("full_unsup (f=",f,").pdf",sep=''))
  plot(g, # http://igraph.org/r/doc/plot.common.html
       edge.width=ewidth,
       edge.color='grey75',
       # vertex.label=vlabels,
       vertex.label.cex=0.5,
       vertex.label.color='black',
       vertex.frame.color=NA,
       vertex.size=5)
  # dev.off()
  # cat("A full hypergraph plot (2-D version) for the unsupervised approach has been saved.\n")
  res <- NULL
  res$elist <- elist
  res$g <- g
  res$folds <- f
  return(res)
}


## INPUT: g, elist from plot_entire_unsup(), and folds should be consistent
plot_top_unsup <- function(g, elist, folds, n_keep=15){
  f = folds
  hdegrees=table(elist[,1:3])
  nokeep=n_keep # number of nodes to keep
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
  V(hg)$size[is.na(hind)]=5
  hvlabels=rep(NA,length(hind))
  hvlabels[!is.na(hind)]=names(V(hg))[!is.na(hind)]
  # hvlabels=as.vector(Term(hvlabels))
  hewgt=as.numeric(hel[,3])
  # hewidth=(hewgt-f)/(mean(hewgt)-f)
  hewidth = log(hewgt-f+1)
  noc=ncol(label)-1 # num of clusters
  labelSum=apply(label[,1:noc],2,sum)
  hvsizes=labelSum[as.numeric(names(V(hg))[!is.na(hind)])]
  hvsizes=(hvsizes-min(hvsizes))/(mean(hvsizes)-min(hvsizes))
  hvsizes=6 * (hvsizes+2)
  V(hg)$size[!is.na(hind)]=hvsizes
  
  # pdf(paste("Top",nokeep,"Connected_unsup (f=",f,").pdf",sep=''))
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
  # dev.off()
  # cat("A hypergraph plot of the",nokeep,"most connected vertices for the unsupervised approach has been saved.\n")
}












