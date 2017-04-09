# setwd("C:\\Users\\DELL\\Dropbox\\Research_Yu\\darkness")
# setwd("C:\\Users\\kyccw\\Dropbox\\Research_Yu\\darkness")

dyn.load('csupp.so')
# library(pheatmap)

# load("triplet_unselected.bin")
# cat("\nTriplets unselected loaded.\n")
# load('elist_GO.bin')

count<-function(triplet,N){
	n=nrow(triplet)
	flow=rep(0,N*N)
	.C("count",tri=as.integer(as.vector(t(triplet))),
		n=as.integer(n),
		flow=as.integer(flow),
		N=as.integer(N))$flow
}

createFlow<-function(triplet,save=F){
    k=length(unique(as.vector(triplet))) 
    flow_vec=count(triplet,k)
    flow=matrix(flow_vec,k,k)
    corFlow=cor(flow)
    if (save==TRUE)
        save(flow,corFlow,file="flow.bin")
    return(list("flow"=flow,"corFlow"=corFlow))
}

createLabel<-function(distMat=NULL,method='Hierarchical',save=F){ 
# return the one-hot label matrix, with one additional column recording the sum of each row
    if (method=='Hierarchical'){
        if (distMat==NULL){
            cat("Error: no distance matrix provided.\n")
            return(NULL)
        }
        cat("Hierarchical clustering with dynamic tree cut.\n")
        library(dynamicTreeCut)
        N=ncol(distMat)
        d=as.dist(distMat)
        h=hclust(d,method = 'average')
        label_vec=cutreeDynamic(h,distM = distMat,verbose=4,
                                        minClusterSize = 20,
                                        deepSplit = 1,
                                        pamStage =F)
		K=length(unique(label_vec))-1		
		cat("Number of clusters: ",K,"\n")
        label=matrix(0,ncol=K,nrow=N)
        for (i in 1:K){
			label[which(label_vec==i),i]=1
        }
        label=cbind(label,apply(label,1,sum))
        colnames(label)=c(1:(ncol(label)-1),'rowSum')
        if (save==TRUE)
            save(label,file='label.bin')
		return(label)
    }
    else if (method=='GO'){
        load("spellman_73_filled.bin")
		glist=rownames(array)
		load("GO_select_yeast_ 0.7 0.4 1000 20 .bin")
        N=length(glist)
		K=length(GO.select)
		label=matrix(0,ncol=K,nrow=N)
		for (i in 1:K){
			label[match(GO.select[[i]],glist),i]=1
			
			# len_act=length(which(label[,i]!=0))
			# len_tgt=length(GO.select[[i]])
			# if (len_act!=len_tgt)
				# cat("WARNING: length inconsistent. Module:",i,len_act,len_tgt,"\n")
		}
        label=cbind(label,apply(label,1,sum))
        colnames(label)=c(1:(ncol(label)-1),'rowSum')
        if (save==TRUE)
            save(label,file='label.bin')
		return(label)
    }
    else {
        return(NULL)
    }
}

# corFlow=createFlow(triplet[,1:3])$corFlow
# label=createLabel(1-corFlow,method='Hierarchical')
label=createLabel(method='GO')
# dim(which(label!=0,arr.ind=T))[1]




