
# dyn.load(paste0(filepath, "csupp.dll")) ## for Windows only
# source(paste0(filepath,"utils.r"))
# source(paste0(filepath,"visualize.r")) 

dyn.load("csupp.dll")
source("utils.r")
source("visualize.r")

## Given data matrix "array" and supervised labels "GO.select"

hypergraph_unsup <- function(array, fdr=0.2, save_triplets=FALSE, save_glist=FALSE){

    dat <- deleteZeros(array)
    dat <- normalizeData(dat)

    nog=nrow(dat) ## number of genes
    nos=ncol(dat) ## number of samples
    glist=rownames(dat) ## gene name list
    clist=colnames(dat) ## sample name list
    rownames(dat)=1:nog ## re-name genes as numbers
    colnames(dat)=NULL

    corr=cor(t(dat))  ## pairwise correlation matrix
    dict <- penalizeCorr(corr)  ## "dictionary" for low correlated pairs
    diag(dict) <- 0

    triplets <- makeTriplets(nog, nos, dat, dict, fdr=fdr) 
    # dim(triplets)

    A <- createFlow(triplets, nog, corMat=FALSE) ## usually don't calculate corMat here
    glist_new <- which(apply(A,2,sum)!=0)
    # cat(length(glist_new),"out of", dim(A)[1],".\n")
    corA <- cor(A[glist_new,glist_new]) ## calculate the correlation matrix for A here
    label_sub <- createLabels_unsup(1-corA) 
    label <- matrix(0,ncol=ncol(label_sub),nrow=nrow(A))
    label[glist_new,] <- label_sub

    noc=ncol(label)-1 ## num of clusters
    labelSum <- apply(label[,1:noc],2,sum)
    net_vec<-shrink(triplets,label)
    fullNet_vec<-fullyConnect(labelSum)
    net<-normalize(noc,net_vec,fullNet_vec)
    r=sum(as.numeric(net_vec))/sum(as.numeric(fullNet_vec))
    net <- net/r
    net <- array(net,c(noc,noc,noc)) 
    
    res_unsup <- NULL
    res_unsup$net <- net
    res_unsup$label <- label
#     res_unsup$enrichment <- enrichment    
    if (isTRUE(save_triplets)){
        res_unsup$triplets <- triplets
    }
    if (isTRUE(save_glist)){
        res_unsup$glist <- glist
    }
    return(res_unsup)
}

hypergraph_sup <- function(array, label_list, fdr=0.2, save_triplets=FALSE, save_glist=FALSE){

    dat <- deleteZeros(array)
    dat <- normalizeData(dat)

    nog=nrow(dat) ## number of genes
    nos=ncol(dat) ## number of samples
    glist=rownames(dat) ## gene name list
    clist=colnames(dat) ## sample name list
    rownames(dat)=1:nog ## re-name genes as numbers
    colnames(dat)=NULL

    corr=cor(t(dat))  ## pairwise correlation matrix
    dict <- penalizeCorr(corr)  ## "dictionary" for low correlated pairs
    diag(dict) <- 0

    triplets <- makeTriplets(nog, nos, dat, dict, fdr=fdr)
    
    label <- createLabels_sup(label_list, glist)

    noc=ncol(label)-1 ## num of clusters
    labelSum=apply(label[,1:noc],2,sum)
    net_vec<-shrink(triplets,label)
    fullNet_vec<-fullyConnect(labelSum)
    net<-normalize(noc,net_vec,fullNet_vec)
    r <- sum(as.numeric(net_vec))/sum(as.numeric(fullNet_vec))
    net <- net/r
    net <- array(net,c(noc,noc,noc))
    
    res_sup <- NULL
    res_sup$net <- net
    res_sup$label <- label
    if (isTRUE(save_triplets)){
        res_sup$triplets <- triplets
    }
    if (isTRUE(save_glist)){
        res_sup$glist <- glist
    }
    return(res_sup)
}    

