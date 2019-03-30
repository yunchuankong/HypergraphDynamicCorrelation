## Author: Yunchuan Kong
## Copyright Reserved 2018

manyZeros<-function(x){
    zeroRate=0.1
    if (length(which(x==min(x)))/length(x)>zeroRate)
    # if (length(which(x==0))/length(x)>zeroRate)
        return(1)
    else
        return(0)
}

deleteZeros <- function(array){
  z=apply(array,1,manyZeros)
  zeroRows=which(z==1)
  if (length(zeroRows)!=0){
    dat=array[-zeroRows,] 
  } else {
    dat=array
  }
  return(dat)
}

normalizeData <- function(dat){
  library(coin)
  rowTrans = function(x)normal_trafo(rank(x, ties.method="random"))
  return(t(apply(dat,1,rowTrans)))
}


penalizeCorr <- function(corr, thres=1){
  v=corr[lower.tri(corr,diag=F)]
  mu_v=mean(v)
  sd_v=sd(v)
  dict=ifelse((corr>mu_v-1*sd_v & corr<mu_v+1*sd_v) , 1, 0)
#   cat("Correlation thresholds:",c(mu_v-thres*sd_v,mu_v+thres*sd_v),"\n")
  return(dict)
}

one_gene_LA<-function(array, x){ # x is a row of the array, i.e. nog*nos
    b<-x*t(array)
    d<-array%*%b
    d 
    # diaganol elements and elements of x.(index) in array are of no use
    # the matrix is symmetric
    # the result should be divided by nos to get LAs
}
 
extractVector<-function(mat){ # extracting the upper triangular into a vector
    return( as.vector(mat[upper.tri(mat, diag=F)]) )
}

pfdrBounds<-function(i,dat,nos,one_vec,naRec, fdr=0.1){
    library(ks)
	  threshold_pfdr = fdr
    n_sample_pairs = 1000
    n_perm = 1
    pairID = sample(which(naRec==FALSE),n_sample_pairs)
    ## Below is not fast. Use the for loop
    # perm = replicate(n_perm,one_gene_LA(dat,sample(dat[i,]))/nos)
    # perm_vec = matrix(unlist(alply(perm,3,extractVector)),ncol=n_perm)
    # perm_vec = perm_vec[pairID,]
    perm_vec=NULL
    for (j in 1:n_perm){ # use this loop when n_perm>1
        perm=one_gene_LA(dat,sample(dat[i,]))/nos
        perm_vec=extractVector(perm)[pairID]
        perm_vec=c(perm_vec,perm_vec)
    }
    est_true = kde(one_vec[pairID]
                   ,eval.points = sort(one_vec[pairID])
                   )
    dens_true = est_true$estimate
    ep = est_true$eval.points
    est_perm = kde(perm_vec
                   ,eval.points = ep
                   )
    dens_perm = est_perm$estimate

    ratio = dens_perm/dens_true
    idx_low=which(ratio[1:(n_sample_pairs/2)]<threshold_pfdr)
    idx_upp=which(ratio[-(1:(n_sample_pairs/2))]<threshold_pfdr)
    if (length(idx_low)==0){
        lb = NA
    } else
        lb = ep[max(idx_low)]
    if (length(idx_upp)==0){
        ub = NA
    } else
        ub = ep[min(idx_upp+n_sample_pairs/2)]
    perm=NULL;perm_vec=NULL;est_perm=NULL;dens_perm=NULL
    est_true=NULL;dens_true=NULL;ratio=NULL
    return(c(lb,ub))
}

make3<-function(x,i,nog){ 
# LA will not be recorded in triplet list any more for saving space (can be modified)
    tri=NULL
    not=length(x[!(is.na(x))]) # Number Of Triplets 
    x[!(is.na(x))]=1
	  x[is.na(x)]=0
    tuples=rep(0,2*not)
    tuples=.C('make3',tuples=as.integer(tuples)
                     ,x=as.integer(x)
                     ,theOne=as.integer(i)
                     ,nog=as.integer(nog)
					 ,not=as.integer(not))$tuples
    tuples=matrix(tuples,ncol=2,nrow=not,byrow=T) 
    tuples=tuples[apply(tuples,1,sum)!=0,]    
    if (i %% 100==0)
        cat("Number of triplets:",not,'.    ')
    if ( length(tuples)!=2 ){
        not=nrow(tuples)
        tri=cbind(rep(i,not),tuples)
    }
    else {
        # cat("Hey! I found you!\n")
        tri=c(i,tuples)
    }
    # free the memory
    tuples=NULL
    return(tri)
}

makeTriplets <- function(nog, nos, dat, dict, fdr){
  triplets <- NULL
  for (i in 1:nog ){
  # for (i in 1:15){
    dicti=dict[i,]*t(dict[i,]*dict)
    one=one_gene_LA(dat,dat[i,])/nos
    one[dicti==0]=NA
    one_vec = extractVector(one)    
    naRec = is.na(one_vec)
    bounds = pfdrBounds(i,dat,nos,one_vec,naRec, fdr=0.1)
    
    if (is.na(bounds[1]) & !is.na(bounds[2])){
      one_vec[ one_vec<bounds[2] ] = NA
    }
    else if (!is.na(bounds[1]) & is.na(bounds[2])){
      one_vec[ one_vec>bounds[1] ] = NA
    } 
    else if (!is.na(bounds[1]) & !is.na(bounds[2])) {
      one_vec[ one_vec>bounds[1] & one_vec<bounds[2] ]=NA
    }
    else {
      # cat("Gene:", i, "No effective triplets - passed.\n")
      next
    }
    
    triplets=rbind(triplets,make3(one_vec,i,nog))
    if (i %% 100==0)
      cat("*** Gene ",i," calculated. ***\n")
    # free the memory
    one=NULL;one_vec=NULL;naRec=NULL;
    # subOne=NULL;x=NULL
  } # for
  cat("\n**** Finished calculating assigned genes. ****\n")
  return(triplets)
}

count<-function(triplet,N){
	n=nrow(triplet)
	flow=rep(0,N*N)
	.C("count",tri=as.integer(as.vector(t(triplet))),
		n=as.integer(n),
		flow=as.integer(flow),
		N=as.integer(N))$flow
}

createFlow<-function(triplets, nog, save=F, corMat=F){
    k=nog    
    flow_vec=count(triplets,k)
    flow=matrix(flow_vec,k,k)
    if (save==TRUE)
      save(flow,file="flow.bin")
    return(flow)
}

createLabels_unsup <- function(distMat, save=F){ 
# return the one-hot label matrix, with one additional column recording the sum of each row
  library(dynamicTreeCut)
  N=ncol(distMat)
  d=as.dist(distMat)
  h=hclust(d, method="average")
  label_vec=cutreeDynamic(h,distM = distMat,verbose=4,
                                    minClusterSize = 20, # let's do 20 for yeast, 100 for human
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
  if (save==TRUE){
    save(label,file='label_unsup.bin')
  }
	return(label)
} 
  
createLabels_sup <- function(GO.select, glist, save=F){ 
## input a list object GO.select 
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
  if (save==TRUE){
    save(label,file='label_sup.bin')
  }
  return(label)
}

shrink<-function(triplet,label){
  n=dim(triplet)[1]    
	K=ncol(label) # num of clusters + 1
	adjacency=rep(0,(K-1)^3)
  .C('shrink',tri=as.integer(as.vector(t(triplet))),
              n=as.integer(n),
              adj=as.integer(adjacency),
              lab=as.integer(as.vector(t(label))),
              k=as.integer(K))$adj
}

fullyConnect<-function(labelSum){
	K=length(labelSum) # num of clusters
	adjacency=rep(0,(K)^3)
	.C("fullyConnect",adj=as.integer(adjacency),
	          			  labSum=as.integer(labelSum),
					          k=as.integer(K))$adj
}

normalize<-function(noc,net_vec,fullNet_vec){
  K=noc # num of clusters
  adjacency=rep(0,(K)^3)
  .C("normalize", net=as.double(net_vec),
                  full=as.double(fullNet_vec),
                  k=as.integer(K),
                  adj=as.double(adjacency))$adj
}

GOen<-function(glist,label){
    library(GOstats)
    all.entrez=glist
    res=list()
    k=ncol(label)-1
    for (i in 1:k){
        sel.entrez=glist[which(label[,i]==1)]
        params <- new("GOHyperGParams", geneIds=sel.entrez
                        ,universeGeneIds=all.entrez, ontology="BP"
                        ,pvalueCutoff=0.01,conditional=F 
                        ,testDirection="over"
                        ,annotation="org.Sc.sgd.db") # yeast
                        # ,annotation="org.Hs.eg.db") # human/skin
        over.pres<-hyperGTest(params)
        sum = summary(over.pres)
        id = sum$Size %in% 5:500
        res[[i]] <- sum[id,]
        cat(i,"out of",k,'done.\n')
    }
    return(res)
}

