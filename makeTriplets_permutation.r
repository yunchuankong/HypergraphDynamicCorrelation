# setwd("C:\\Users\\DELL\\Dropbox\\Research_Yu\\elves")
# setwd("C:\\Users\\kyccw\\Dropbox\\Research_Yu\\elves")
dyn.load('csupp.so')
load("data_hum.bin") # Human gene dataset
# load("data_yea.bin")
# library(fdrtool)
# library(plyr)
library(ks)

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

# lfdrBounds<-function(x){
#     threshold_lfdr = 0.001
#     x_sample=x[ sample(length(x),5000) ]
#     lfdr_sample=fdrtool(x=x_sample,plot=F,verbose=F)
#     y_sample=x_sample[lfdr_sample$lfdr<threshold_lfdr]
#     
#     length(y_sample)
#     range(lfdr_sample$lfdr)
#     
#     lb=max(y_sample[y_sample<0])
#     ub=min(y_sample[y_sample>0])
#     return(c(lb,ub))
# }

pfdrBounds<-function(i,dat,nos,one_vec,naRec){
	threshold_pfdr = 0.01
    n_sample_pairs = 1000
    n_perm = 1
    pairID = sample(which(naRec==FALSE),n_sample_pairs)
    ## Below is not fast. Use the for loop
    # perm = replicate(n_perm,one_gene_LA(dat,sample(dat[i,]))/nos)
    # perm_vec = matrix(unlist(alply(perm,3,extractVector)),ncol=n_perm)
    # perm_vec = perm_vec[pairID,]
    perm_vec=NULL
    # perm=one_gene_LA(dat,sample(dat[i,]))/nos
    # perm_vec=extractVector(perm)[pairID]
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
    
    # plot(est_perm,cont=T)
    # lines(ep,dens_true,type='l',col='red')
    
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
    # view=cbind(ep,dens_perm,dens_true,ratio)
    # View(round(view,3))
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

triplet1=NULL
bin=(nog-2)/5
for (i in 1:as.integer(bin) ){
# for (i in 1:nog ){
# for (i in 1:15){
    dicti=dict[i,]*t(dict[i,]*dict)
    one=one_gene_LA(dat,dat[i,])/nos
    one[dicti==0]=NA
    one_vec = extractVector(one)    
    naRec = is.na(one_vec)
    bounds = pfdrBounds(i,dat,nos,one_vec,naRec)
    # cat("Bounds:",bounds,"\n")

#   # screen previous i to avoid duplicated triplet
#     subOne=one[(i+1):nog,(i+1):nog]
#     x=as.vector(subOne[upper.tri(subOne, diag=F)])
# 	# cat("range x:",range(x[!(is.na(x))]),"\n")
#     x[ !(x<bounds[1] | x>bounds[2]) ]=NA
# 	# cat("not of x:",length(x[!(is.na(x))]),'\n')
#     if( all( is.na(x) ) ){
#         cat("Gene:", i, "No effective triplets - passed.\n")
#         next
#     }
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
        cat("Gene:", i, "No effective triplets - passed.\n")
        next
    }
    
    triplet1=rbind(triplet1,make3(one_vec,i,nog))
	if (i %% 100==0)
        cat("*** Gene ",i," calculated. ***\n")
    # free the memory
    one=NULL;one_vec=NULL;naRec=NULL;
    # subOne=NULL;x=NULL
} # for
cat("\n**** Finished calculating assigned genes. ****\n")
save(triplet1,file="GMG/ykong24/triplet1_hum.bin")




