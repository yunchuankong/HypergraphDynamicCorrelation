# setwd("C:\\Users\\DELL\\Dropbox\\Research_Yu\\elves")
# setwd("C:\\Users\\kyccw\\Dropbox\\Research_Yu\\elves")
dyn.load('csupp.so')
load("data.bin") # Human gene dataset
library(fdrtool)

one.gene.LA<-function(array, x){ # x is a row of the array, i.e. nog*nos
    b<-x*t(array)
    d<-array%*%b
    d 
    # diaganol elements and elements of x.(index) in array are of no use
    # the matrix is symmetric
    # the result should be divided by nos to get LAs
}

lfdrBounds<-function(x){
    threshold_lfdr = 0.001
    x_sample=x[ sample(length(x),5000) ]
    lfdr_sample=fdrtool(x_sample,statistic="normal",cutoff.method="fndr",plot=F,verbose=F)
    y_sample=x_sample[lfdr_sample$lfdr<threshold_lfdr]
    lb=max(y_sample[y_sample<0])
    ub=min(y_sample[y_sample>0])
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
	return(tri)
}

triplet=NULL
for (i in 1:(nog-2)){
# for (i in 1:15){
    one=one.gene.LA(dat,dat[i,])/nos
	one[dict==0]=NA
    one_vec=as.vector(one[upper.tri(one, diag=F)])
	naRec=is.na(one_vec)
    bounds=lfdrBounds(one_vec[naRec==F])
    subOne=one[(i+1):nog,(i+1):nog]
    x=as.vector(subOne[upper.tri(subOne, diag=F)])
    x[ !(x<bounds[1] | x>bounds[2]) ]=NA
	# cat("not of x:",length(x[!(is.na(x))]),'\n')
    if( all( is.na(x) ) ){
        cat("Gene:", i, "No effective triplets - passed.\n")
        next
    }
    one_vec[ !(one_vec<bounds[1] | one_vec>bounds[2]) ]=NA
    triplet=rbind(triplet,make3(one_vec,i,nog))
	if (i %% 100==0)
        cat("*** Gene ",i," calculated. ***\n")
} # for
cat("\n**** Finished calculating assigned genes. ****\n")
save(triplet
    # ,nog,glist
    ,file="triplet_hum.bin")

# cat("\nTriplet size before unique: ", nrow(triplet1),"\n")
# triplet=unique(triplet1)
# cat("\nTriplet size after unique: ", nrow(triplet),"\n")
# save(triplet1,file='triplet_unique.bin')
# dup=duplicated(triplet1)
# print(triplet1[dup,])
# print(triplet1[1:20,])

# cat(range(triplet),'\n')
# print(length(unique(as.vector(triplet))))




