# setwd("C:\\Users\\DELL\\Dropbox\\Research_Yu\\darkness")
# setwd("C:\\Users\\kyccw\\Dropbox\\Research_Yu\\darkness")
dyn.load('csupp.so')
load("triplet_unselected.bin")
load("label_hclus.bin")

shrink<-function(triplet,label){
    n=dim(triplet)[1]    
	K=ncol(label) # num of clusters + 1
	adjacency=rep(0,(K-1)^3)
    # cat("length of adj:",length(adjacency),'\n')
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


noc=ncol(label)-1 # num of clusters
labelSum=apply(label[,1:noc],2,sum)
net_vec<-shrink(triplet[,1:3],label)
fullNet_vec<-fullyConnect(labelSum)
net<-normalize(noc,net_vec,fullNet_vec)
r=sum(as.numeric(net_vec))/sum(as.numeric(fullNet_vec))
net=net/r
net=array(net,c(noc,noc,noc))
save(net_vec,fullNet_vec,net,file='net_hclus.bin')


