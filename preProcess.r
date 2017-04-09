# setwd("C:\\Users\\DELL\\Dropbox\\Research_Yu\\elves")
# setwd("C:\\Users\\kyccw\\Dropbox\\Research_Yu\\elves")

# library(fdrtool)
library(coin)
load("BRCA.tumor.bin") # 20532 by 762

manyZeros<-function(x){
    zeroRate=0
    if (length(which(x==0))/length(x)>zeroRate)
        return(1)
    else
        return(0)
}

arctanh<-function(x){
	return( 1/2*log((1+x)/(1-x)) )
}

# logNorm<-function(x){
    # lnX=log(x)
    # return( (lnX-mean(lnX))/sd(lnX) )
# }

z=apply(array,1,manyZeros)
zeroRows=which(z==1)
dat=array[-zeroRows,] # 14760 by 762
# dat=dat[1:333,]
dat=apply(dat,2,normal_trafo) # normal_trafo is preferred over logNorm
# dat=apply(dat,2,logNorm)
nog=nrow(dat)
nos=ncol(dat)
glist=rownames(dat)
clist=colnames(dat)
rownames(dat)=1:nog
colnames(dat)=NULL

dict=cor(t(dat))


# v=dict[lower.tri(dict,diag=F)]
# for (i in 1:15)
	# cat(i,':',quantile(v,(50-i)/100),quantile(v,(50+i)/100),'\n')
# quantile 45% - 55%: -0.0225599 0.01799315

# dict_z=arctanh(dict)
# pdf('histograms.pdf')
# hist(v)
# hist(dict_z[lower.tri(dict_z,diag=F)])
# dev.off()



# threshold_cor = 0.1
# dict=ifelse(abs(dict)<threshold_cor, 1, 0)
diag(dict)=0
dict=ifelse((dict>-0.023 & dict<0.018) , 1, 0)
save(dat
    ,clist
    ,dict
    ,glist
    ,nog
    ,nos
    ,file="data.bin")

                 


