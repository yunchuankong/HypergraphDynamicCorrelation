## modified by Yunchuan

example.pdf<-function(filename, new.array, xyz, triplets.to.plot=20, cex=0.5)
{
    library(mclust)
    library(locfdr)
    
    pdf(filename, width=8, height=8)

    
    # for(n in 1:min(nrow(xyz),triplets.to.plot))
	for(n in 1:nrow(xyz))
    {
				par(mfrow=c(2,2))
                r<-plot.la(x=new.array[xyz[n,1],], y=new.array[xyz[n,2],], z=new.array[xyz[n,3],], xlab=rownames(new.array)[xyz[n,1]], ylab=rownames(new.array)[xyz[n,2]], use.mclust=T, num.grp=2, cex=cex)
				plot(0,0,type="n", axes=F, xlab="", ylab="", ylim=c(0,1), xlim=c(0,1))
                text(0, 0.9, paste(rownames(new.array)[xyz[n,1]], rownames(new.array)[xyz[n,2]]), pos=4)
                text(0, 0.8, "locfdr grouping", pos=4)
                text(0, 0.7, paste("z:", rownames(new.array)[xyz[n,3]]), pos=4)
                for(i in 1:length(r)) text(0, 0.6-0.1*i, r[i], pos=4)


                 r<-plot.la(x=new.array[xyz[n,1],], y=new.array[xyz[n,2],], z=new.array[xyz[n,3],], xlab=rownames(new.array)[xyz[n,1]], ylab=rownames(new.array)[xyz[n,2]], use.mclust=F, cex=cex)
                 plot(0,0,type="n", axes=F, xlab="", ylab="", ylim=c(0,1), xlim=c(0,1))
                 text(0, 0.9, paste(rownames(new.array)[xyz[n,1]], rownames(new.array)[xyz[n,2]]), pos=4)
                 text(0, 0.8, "z quantiles", pos=4)
                 text(0, 0.7, paste("z:", rownames(new.array)[xyz[n,3]]), pos=4)
                 for(i in 1:length(r)) text(0, 0.6-0.1*i, r[i], pos=4)
    }

    dev.off()
}



addTrans<-function(color,trans)
{
    # This function adds transparancy to a color.
    # Define transparancy with an integer between 0 and 255
    # 0 being fully transparant and 255 being fully visable
    # Works with either color and trans a vector of equal length,
    # or one of the two of length 1.
    
    if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
    if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
    if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
    
    num2hex <- function(x)
    {
        hex <- unlist(strsplit("0123456789ABCDEF",split=""))
        return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
    }
    rgb <- rbind(col2rgb(color),trans)
    res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
    return(res)
}

plot.la<-function(x,y,z, use.mclust=FALSE,xlab="x",ylab="y", cols=c("red","green","blue","purple","cyan"),num.grp=NULL, cex=0.5)
{
    
    if(use.mclust)
    {
        #l<-Mclust(z, modelNames ="V", G=num.grp, prior = priorControl(functionName="defaultPrior"))
        l<-locfdr(z, plot=0)
        grp<-rep(2, length(z))
        grp[z<=max(z[z<median(z) & l$fdr<=0.5])]<-1
        grp[z>=min(z[z>median(z) & l$fdr<=0.5])]<-3
    }else{
        cuts<-quantile(z, c(0.3333, 0.6666))
        grp<-rep(2, length(z))
        grp[z<cuts[1]]<-1
        grp[z>cuts[2]]<-3
    }
    
    plot(x,y, type="n", xlab=xlab, ylab=ylab)
    
    all.grp<-unique(grp)
    all.grp<-all.grp[order(all.grp)]
    
    r<-NULL

#   message(c("LA score: ", sum(x*y*z)))
    r<-c(r,paste("LA score: ", signif(sum(x*y*z),3)))
    
    for(i in 1:length(all.grp))
    {
        s<-which(grp == all.grp[i])
        points(x[s],y[s],col=addTrans(cols[i],120), pch=19, cex=cex)
        
        
        #message(c("group: ", i,"  ", cols[i], "  corr=", round(cor(x[s],y[s]), 3)))
        r<-c(r, paste("group:", i,"(", signif(median(z[s]),2), ")", cols[i], ", corr=", round(cor(x[s],y[s]), 3)))
    }
    r
}


plot.la.fixed<-function(x,y,grp,cols=c("green","red","blue","purple","cyan"))
# the grouping is given
{
    plot(x,y)
    
    all.grp<-unique(grp)
    
    #message(c("LA score: ", sum(x*y*z)))
    
    for(i in 1:length(all.grp))
    {
        s<-which(grp == all.grp[i])
        points(x[s],y[s],col=cols[i], pch=19)
        
        message(c("group: ", i,"  ", cols[i], "  corr=", round(cor(x[s],y[s]), 3)))
    }
}
