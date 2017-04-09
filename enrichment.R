# setwd("C:\\Users\\DELL\\Dropbox\\Research_Yu\\darkness")
# setwd("C:\\Users\\kyccw\\Dropbox\\Research_Yu\\darkness")

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
                        ,testDirection="over", annotation="org.Sc.sgd.db")
        over.pres<-hyperGTest(params)
        res[[i]] <-summary(over.pres)
        cat(i,"out of",k,'done.\n')
    }
    return(res)
}

load('label_hclus.bin')
load("spellman_73_filled.bin")
glist=rownames(array)
enrichment=GOen(glist,label)
save(enrichment,file='enrichment_hclus.bin')
