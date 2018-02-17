# setwd("C:\\Users\\DELL\\Dropbox\\Research_Yu\\fortress")
# setwd("C:\\Users\\kyccw\\Dropbox\\Research_Yu\\fortress")

# dyn.load('csupp.so')
# setwd("/home/ykong24/GMG/ykong24")

setwd("C:\\Users\\yunchuan\\Dropbox\\Research_Yu\\PROJECT1\\fortress")
dyn.load("csupp.dll")
setwd("C:\\Users\\yunchuan\\Dropbox\\Research_Yu\\PROJECT1\\goblin")
library(GO.db)

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

load("triplet_ski.bin")
load("net_ski_GO.bin")
load("label_ski_GO.bin")
load("data_ski.bin")
load("GO_select_skin_ 0.7 0.4 1000 100 .bin")
# load("triplet_yea.bin")
# load("net_yea_GO.bin")
# load("label_yea_GO.bin")
# load("data_yea.bin")
# load("GO_select_yeast_ 0.7 0.4 1000 20 .bin")
rownames(dat) <- glist
setwd("C:\\Users\\yunchuan\\Dropbox\\Research_Yu\\PROJECT1\\Revision1_1")

module_names <- NULL
GO <- names(GO.select)
for (i in 1:length(GO.select)){
  module_names <- c(module_names, Term(GO[i]))
}

# which(module_names == "DNA damage response, signal transduction by p53 class mediator")
# # 382
# which(module_names == "sphingolipid metabolic process")
# # 130
# which(module_names == "glycolipid metabolic process")
# # 274
# 
# which(module_names == "DNA???dependent DNA replication")
# # 358
# which(module_names == "visual perception")
# # 388
# which(module_names == "coenzyme biosynthetic process")
# # 268
# which(module_names == "vacuole organization")
# # 154
# which(module_names == "Golgi vesicle transport")
# # 343
# which(module_names == "membrane lipid biosynthetic process")
# # 136
# 
h_list <- c(382, 382, 130,
            382, 382, 274,
            358, 388, 268,
            358, 388, 154,
            358, 388, 343,
            358, 388, 136)

# which(module_names == "sulfur compound transport")
# # 42
# which(module_names == "cell redox homeostasis")
# # 34
# which(module_names == "Golgi to vacuole transport")
# # 176
# which(module_names == "plasma membrane organization")
# # 44
# which(module_names == "amino sugar biosynthetic process")
# # 137
# 
# h_list <- c(42, 42, 34,
#             42, 176, 44,
#             42, 42, 137)

h_list <- matrix(h_list, ncol=3, byrow=T)
# write.table(h_list, "h_list_ski.txt", col.names=F, row.names=F)

f = 2 ## corresponding fold changes
net_sel <- net
net_sel[net<f] = 0

## inquiries from a hyperedge list
source("C:\\Users\\yunchuan\\Dropbox\\Research_Yu\\PROJECT1\\Revision1_1\\plot examples_Yu_modified.r")
# h_list <- as.matrix(read.table("h_list_ski.txt", header=F))
for (i in 1:nrow(h_list)){
  hyperedge <- h_list[i,]
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

  sel_list <- query(triplet, c1, c2, c3)
  sel_list <- sel_list[sel_list!=0]
  sel_list <- matrix(sel_list, ncol=3, byrow=T)
  
  ## sort triplets by LA decending order
  la <- rep(0, nrow(sel_list))
  calculate_LA <- function(x,y,z){sum(x*y*z)/length(x)}
  for (j in 1:nrow(sel_list)){
    la[j] <- calculate_LA(dat[sel_list[j,1],], dat[sel_list[j,2],], dat[sel_list[j,3],])
  }
  if (nrow(sel_list)==1){
    sel_list <- t(as.matrix(sel_list[order(la, decreasing=T),],ncol=3,nrow=1))
  }
  else{
    sel_list <- as.matrix(sel_list[order(la, decreasing=T),],ncol=3,nrow=1)
  }
  
  file <- paste0("hyperedge_yea_", i,".pdf")
  example.pdf(file, dat, sel_list)
  # write.table(sel_list, file=file, col.names=F, row.names=F)
  print(paste0(file," finished."))
}

#################################################################################################################
## enquiry for a single hyperedge
# hyperedge <- c(353,  291,   51)
# hyperedge <- sort(hyperedge, decreasing=T)
# if (net_sel[hyperedge[1],hyperedge[2],hyperedge[3]] == 0){
#   cat("The hyperedge does not exist under the current fold change.\n")
# }
# 
# c1 <- label[,hyperedge[1]]
# c1 <- which(c1==1)
# c2 <- label[,hyperedge[2]]
# c2 <- which(c2==1)
# c3 <- label[,hyperedge[3]]
# c3 <- which(c3==1)
# 
# sel_list <- query(triplet, c1, c2, c3)
# sel_list <- sel_list[sel_list!=0]
# sel_list <- matrix(sel_list, ncol=3, byrow=T)
# print(sel_list)
# write.table(sel_list, file=paste0(hyperedge,".txt"), col.names=F, row.names=F)
