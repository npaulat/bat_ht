##script to count descendants assigned to each tip

##libraries
library(phytools)
library(geiger)

##bigtree
tbig<-read.nexus("T87680_mod.nex")

##small tree
tree<-read.nexus("bat_ht_tree_timed.trees.nex")

##cleaned up names to match
##check for mismatches in names, returns empty now
setdiff(tree$tip,tbig$tip)

##get matching names in big tree
bN<-as.numeric(sapply(tree$tip.label,function(x) which(tbig$tip.label %in% x)))

##nodes
nodes<-1:tbig$Nnode+Ntip(tbig) ## all nodes

##all subtrees
subtrees<-list()

#extract all clades
for(i in 1:tbig$Nnode) subtrees[[i]]<-extract.clade(tbig,nodes[i])

##names these subtrees
names(subtrees)<-nodes ## all subtrees

##get all the subtrees that have target tips
len.sub<-sapply(subtrees,function(x) length(which(tree$tip.label %in% x$tip.label  )))

##make dataframe of subtree content of target taxa
subt<-as.data.frame(cbind( nodes, len.sub))

##make sure only subtrees with one are kept
subt<-subt[subt$len.sub ==1,]

##descendants of each node on the big tree
desc<-sapply(tbig$edge[,2],function(x) tips(tbig,x))

##now count
Ndesc<-sapply(desc,function(x) length(x))

##add root descendents internal value
Ndesc <- append(815, Ndesc)

##make data frame of counts and nodes they came from
dat<-as.data.frame(cbind(Ndesc, c(tbig$edge[1,1], tbig$edge[,2])))

##merge
subt<-merge(subt,dat,  by.y="V2", by.x="nodes")

pdf("bigtree_v2.pdf", width=36, height=36)
plotTree(tbig,type="fan", ftype="off",lwd=1)
##use same order that they came from
labelnodes(subt$Ndesc, node=subt$nodes, interactive=F, cex=2)
tiplabels(tbig$tip.label[bN], bN, adj = 0, cex=1.5)
dev.off()
