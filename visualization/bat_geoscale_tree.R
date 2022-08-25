library("strap")

setwd("C:/Users/Nikki/OneDrive - Texas Tech University")

bat.ages <- as.matrix(read.table("bat_fad_lad.tsv", sep = "\t"))
bat.nex <- read.nexus("bat_ht_tree_timed.trees.nex")
bat.tree.test <- plot.phylo(bat.nex)
#bat.tree <- DatePhylo(bat.nex, bat.ages, method="equal", rlen=19)
bat.tree <- DatePhylo(bat.nex, bat.ages, method="ruta", rlen=19)

geoscalePhylo(ladderize(bat.tree), bat.ages, units = c("Period","Epoch"),
              boxes = "Epoch", cex.age = 0.7, cex.ts = 0.7, cex.tip = 0.8, width = 0.5, label.offset = 0.1)

