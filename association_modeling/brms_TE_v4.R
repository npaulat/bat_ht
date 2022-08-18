##library
library(brms)
library(ape)

##clear out cache
rm(list=ls())

##get the data
bats<-read.csv("bats_PA_tables_edited_v3.csv")

##get the tree
tree<-read.nexus("bat_ht_tree_timed.trees.nex")

##get phylogenetic matrix
A <- ape::vcv.phylo(tree)

##make a proportion variable
bats$pro<-bats$Nspp/sum(bats$Nspp)

##increase iterations, raise adapt_delta and treedepth and rerun indivdidual models to eliminate or reduce divergent transitions 

##species proportion to total
pro.t <- bf(pro ~ Total + (1|gr(species, cov = A)))

##number total TEs
tot <- bf(Total ~ (1|gr(species, cov = A)))

##run all models together at a time
m.tot <- brm(pro.t +
        	tot +
            set_rescor(FALSE), 
            family = list("beta",
            negbinomial(link = "log")), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 20, 
                  control = list(adapt_delta = 0.9999, max_treedepth=20),
                  data = bats,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_TE_v4.RData")

##species proportion to hAT
pro.h <- bf(pro ~ hAT + (1|gr(species, cov = A)))

##number hAT
hat <- bf(hAT ~ (1|gr(species, cov = A)))

##run all models together at a time
m.hat <- brm(pro.h +
        	hat +
            set_rescor(FALSE), 
            family = list("beta",
            negbinomial(link = "log")), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 20, 
                  control = list(adapt_delta = 0.999, max_treedepth=15),
                  data = bats,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_TE_v4.RData")

##species proportion to piggy
pro.p <- bf(pro ~ piggyBac + (1|gr(species, cov = A)))

##number piggyBac
pig <- bf(piggyBac ~ (1|gr(species, cov = A)))

##run all models together at a time
m.pig <- brm(pro.p +
        	pig +
            set_rescor(FALSE), 
            family = list("beta",
            negbinomial(link = "log")), 
                  data2 = list(A = A),
                  iter  = 20000, warmup = 4000, thin = 40, 
                  control = list(adapt_delta = 0.9999),
                  data = bats,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_TE_v4.RData")

##species proportion to mariner
pro.m <- bf(pro ~ TcMariner + (1|gr(species, cov = A)))

##number TcMariner
mar <- bf(TcMariner ~ (1|gr(species, cov = A)))

##run all models together at a time
m.mar <- brm(pro.m +
        	mar +
            set_rescor(FALSE), 
            family = list("beta",
            negbinomial(link = "log")), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 4000, thin = 20, 
                  control = list(adapt_delta = 0.9999),
                  data = bats,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_TE_v4.RData")

##species proportion to helitron
pro.r <- bf(pro ~ Helitron + (1|gr(species, cov = A)))

##number Helitron
hel <- bf(Helitron ~ (1|gr(species, cov = A)))

##run all models together at a time
m.hel <- brm(pro.r +
        	hel +
            set_rescor(FALSE), 
            family = list("beta",
            negbinomial(link = "log")), 
                  data2 = list(A = A),
                  iter  = 10000, warmup = 2000, thin = 20, 
                  control = list(adapt_delta = 0.999),
                  data = bats,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_TE_v4.RData")
