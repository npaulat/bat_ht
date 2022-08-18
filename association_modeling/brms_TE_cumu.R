##library
library(brms)
library(ape)

##clear out cache
rm(list=ls())

##get the data
bats<-read.csv("bats_accumulation.csv")

##get the tree
tree<-read.nexus("bat_ht_tree_timed.trees.nex")

##get phylogenetic matrix
A <- ape::vcv.phylo(tree)

##make a proportion variable
bats$pro<-bats$Nspp/sum(bats$Nspp)

##count values are too high, messing up the beta regression
##predictors need to be scaled
##they no longer will be counts though
##save the original data
bat<-bats

##scale 
bats[,c(2:12)]<-log10(bats[,c(2:12)])
bats$Helitron[is.infinite(bats$Helitron)]<-0
bats[,c(2:12)]<-scale(bats[,c(2:12)])

##the time has comew to stop being fancy and regress like we know how to
##assume TE values known without error

##increase iterations, raise adapt_delta and treedepth and rerun indivdidual models to eliminate or reduce divergent transitions 

##species proportion to total DNA
pro.t <- bf(pro ~ Total.DNA.transposons + (1|gr(species, cov = A)))

##run all models together at a time
m.tot <- brm(pro.t, 
            family = "beta", 
                  data2 = list(A = A),
                  iter  = 50000, warmup = 20000, thin = 100, 
                  control = list(adapt_delta = 0.9999, max_treedepth=20),
                  data = bats,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_TE_cumu.RData")

##species proportion to total retro
pro.b <- bf(pro ~ Total.Retrotransposons + (1|gr(species, cov = A)))

##number total TEs
ttt <- bf(Total.Retrotransposons ~ (1|gr(species, cov = A)))

##run all models together at a time
m.ttt <- brm(pro.b, 
            family = "beta", 
                  data2 = list(A = A),
                  iter  = 20000, warmup = 4000, thin = 40, 
                  control = list(adapt_delta = 0.9999, max_treedepth=20),
                  data = bats,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_TE_cumu.RData")

##species proportion to hAT
pro.h <- bf(pro ~ hAT + (1|gr(species, cov = A)))

##number hAT
hat <- bf(hAT ~ (1|gr(species, cov = A)))

##run all models together at a time
m.hat <- brm(pro.h, 
            family = "beta", 
                  data2 = list(A = A),
                  iter  = 50000, warmup = 20000, thin = 100, 
                  control = list(adapt_delta = 0.9999, max_treedepth=20),
                  data = bats,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_TE_cumu.RData")

##species proportion to piggy
pro.p <- bf(pro ~ piggyBac + (1|gr(species, cov = A)))

##number piggyBac
pig <- bf(piggyBac ~ (1|gr(species, cov = A)))

##run all models together at a time
m.pig <- brm(pro.p, 
            family = "beta", 
                  data2 = list(A = A),
                  iter  = 20000, warmup = 4000, thin = 40, 
                  control = list(adapt_delta = 0.9999),
                  data = bats,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_TE_cumu.RData")

##species proportion to mariner
pro.m <- bf(pro ~ TcMariner + (1|gr(species, cov = A)))

##number TcMariner
mar <- bf(TcMariner ~ (1|gr(species, cov = A)))

##run all models together at a time
m.mar <- brm(pro.m, 
            family = "beta", 
                  data2 = list(A = A),
                  iter  = 20000, warmup = 4000, thin = 40, 
                  control = list(adapt_delta = 0.999),
                  data = bats,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_TE_cumu.RData")

##species proportion to dna
pro.d <- bf(pro ~ DNA + (1|gr(species, cov = A)))

##number DNA
dna <- bf(DNA ~ (1|gr(species, cov = A)))

##run all models together at a time
m.dna <- brm(pro.d, 
            family = "beta", 
                  data2 = list(A = A),
                  iter  = 50000, warmup = 20000, thin = 100, 
                  control = list(adapt_delta = 0.999),
                  data = bats,
                  cores = 3, chains = 3)

##save
save.image("multiresponse_TE_cumu.RData")

##species proportion to helitron
pro.r <- bf(pro ~ Helitron + (1|gr(species, cov = A)))

##number Helitron
hel <- bf(Helitron ~ (1|gr(species, cov = A)))

##run all models together at a time
m.hel <- brm(pro.r, 
            family = "beta", 
                  data2 = list(A = A),
                  iter  = 20000, warmup = 4000, thin = 40, 
                  control = list(adapt_delta = 0.999),
                  data = bats,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_TE_cumu.RData")

##species proportion to line
pro.l <- bf(pro ~ LINE + (1|gr(species, cov = A)))

##number line
lin <- bf(LINE ~ (1|gr(species, cov = A)))

##run all models together at a time
m.lin <- brm(pro.l, 
            family = "beta", 
                  data2 = list(A = A),
                  iter  = 50000, warmup = 20000, thin = 100, 
                  control = list(adapt_delta = 0.999),
                  data = bats,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_TE_cumu.RData")

##species proportion to sine
pro.s <- bf(pro ~ SINE + (1|gr(species, cov = A)))

##number sine
sin <- bf(SINE ~ (1|gr(species, cov = A)))

##run all models together at a time
m.sin <- brm(pro.s, 
            family = "beta", 
                  data2 = list(A = A),
                  iter  = 20000, warmup = 4000, thin = 40, 
                  control = list(adapt_delta = 0.999),
                  data = bats,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_TE_cumu.RData")

##species proportion to ltr
pro.a <- bf(pro ~ LTR + (1|gr(species, cov = A)))

##number ltr
ltr <- bf(LTR ~ (1|gr(species, cov = A)))

##run all models together at a time
m.ltr <- brm(pro.a, 
            family = "beta", 
                  data2 = list(A = A),
                  iter  = 20000, warmup = 4000, thin = 40, 
                  control = list(adapt_delta = 0.999),
                  data = bats,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_TE_cumu.RData")
