##library
library(brms)
library(ape)

##remove prior data
rm(list=ls()) 

##getting rid of divergent transitions
##get prior results
load("multiresponse_TE_cumu.RData")

##the time has come to stop being fancy and regress like we know how to
##assume TE values known without error

##increase iterations, raise adapt_delta and treedepth and rerun indivdidual models to eliminate or reduce divergent transitions 

##species proportion to helitron
pro.r <- bf(pro ~ Helitron + (1|gr(species, cov = A)))

##number Helitron
hel <- bf(Helitron ~ (1|gr(species, cov = A)))

##run all models together at a time
m.hel <- brm(pro.r, 
            family = "beta", 
                  data2 = list(A = A),
                  iter  = 100000, warmup = 50000, thin = 50, 
                  control = list(adapt_delta = 0.99999),
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
                  iter  = 100000, warmup = 50000, thin = 50, 
                  control = list(adapt_delta = 0.99999),
                  data = bats,
                  cores = 4, chains = 4)

##save
save.image("multiresponse_TE_cumu.RData")
