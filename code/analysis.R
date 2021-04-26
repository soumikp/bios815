library(matrixStats)
library(vegan)
library(GUniFrac)
library(MedTest)
require(dirmult)

load("~/data/aofib.rda")
load("~/data/bmi.rda")
load("~/data/otu.tab.rda")
load("~/data/tree.rooted.rda")

unifracs <- GUniFrac(otu.tab, tree.rooted)$unifracs
m.list <- list(BC=vegdist(otu.tab, method="bray"),
               JAC=as.matrix(vegdist(otu.tab, 'jaccard', binary=TRUE)),
               UniFrac=unifracs[, , c('d_UW')],
               GUniFrac=unifracs[, , c('d_0.5')],
               WUniFrac=unifracs[, , c('d_1')])

set.seed(12345)
lapply(permutation_test(x= aofib, y = bmi,m.list), round, 3)
