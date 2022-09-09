# Library Packages -----
library(gbm)
library(mlr)
library(mice)
library(ROCR)
library(pROC)
library(DMwR)
library(ROSE)
library(caret)
library(Mfuzz)
library(ggpubr)
library(writexl)
library(ggplot2)
library(biomaRt)
library(magrittr)
library(gg.layers)
library(homologene)
library(VennDiagram)
library(philentropy)
library(randomForest)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(caretEnsemble)
library(clusterProfiler)

rm(list = ls(all.names = TRUE))
gc()

# Prepare Functions -----
modulus <- function(v) {
  sqrt(sum(v^2))
}
outHigh <- function(x) {
  x[x > quantile(x, 0.95)] <- quantile(x, 0.75)
  x
}
outLow <- function(x) {
  x[x < quantile(x, 0.05)] <- quantile(x, 0.25)
  x
}

# Load Public Data -----
load("../../data/MultiSpecies.Testis.RData")
load("../../data/MultiSpecies.RData")
load('../../data/GTEx.Testis.RData')
load("../../data/GTEx.All.RData")
load("../../data/SCSeq.Testis.RData")