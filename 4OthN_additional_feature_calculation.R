source("./0_presets.R")
load("./All_Merged, Imputated, Calculated Features.RData")

OthAb_normal <- read.csv("../../data/genelist_Normal.csv") %>% .[!.$Phenotype.Other == "normal",] %>% .[,c(1:3)]
sn_data <- sn_data[(!sn_data$GeneName %in% OthAb_normal$GeneName) | (!sn_data$Ensembl.ID %in% OthAb_normal$Ensembl.ID),]

dim(sn_data)
table(sn_data$Phenotype.Male)

save(sn_data,file = "OthN_Merged, Imputated, Calculated Features.RData")
