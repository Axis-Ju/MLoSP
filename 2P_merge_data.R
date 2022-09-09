source("./0_presets.R")
load("./All_Unmerged, Unimputated, Uncalculated Features.RData")

# 00.Genelists to Predict -----
pooled_genes.GeneName <- Reduce(
  union,
  list(
    mouse.data.GeneName = bitr(mouse.data$Ensembl.ID,fromType = 'ENSEMBL',toType = c('SYMBOL'),OrgDb = 'org.Mm.eg.db') %>% .[,2],
    human.data$GeneName,gtex.human.data$GeneName,
    gse107644$GeneName,gse112393$GeneName,
    gse136220.GeneName = bitr(gse136220$Ensembl.ID,fromType = 'ENSEMBL',toType = c('SYMBOL'),OrgDb = 'org.Mm.eg.db') %>% .[,2],
    gse148032$GeneName,gse97519$GeneName,gse108974$GeneName,gse108970$GeneName,gse82174$GeneName,gse107711$GeneName,
    gse106487$GeneName,lab.SCSeq$GeneName,gse134144$GeneName,gse124263$GeneName,gse112013$GeneName,gse142585$GeneName,gse86146$GeneName,gse108977$GeneName)) %>% na.omit()
to_pred <- data.frame(
  GeneName = pooled_genes.GeneName,
  Phenotype.Male = rep("Others",length(pooled_genes.GeneName)));dim(to_pred)
# Translate ID
pooled_genes.Ensembl <- bitr(pooled_genes.GeneName,fromType = 'SYMBOL',toType = c('ENSEMBL'),OrgDb = 'org.Mm.eg.db',drop = TRUE)
names(pooled_genes.Ensembl) <- c("GeneName","Ensembl.ID")
to_pred <- merge(to_pred,pooled_genes.Ensembl,all.X = TRUE) %>% 
  .[,c(1,3,2)]
# Remove Known
to_pred <- to_pred[!to_pred$GeneName %in% gl$GeneName,];dim(to_pred)
to_pred <- to_pred[!to_pred$Ensembl.ID %in% gl$Ensembl.ID,];dim(to_pred)
# Remove Non-Coding
pro_coding_list <- read.table("../../data/gencode.vM15.annotation.txt") %>% .[.$gene_type == "protein_coding",15] %>% .[!duplicated(.)]
to_pred <- to_pred[to_pred$GeneName %in% pro_coding_list,];dim(to_pred)

# 01.Merge with Genomic Features -----
### 1.1 RVIS Score
pred_data <- merge(to_pred,RVIS,all.x = TRUE)
### 1.2 Shet Score
pred_data <- merge(pred_data,Shet,all.x = TRUE)
### 1.3 Constraint
pred_data <- merge(pred_data,constraint,all.x = TRUE)

# 02.Merge with mRNA Expression Features -----
### 2.1 ZMBH Mouse Testis
pred_data <- merge(pred_data,mouse.data,all.x = TRUE)
### 2.2 ZMBH Human Testis
pred_data <- merge(pred_data,human.data,all.x = TRUE)
### 2.3 GTEx Human Testis
pred_data <- merge(pred_data,gtex.human.data,all.x = TRUE)

# 03.Merge with SCS of Spermatogenesis -----
### 3.1 (Mus Spermatogenesis) GSE107644
pred_data <- merge(pred_data,gse107644,all.x = TRUE)
### 3.2 (Mus Somat) GSE112393
pred_data <- merge(pred_data,gse112393,all.x = TRUE)
### 3.3 (Mus Gonad) GSE136220
pred_data <- merge(pred_data,gse136220,all.x = TRUE)
### 3.4 (Mus Gonad) GSE148032
pred_data <- merge(pred_data,gse148032,all.x = TRUE)
### 3.5 (Mus E10.5-E16.5 Gonad) GSE97519
pred_data <- merge(pred_data,gse97519,all.x = TRUE)
### 3.6 (Mus Adult SPG) GSE108974
pred_data <- merge(pred_data,gse108974,all.x = TRUE)
### 3.7 (Mus P6 SPG) GSE108970
pred_data <- merge(pred_data,gse108970,all.x = TRUE)
### 3.8 (Mus P3&7 SPG) GSE82174
pred_data <- merge(pred_data,gse82174,all.x = TRUE)
### 3.9 (Mus P5.5 SCS) GSE107711
pred_data <- merge(pred_data,gse107711,all.x = TRUE)
### 3.10 (Hom Spermatogenesis) GSE106487
pred_data <- merge(pred_data,gse106487,all.x = TRUE)
### 3.11 (Hom Spermatogenesis) Lab
pred_data <- merge(pred_data,lab.SCSeq,all.x = TRUE)
### 3.12 (Hom Puberty) GSE134144
pred_data <- merge(pred_data,gse134144,all.x = TRUE)
### 3.13 (Hom Somat) GSE124263
pred_data <- merge(pred_data,gse124263,all.x = TRUE)
### 3.14 (Hom Somat) GSE112013
pred_data <- merge(pred_data,gse112013,all.x = TRUE)
### 3.15 (Hom Somat) GSE142585
pred_data <- merge(pred_data,gse142585,all.x = TRUE)
### 3.16 (Hom Gonad) GSE86146
pred_data <- merge(pred_data,gse86146,all.x = TRUE)
### 3.17 (Hom Adult SPG) GSE108977
pred_data <- merge(pred_data,gse108977,all.x = TRUE)

# Display
dim(pred_data)

# 04.Save -----
save(pred_data,file = "Merged, Unimputated, Uncalculated Features (Prediction).RData")
