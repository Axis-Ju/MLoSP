load("./All_Unmerged, Unimputated, Uncalculated Features.RData")

# 01.Merge with Genomic Features -----
### 1.1 RVIS Score
sn_data <- merge(gl,RVIS,all.x = TRUE)
### 1.2 Shet Score
sn_data <- merge(sn_data,Shet,all.x = TRUE)
### 1.3 Constraint
sn_data <- merge(sn_data,constraint,all.x = TRUE)

# 02.Merge with mRNA Expression Features -----
### 2.1 ZMBH Mouse Testis
sn_data <- merge(sn_data,mouse.data,all.x = TRUE)
### 2.2 ZMBH Human Testis
sn_data <- merge(sn_data,human.data,all.x = TRUE)
### 2.3 GTEx Human Testis
sn_data <- merge(sn_data,gtex.human.data,all.x = TRUE)

# 03.Merge with SCS of Spermatogenesis -----
### 3.1 (Mus Spermatogenesis) GSE107644
sn_data <- merge(sn_data,gse107644,all.x = TRUE)
### 3.2 (Mus Somat) GSE112393
sn_data <- merge(sn_data,gse112393,all.x = TRUE)
### 3.3 (Mus Gonad) GSE136220
sn_data <- merge(sn_data,gse136220,all.x = TRUE)
### 3.4 (Mus Gonad) GSE148032
sn_data <- merge(sn_data,gse148032,all.x = TRUE)
### 3.5 (Mus E10.5-E16.5 Gonad) GSE97519
sn_data <- merge(sn_data,gse97519,all.x = TRUE)
### 3.6 (Mus Adult SPG) GSE108974
sn_data <- merge(sn_data,gse108974,all.x = TRUE)
### 3.7 (Mus P6 SPG) GSE108970
sn_data <- merge(sn_data,gse108970,all.x = TRUE)
### 3.8 (Mus P3&7 SPG) GSE82174
sn_data <- merge(sn_data,gse82174,all.x = TRUE)
### 3.9 (Mus P5.5 SCS) GSE107711
sn_data <- merge(sn_data,gse107711,all.x = TRUE)
### 3.10 (Hom Spermatogenesis) GSE106487
sn_data <- merge(sn_data,gse106487,all.x = TRUE)
### 3.11 (Hom Spermatogenesis) Lab
sn_data <- merge(sn_data,lab.SCSeq,all.x = TRUE)
### 3.12 (Hom Puberty) GSE134144
sn_data <- merge(sn_data,gse134144,all.x = TRUE)
### 3.13 (Hom Somat) GSE124263
sn_data <- merge(sn_data,gse124263,all.x = TRUE)
### 3.14 (Hom Somat) GSE112013
sn_data <- merge(sn_data,gse112013,all.x = TRUE)
### 3.15 (Hom Somat) GSE142585
sn_data <- merge(sn_data,gse142585,all.x = TRUE)
### 3.16 (Hom Gonad) GSE86146
sn_data <- merge(sn_data,gse86146,all.x = TRUE)
### 3.17 (Hom Adult SPG) GSE108977
sn_data <- merge(sn_data,gse108977,all.x = TRUE)

# Display
dim(sn_data)
table(sn_data$Phenotype.Male)

# 04.Save -----
save(sn_data,file = "All_Merged, Unimputated, Uncalculated Features.RData")
