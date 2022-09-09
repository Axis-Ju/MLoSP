source("./0_presets.R")
load("./All_Merged, Imputated, Uncalculated Features.RData")





#####
##### mRNA Expression Features -------------------------------------------------
#####





# 1.1 ZMBH Mouse Testis -----
data <- sn_data[,c(2,grep("ZMBH_Mouse",names(sn_data)))]
data_all <- mouse
# Timepoint Mean
time_list <- c("e10.5","e11.5","e12.5","e13.5","e14.5","e15.5","e16.5","e17.5","e18.5","P0","P3","P14","P28","P63")
timepoint_mean <- as.data.frame(array(dim = c(nrow(data),0)))
for (time in time_list) {
  time.idx <- grep(time,names(data)) ### Which timepoint to calculate mean?
  ifelse(length(time.idx) > 1,tmp_mean <- rowMeans(data[,time.idx]),tmp_mean <- data[,time.idx]) ### In case some timepoint e.g. e18.5 has only one sample
  timepoint_mean <- cbind(timepoint_mean,tmp_mean)
  names(timepoint_mean)[ncol(timepoint_mean)] <- time
}
# Alltime Mean
alltime_mean <- rowMeans(data[,2:ncol(data)])
# Prenatal Mean
prenatal_mean <- rowMeans(data[,2:(grep("P0",names(data))[1]-1)])
# Postnatal Mean
postnatal_mean <- rowMeans(data[,(grep("P0",names(data))[1]):ncol(data)])
# Temporal SPM
time_list <- c("e10.5","e11.5","e12.5","e13.5","e14.5","e15.5","e16.5","e17.5","e18.5","P0","P3","P14","P28","P63")
temporal_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (time in time_list) {
  time.idx <- grep(time,names(data)) ### Which timepoint to calculate SPM?
  spm <- c()
  for (i in 1:nrow(data)) {
    v1 <- as.numeric(data[i,2:ncol(data)]) ### From col 2 to the end
    v2 <- rep(0,length(v1))
    v2[time.idx-1] <- v1[time.idx-1] ### 1 col before actual data
    tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
    spm <- c(spm,tmp.spm)
  }
  spm[is.na(spm)] <- 0 ### Replace NA with 0
  temporal_spm <- cbind(temporal_spm,spm)
  names(temporal_spm)[ncol(temporal_spm)] <- paste("Temporal_SPM",time,sep = "_")
}
# Temporal Dynamics
temporal_dynamics <- c()
for (i in 1:nrow(data)){
  x <- c(1:(ncol(data)-1))
  y <- as.numeric(data[i,2:ncol(data)])
  tmp_lm <- lm(y~x)
  tmp.beta <- tmp_lm$coefficients[2]
  temporal_dynamics <- c(temporal_dynamics,tmp.beta)
}
# All Tissue SPM
tissue_list <- c("Brain","Cerebellum","Heart","Kidney","Liver","Ovary","Testis")
time_list <- c("e10.5","e11.5","e12.5","e13.5","e14.5","e15.5","e16.5","e17.5","e18.5","P0","P3","P14","P28","P63")
all_tissue_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (time in time_list) {
  idx <- c()
  for (tissue in tissue_list) {
    tissue.idx <- grep(tissue,names(data_all))
    time.idx <- grep(time,names(data_all))
    idx <- c(idx,intersect(tissue.idx,time.idx))
  }
  spm <- c()
  for (i in data$Ensembl.ID) {
    if (i %in% rownames(data_all)) {
      timepoint_data <- data_all[rownames(data_all) == i,idx]
      testis.idx <- grep("Testis",names(timepoint_data))
      v1 <- as.numeric(timepoint_data[,1:ncol(timepoint_data)]) ### From col 2 to the end
      v2 <- rep(0,length(v1))
      v2[testis.idx] <- v1[testis.idx] ### 1 col before actual data
      tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
      tmp.spm <- ifelse(sum(v1[testis.idx]) == 0, 0, tmp.spm) ### Replace NA with 0
      spm <- c(spm,tmp.spm)
    } else {
      spm <- c(spm,NA)
    }
  }
  all_tissue_spm <- cbind(all_tissue_spm,spm)
  names(all_tissue_spm)[ncol(all_tissue_spm)] <- paste("All_Tissue_SPM",time,sep = "_")
}
# Selective Tissue SPM
tissue_list <- c("Heart","Kidney","Liver","Testis")
time_list <- c("e10.5","e11.5","e12.5","e13.5","e14.5","e15.5","e16.5","e17.5","e18.5","P0","P3","P14","P28","P63")
selective_tissue_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (time in time_list) {
  idx <- c()
  for (tissue in tissue_list) {
    tissue.idx <- grep(tissue,names(data_all))
    time.idx <- grep(time,names(data_all))
    idx <- c(idx,intersect(tissue.idx,time.idx))
  }
  spm <- c()
  for (i in data$Ensembl.ID) {
    if (i %in% rownames(data_all)) {
      timepoint_data <- data_all[rownames(data_all) == i,idx]
      testis.idx <- grep("Testis",names(timepoint_data))
      v1 <- as.numeric(timepoint_data[,1:ncol(timepoint_data)]) ### From col 2 to the end
      v2 <- rep(0,length(v1))
      v2[testis.idx] <- v1[testis.idx] ### 1 col before actual data
      tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
      tmp.spm <- ifelse(sum(v1[testis.idx]) == 0, 0, tmp.spm) ### Replace NA with 0
      spm <- c(spm,tmp.spm)
    } else {
      spm <- c(spm,NA)
    }
  }
  selective_tissue_spm <- cbind(selective_tissue_spm,spm)
  names(selective_tissue_spm)[ncol(selective_tissue_spm)] <- paste("Selective_Tissue_SPM",time,sep = "_")
}
# Imputation of Tissue SPM
to_imputate_data <- cbind(
  timepoint_mean,
  Alltime_Mean = alltime_mean,
  Prenatal_Mean = prenatal_mean,
  Postnatal_Mean = postnatal_mean,
  temporal_spm,
  Beta = temporal_dynamics,
  all_tissue_spm,
  selective_tissue_spm)
### Rename
names(to_imputate_data) <- paste("ZMBH_Mouse",names(to_imputate_data),sep = "_")
to_imputate_data$Response <- sn_data$Phenotype.Male %>% 
  gsub("Normal",0,.) %>% 
  gsub("Abnormal",1,.) %>% 
  as.numeric()
### Multiple Imputation of Tissue SPM
md.pattern(to_imputate_data)
imputation_data <- mice(
  to_imputate_data,
  m = 3,
  maxit = 10,
  seed = 777,
  method = "pmm",
  remove.collinear = FALSE,
  remove.constant = FALSE)
### Visualization of Imputation
#summary(imputation_data)
#stripplot(imputation_data,col = c("grey",mdc(2)),pch = c(1,20))
### Fitting to Diagnose
#imputation_data_fit <- lm.mids(Response ~ .,data = imputation_data)
#imputation_data_pooled <- pool(imputation_data_fit)
### Diagnosis of Pooled Imputation
#res.1 <- imputation_data_pooled$glanced;res.1
#res.2 <- cbind(imputation_data_pooled$pooled,summary(imputation_data_pooled))
#pool.r.squared(imputation_data_fit)
### Completion
imputed_data <- complete(imputation_data,action = 3) %>% .[,-ncol(.)]
#write.csv(cbind(data,imputed_data),"./ZMBH_Mouse.csv")
# Merge
sn_data <- sn_data[,-grep("ZMBH_Mouse",names(sn_data))] %>% cbind(.,imputed_data)

# 1.2 ZMBH Human Testis -----
data <- sn_data[,c(1,grep("ZMBH_Human",names(sn_data)))]
data_all <- human
data_all$Ensembl.ID <- rownames(data_all)
### From Human Ensembl to Human Symbol
genoinfo1 <- bitr(data_all$Ensembl.ID,fromType = 'ENSEMBL',toType = c('SYMBOL'),OrgDb = 'org.Hs.eg.db')
names(genoinfo1) <- c("Ensembl.ID","humanGene")
### From Human Symbol to Mouse Symbol
genoinfo2 <- human2mouse(genoinfo1$humanGene)
genoinfo <- merge(genoinfo1,genoinfo2,all.x = TRUE) %>% .[,2:3]
### Merge with Symbol Annotation
data_all <- merge(data_all,genoinfo,all.x = TRUE) %>%
  .[!is.na(.$mouseGene),] %>% 
  .[,c(299,2:298)]
names(data_all)[1] <- "GeneName"
### Tackle Dplications
dp_list <- data_all[duplicated(data_all$GeneName),] %>% 
  .[!duplicated(.$GeneName),]
for (i in dp_list$GeneName) {
  tp <- data_all[data_all$GeneName == i,]
  ##### Calculate Shannon’s Entropy
  tp_H <- apply(tp[,2:ncol(tp)],1,function(x) {
    x <- as.numeric(x)
    if (sum(x) == 0) {
      H = 0
    } else {
      Prob <- x/sum(x)
      H <- H(Prob) 
    }
  })
  if (!sum(tp_H == max(tp_H)) == nrow(tp)) {
    tp_max_H <- data_all[rownames(data_all) == names(tp_H[tp_H == max(tp_H)][1]),]
  } else {
    tp_max_H <- data_all[rownames(data_all) == names(tp_H[1]),]
  }
  data_all[data_all$GeneName == i,] <- tp_max_H
}
data_all <- data_all[!duplicated(data_all$GeneName),]
# Timepoint Mean
time_list <- c("4wpc","5wpc","6wpc","7wpc","8wpc","9wpc","10wpc","11wpc","12wpc","13wpc","16wpc","18wpc","19wpc","infant","toddler","youngTeenager","oldTeenager","youngAdult","youngMidAge","olderMidAge","Senior")
timepoint_mean <- as.data.frame(array(dim = c(nrow(data),0)))
for (time in time_list) {
  time.idx <- grep(time,names(data)) ### Which timepoint to calculate mean?
  ifelse(length(time.idx) > 1,tmp_mean <- rowMeans(data[,time.idx]),tmp_mean <- data[,time.idx]) ### In case some timepoint has only one sample
  timepoint_mean <- cbind(timepoint_mean,tmp_mean)
  names(timepoint_mean)[ncol(timepoint_mean)] <- time
}
# Alltime Mean
alltime_mean <- rowMeans(data[,2:ncol(data)])
# Prenatal Mean
prenatal_mean <- rowMeans(data[,2:(grep("infant",names(data))[1]-1)])
# Postnatal Mean
postnatal_mean <- rowMeans(data[,(grep("infant",names(data))[1]):ncol(data)])
# Temporal SPM
time_list <- c("4wpc","5wpc","6wpc","7wpc","8wpc","9wpc","10wpc","11wpc","12wpc","13wpc","16wpc","18wpc","19wpc","infant","toddler","youngTeenager","oldTeenager","youngAdult","youngMidAge","olderMidAge","Senior")
temporal_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (time in time_list) {
  time.idx <- grep(time,names(data)) ### Which timepoint to calculate SPM?
  spm <- c()
  for (i in 1:nrow(data)) {
    v1 <- as.numeric(data[i,2:ncol(data)]) ### From col 2 to the end
    v2 <- rep(0,length(v1))
    v2[time.idx-1] <- v1[time.idx-1] ### 1 col before actual data
    tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
    spm <- c(spm,tmp.spm)
  }
  spm[is.na(spm)] <- 0 ### Replace NA with 0
  temporal_spm <- cbind(temporal_spm,spm)
  names(temporal_spm)[ncol(temporal_spm)] <- paste("Temporal_SPM",time,sep = "_")
}
# Temporal Dynamics
temporal_dynamics <- c()
for (i in 1:nrow(data)){
  x <- c(1:(ncol(data)-1))
  y <- as.numeric(data[i,2:ncol(data)])
  tmp_lm <- lm(y~x)
  tmp.beta <- tmp_lm$coefficients[2]
  temporal_dynamics <- c(temporal_dynamics,tmp.beta)
}
# All Tissue SPM
tissue_list <- c("Brain","Cerebellum","Heart","Kidney","Liver","Ovary","Testis")
time_list <- c("4wpc","5wpc","6wpc","7wpc","8wpc","9wpc","10wpc","11wpc","12wpc","13wpc","16wpc","18wpc","19wpc","infant","toddler","youngTeenager","oldTeenager","youngAdult","youngMidAge","olderMidAge")
all_tissue_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (time in time_list) {
  idx <- c()
  for (tissue in tissue_list) {
    tissue.idx <- grep(tissue,names(data_all))
    time.idx <- grep(time,names(data_all))
    idx <- c(idx,intersect(tissue.idx,time.idx))
  }
  spm <- c()
  for (i in data$GeneName) {
    if (i %in% data_all$GeneName) {
      timepoint_data <- data_all[data_all$GeneName == i,idx]
      testis.idx <- grep("Testis",names(timepoint_data))
      v1 <- as.numeric(timepoint_data[,1:ncol(timepoint_data)]) ### From col 2 to the end
      v2 <- rep(0,length(v1))
      v2[testis.idx] <- v1[testis.idx] ### 1 col before actual data
      tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
      tmp.spm <- ifelse(sum(v1[testis.idx]) == 0, 0, tmp.spm) ### Replace NA with 0
      spm <- c(spm,tmp.spm)
    } else {
      spm <- c(spm,NA)
    }
  }
  all_tissue_spm <- cbind(all_tissue_spm,spm)
  names(all_tissue_spm)[ncol(all_tissue_spm)] <- paste("All_Tissue_SPM",time,sep = "_")
}
# Selective Tissue SPM
tissue_list <- c("Heart","Kidney","Liver","Testis")
time_list <- c("4wpc","5wpc","6wpc","7wpc","8wpc","9wpc","10wpc","11wpc","12wpc","13wpc","16wpc","18wpc","19wpc","infant","toddler","youngTeenager","oldTeenager","youngAdult","olderMidAge")
selective_tissue_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (time in time_list) {
  idx <- c()
  for (tissue in tissue_list) {
    tissue.idx <- grep(tissue,names(data_all))
    time.idx <- grep(time,names(data_all))
    idx <- c(idx,intersect(tissue.idx,time.idx))
  }
  spm <- c()
  for (i in data$GeneName) {
    if (i %in% data_all$GeneName) {
      timepoint_data <- data_all[data_all$GeneName == i,idx]
      testis.idx <- grep("Testis",names(timepoint_data))
      v1 <- as.numeric(timepoint_data[,1:ncol(timepoint_data)]) ### From col 2 to the end
      v2 <- rep(0,length(v1))
      v2[testis.idx] <- v1[testis.idx] ### 1 col before actual data
      tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
      tmp.spm <- ifelse(sum(v1[testis.idx]) == 0, 0, tmp.spm) ### Replace NA with 0
      spm <- c(spm,tmp.spm)
    } else {
      spm <- c(spm,NA)
    }
  }
  selective_tissue_spm <- cbind(selective_tissue_spm,spm)
  names(selective_tissue_spm)[ncol(selective_tissue_spm)] <- paste("Selective_Tissue_SPM",time,sep = "_")
}
# Imputation of Tissue SPM
to_imputate_data <- cbind(
  timepoint_mean,
  Alltime_Mean = alltime_mean,
  Prenatal_Mean = prenatal_mean,
  Postnatal_Mean = postnatal_mean,
  temporal_spm,
  Beta = temporal_dynamics,
  all_tissue_spm,
  selective_tissue_spm)
### Rename
names(to_imputate_data) <- paste("ZMBH_Human",names(to_imputate_data),sep = "_")
to_imputate_data$Response <- sn_data$Phenotype.Male %>% 
  gsub("Normal",0,.) %>% 
  gsub("Abnormal",1,.) %>% 
  as.numeric()
### Multiple Imputation of Tissue SPM
md.pattern(to_imputate_data)
imputation_data <- mice(
  to_imputate_data,
  m = 3,
  maxit = 10,
  seed = 777,
  method = "pmm",
  remove.collinear = FALSE,
  remove.constant = FALSE)
### Visualization of Imputation
#summary(imputation_data)
#stripplot(imputation_data,col = c("grey",mdc(2)),pch = c(1,20))
### Fitting to Diagnose
#imputation_data_fit <- lm.mids(Response ~ .,data = imputation_data)
#imputation_data_pooled <- pool(imputation_data_fit)
### Diagnosis of Pooled Imputation
#res.1 <- imputation_data_pooled$glanced;res.1
#res.2 <- cbind(imputation_data_pooled$pooled,summary(imputation_data_pooled))
#pool.r.squared(imputation_data_fit)
### Completion
imputed_data <- complete(imputation_data,action = 2) %>% .[,-ncol(.)]
#write.csv(cbind(data,imputed_data),"./ZMBH_Human.csv")
# Merge
sn_data <- sn_data[,-grep("ZMBH_Human",names(sn_data))] %>% cbind(.,imputed_data)

# 1.3 GTEx Human Testis -----
data <- sn_data[,c(1,grep("GTEx_Human",names(sn_data)))]
data_all <- gtex.all.median
names(data_all)[2] <- "humanGene"
# From Human Symbol to Mouse Symbol
genoinfo <- human2mouse(data_all$humanGene)
# Merge with Symbol Annotation
data_all <- merge(data_all,genoinfo,all.x = TRUE) %>%
  .[!is.na(.$mouseGene),] %>% 
  .[,c(57,3:56)]
names(data_all)[1] <- "GeneName"
### Tackle Dplications
dp_list <- data_all[duplicated(data_all$GeneName),] %>% 
  .[!duplicated(.$GeneName),]
for (i in dp_list$GeneName) {
  tp <- data_all[data_all$GeneName == i,]
  ##### Calculate Shannon’s Entropy
  tp_H <- apply(tp[,2:ncol(tp)],1,function(x) {
    x <- as.numeric(x)
    if (sum(x) == 0) {
      H = 0
    } else {
      Prob <- x/sum(x)
      H <- H(Prob) 
    }
  })
  if (!sum(tp_H == max(tp_H)) == nrow(tp)) {
    tp_max_H <- data_all[rownames(data_all) == names(tp_H[tp_H == max(tp_H)][1]),]
  } else {
    tp_max_H <- data_all[rownames(data_all) == names(tp_H[1]),]
  }
  data_all[data_all$GeneName == i,] <- tp_max_H
}
data_all <- data_all[!duplicated(data_all$GeneName),]
# All Tissue SPM
tissue_list <- gsub("GTEx_Human_","",names(data_all)) %>% 
  gsub("\\(.*\\)","",.) %>% 
  .[-1]
all_tissue_spm <- c()
for (i in data$GeneName) {
  if (i %in% data_all$GeneName) {
    idx <- c()
    for (tissue in tissue_list) {
      tissue.idx <- grep(tissue,names(data_all))
      idx <- c(idx,tissue.idx)
    }
    temp_data <- data_all[data_all$GeneName == i,idx]
    testis.idx <- grep("Testis",names(temp_data))
    v1 <- as.numeric(temp_data[,1:ncol(temp_data)]) ### From col 2 to the end
    v2 <- rep(0,length(v1))
    v2[testis.idx] <- v1[testis.idx] ### 1 col before actual data
    tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
    tmp.spm <- ifelse(sum(v1[testis.idx]) == 0, 0, tmp.spm) ### Replace NA with 0
    all_tissue_spm <- c(all_tissue_spm,tmp.spm)
  } else {
    all_tissue_spm <- c(all_tissue_spm,NA)
  }
}
# Selective Tissue SPM
remove_list <- c(
  "Uterus","Vagina","Ovary","Bladder","Cervix-Ectocervix","Cervix-Endocervix","Prostate",
  "Nerve-Tibial",
  "Brain-Cortex","Brain-Cerebellum","Brain-Frontal Cortex ","Brain-Caudate ","Brain-Nucleus accumbens ","Brain-Putamen ","Brain-Hypothalamus","Brain-Spinal cord ","Brain-Hippocampus","Brain-Anterior cingulate cortex ","Brain-Cerebellar Hemisphere","Brain-Substantia nigra","Brain-Amygdala",
  "Cells-Cultured fibroblasts","Cells-EBV-transformed lymphocytes")
tissue_list <- gsub("GTEx_Human_","",names(data_all)) %>% 
  gsub("\\(.*\\)","",.) %>% 
  .[-1] %>% 
  .[-which(. %in% remove_list)]
selective_tissue_spm <- c()
for (i in data$GeneName) {
  if (i %in% data_all$GeneName) {
    idx <- c()
    for (tissue in tissue_list) {
      tissue.idx <- grep(tissue,names(data_all))
      idx <- c(idx,tissue.idx)
    }
    temp_data <- data_all[data_all$GeneName == i,idx]
    testis.idx <- grep("Testis",names(temp_data))
    v1 <- as.numeric(temp_data[,1:ncol(temp_data)]) ### From col 2 to the end
    v2 <- rep(0,length(v1))
    v2[testis.idx] <- v1[testis.idx] ### 1 col before actual data
    tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
    tmp.spm <- ifelse(sum(v1[testis.idx]) == 0, 0, tmp.spm) ### Replace NA with 0
    selective_tissue_spm <- c(selective_tissue_spm,tmp.spm)
  } else {
    selective_tissue_spm <- c(selective_tissue_spm,NA)
  }
}
# Imputation of Tissue SPM
to_imputate_data <- data.frame(
  All_Tissue_SPM = all_tissue_spm,
  Selective_Tissue_SPM = selective_tissue_spm)
### Rename
names(to_imputate_data) <- paste("GTEx_Human",names(to_imputate_data),sep = "_")
to_imputate_data$Response <- sn_data$Phenotype.Male %>% 
  gsub("Normal",0,.) %>% 
  gsub("Abnormal",1,.) %>% 
  as.numeric()
### Multiple Imputation of Tissue SPM
md.pattern(to_imputate_data)
imputation_data <- mice(
  to_imputate_data,
  m = 3,
  maxit = 10,
  seed = 777,
  method = "pmm",
  remove.collinear = FALSE,
  remove.constant = FALSE)
### Visualization of Imputation
#summary(imputation_data)
#stripplot(imputation_data,col = c("grey",mdc(2)),pch = c(1,20))
### Fitting to Diagnose
#imputation_data_fit <- lm.mids(Response ~ .,data = imputation_data)
#imputation_data_pooled <- pool(imputation_data_fit)
### Diagnosis of Pooled Imputation
#res.1 <- imputation_data_pooled$glanced;res.1
#res.2 <- cbind(imputation_data_pooled$pooled,summary(imputation_data_pooled))
#pool.r.squared(imputation_data_fit)
### Completion
imputed_data <- complete(imputation_data,action = 2) %>% .[,-ncol(.)]
#write.csv(cbind(data,imputed_data),"./GTEx_Human.csv")
# Merge
sn_data <- cbind(sn_data,imputed_data)




#####
##### Single Cell Features -----------------------------------------------------
#####





# 2.1 (Mus Spermatogenesis) GSE107644 -----
data <- sn_data[,c(1,grep("GSE107644_Mouse",names(sn_data)))]
# Cell-Type SPM
celltype_list <- gsub("GSE107644_Mouse_","",names(data))[-1]
celltype_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (celltype in celltype_list) {
  celltype.idx <- which(names(data) == paste("GSE107644_Mouse",celltype,sep = "_"))
  spm <- c()
  for (i in 1:nrow(data)) {
    v1 <- as.numeric(data[i,2:ncol(data)]) ### From col 2 to the end
    v2 <- rep(0,length(v1))
    v2[celltype.idx-1] <- v1[celltype.idx-1] ### 1 col before actual data
    tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
    spm <- c(spm,tmp.spm)
  }
  spm[is.na(spm)] <- 0 ### Replace NA with 0
  celltype_spm <- cbind(celltype_spm,spm)
  names(celltype_spm)[ncol(celltype_spm)] <- paste("CellType_SPM",celltype,sep = "_")
}
# Merge
to_merge_data <- cbind(
  celltype_spm)
### Rename
names(to_merge_data) <- paste("GSE107644_Mouse",names(to_merge_data),sep = "_")
#write.csv(cbind(data,to_merge_data),"./GSE107644_Mouse.csv")
sn_data <- cbind(sn_data,to_merge_data)

# 2.2 (Mus Somat) GSE112393 -----
data <- sn_data[,c(1,grep("GSE112393_Mouse",names(sn_data)))]
# Cell-Type SPM
celltype_list <- gsub("GSE112393_Mouse_","",names(data))[-1]
celltype_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (celltype in celltype_list) {
  celltype.idx <- which(names(data) == paste("GSE112393_Mouse",celltype,sep = "_"))
  spm <- c()
  for (i in 1:nrow(data)) {
    v1 <- as.numeric(data[i,2:ncol(data)]) ### From col 2 to the end
    v2 <- rep(0,length(v1))
    v2[celltype.idx-1] <- v1[celltype.idx-1] ### 1 col before actual data
    tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
    spm <- c(spm,tmp.spm)
  }
  spm[is.na(spm)] <- 0 ### Replace NA with 0
  celltype_spm <- cbind(celltype_spm,spm)
  names(celltype_spm)[ncol(celltype_spm)] <- paste("CellType_SPM",celltype,sep = "_")
}
# Merge
to_merge_data <- cbind(
  celltype_spm)
### Rename
names(to_merge_data) <- paste("GSE112393_Mouse",names(to_merge_data),sep = "_")
#write.csv(cbind(data,to_merge_data),"./GSE112393_Mouse.csv")
sn_data <- cbind(sn_data,to_merge_data)

# 2.3 (Mus Gonad) GSE136220 -----
data <- sn_data[,c(1,grep("GSE136220_Mouse",names(sn_data)))]
# Cell-Type SPM
celltype_list <- gsub("GSE136220_Mouse_","",names(data))[-1]
celltype_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (celltype in celltype_list) {
  celltype.idx <- which(names(data) == paste("GSE136220_Mouse",celltype,sep = "_"))
  spm <- c()
  for (i in 1:nrow(data)) {
    v1 <- as.numeric(data[i,2:ncol(data)]) ### From col 2 to the end
    v2 <- rep(0,length(v1))
    v2[celltype.idx-1] <- v1[celltype.idx-1] ### 1 col before actual data
    tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
    spm <- c(spm,tmp.spm)
  }
  spm[is.na(spm)] <- 0 ### Replace NA with 0
  celltype_spm <- cbind(celltype_spm,spm)
  names(celltype_spm)[ncol(celltype_spm)] <- paste("CellType_SPM",celltype,sep = "_")
}
# Merge
to_merge_data <- cbind(
  celltype_spm)
### Rename
names(to_merge_data) <- paste("GSE136220_Mouse",names(to_merge_data),sep = "_")
#write.csv(cbind(data,to_merge_data),"./GSE136220_Mouse.csv")
sn_data <- cbind(sn_data,to_merge_data)

# 2.4 (Mus Gonad) GSE148032 -----
data <- sn_data[,c(1,grep("GSE148032_Mouse",names(sn_data)))]
# Cell-Type SPM
celltype_list <- gsub("GSE148032_Mouse_","",names(data))[-1]
celltype_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (celltype in celltype_list) {
  celltype.idx <- which(names(data) == paste("GSE148032_Mouse",celltype,sep = "_"))
  spm <- c()
  for (i in 1:nrow(data)) {
    v1 <- as.numeric(data[i,2:ncol(data)]) ### From col 2 to the end
    v2 <- rep(0,length(v1))
    v2[celltype.idx-1] <- v1[celltype.idx-1] ### 1 col before actual data
    tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
    spm <- c(spm,tmp.spm)
  }
  spm[is.na(spm)] <- 0 ### Replace NA with 0
  celltype_spm <- cbind(celltype_spm,spm)
  names(celltype_spm)[ncol(celltype_spm)] <- paste("CellType_SPM",celltype,sep = "_")
}
# Merge
to_merge_data <- cbind(
  celltype_spm)
### Rename
names(to_merge_data) <- paste("GSE148032_Mouse",names(to_merge_data),sep = "_")
#write.csv(cbind(data,to_merge_data),"./GSE148032_Mouse.csv")
sn_data <- cbind(sn_data,to_merge_data)

# 2.5 (Hom Spermatogenesis) GSE106487 -----
data <- sn_data[,c(1,grep("GSE106487_Human",names(sn_data)))]
# Cell-Type SPM
celltype_list <- gsub("GSE106487_Human_","",names(data))[-1]
celltype_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (celltype in celltype_list) {
  celltype.idx <- which(names(data) == paste("GSE106487_Human",celltype,sep = "_"))
  spm <- c()
  for (i in 1:nrow(data)) {
    v1 <- as.numeric(data[i,2:ncol(data)]) ### From col 2 to the end
    v2 <- rep(0,length(v1))
    v2[celltype.idx-1] <- v1[celltype.idx-1] ### 1 col before actual data
    tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
    spm <- c(spm,tmp.spm)
  }
  spm[is.na(spm)] <- 0 ### Replace NA with 0
  celltype_spm <- cbind(celltype_spm,spm)
  names(celltype_spm)[ncol(celltype_spm)] <- paste("CellType_SPM",celltype,sep = "_")
}
# Merge
to_merge_data <- cbind(
  celltype_spm)
### Rename
names(to_merge_data) <- paste("GSE106487_Human",names(to_merge_data),sep = "_")
#write.csv(cbind(data,to_merge_data),"./GSE106487_Human.csv")
sn_data <- cbind(sn_data,to_merge_data)

# 2.6 (Hom Spermatogenesis) Lab -----
data <- sn_data[,c(1,grep("Lab_Human",names(sn_data)))]
# Cell-Type SPM
celltype_list <- gsub("Lab_Human_","",names(data))[-1]
celltype_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (celltype in celltype_list) {
  celltype.idx <- which(names(data) == paste("Lab_Human",celltype,sep = "_"))
  spm <- c()
  for (i in 1:nrow(data)) {
    v1 <- as.numeric(data[i,2:ncol(data)]) ### From col 2 to the end
    v2 <- rep(0,length(v1))
    v2[celltype.idx-1] <- v1[celltype.idx-1] ### 1 col before actual data
    tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
    spm <- c(spm,tmp.spm)
  }
  spm[is.na(spm)] <- 0 ### Replace NA with 0
  celltype_spm <- cbind(celltype_spm,spm)
  names(celltype_spm)[ncol(celltype_spm)] <- paste("CellType_SPM",celltype,sep = "_")
}
# Merge
to_merge_data <- cbind(
  celltype_spm)
### Rename
names(to_merge_data) <- paste("Lab_Human",names(to_merge_data),sep = "_")
#write.csv(cbind(data,to_merge_data),"./Lab_Human.csv")
sn_data <- cbind(sn_data,to_merge_data)

# 2.7 (Hom Puberty) GSE134144 -----
data <- sn_data[,c(1,grep("GSE134144_Human",names(sn_data)))]
# Cell-Type SPM
celltype_list <- gsub("GSE134144_Human_","",names(data))[-1]
celltype_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (celltype in celltype_list) {
  celltype.idx <- which(names(data) == paste("GSE134144_Human",celltype,sep = "_"))
  spm <- c()
  for (i in 1:nrow(data)) {
    v1 <- as.numeric(data[i,2:ncol(data)]) ### From col 2 to the end
    v2 <- rep(0,length(v1))
    v2[celltype.idx-1] <- v1[celltype.idx-1] ### 1 col before actual data
    tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
    spm <- c(spm,tmp.spm)
  }
  spm[is.na(spm)] <- 0 ### Replace NA with 0
  celltype_spm <- cbind(celltype_spm,spm)
  names(celltype_spm)[ncol(celltype_spm)] <- paste("CellType_SPM",celltype,sep = "_")
}
# Merge
to_merge_data <- cbind(
  celltype_spm)
### Rename
names(to_merge_data) <- paste("GSE134144_Human",names(to_merge_data),sep = "_")
#write.csv(cbind(data,to_merge_data),"./GSE134144_Human.csv")
sn_data <- cbind(sn_data,to_merge_data)

# 2.8 (Hom Somat) GSE124263 -----
data <- sn_data[,c(1,grep("GSE124263_Human",names(sn_data)))]
# Cell-Type SPM
celltype_list <- gsub("GSE124263_Human_","",names(data))[-1]
celltype_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (celltype in celltype_list) {
  celltype.idx <- which(names(data) == paste("GSE124263_Human",celltype,sep = "_"))
  spm <- c()
  for (i in 1:nrow(data)) {
    v1 <- as.numeric(data[i,2:ncol(data)]) ### From col 2 to the end
    v2 <- rep(0,length(v1))
    v2[celltype.idx-1] <- v1[celltype.idx-1] ### 1 col before actual data
    tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
    spm <- c(spm,tmp.spm)
  }
  spm[is.na(spm)] <- 0 ### Replace NA with 0
  celltype_spm <- cbind(celltype_spm,spm)
  names(celltype_spm)[ncol(celltype_spm)] <- paste("CellType_SPM",celltype,sep = "_")
}
# Merge
to_merge_data <- cbind(
  celltype_spm)
### Rename
names(to_merge_data) <- paste("GSE124263_Human",names(to_merge_data),sep = "_")
#write.csv(cbind(data,to_merge_data),"./GSE124263_Human.csv")
sn_data <- cbind(sn_data,to_merge_data)

# 2.9 (Hom Somat) GSE112013 -----
data <- sn_data[,c(1,grep("GSE112013_Human",names(sn_data)))]
# Cell-Type SPM
celltype_list <- gsub("GSE112013_Human_","",names(data))[-1]
celltype_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (celltype in celltype_list) {
  celltype.idx <- which(names(data) == paste("GSE112013_Human",celltype,sep = "_"))
  spm <- c()
  for (i in 1:nrow(data)) {
    v1 <- as.numeric(data[i,2:ncol(data)]) ### From col 2 to the end
    v2 <- rep(0,length(v1))
    v2[celltype.idx-1] <- v1[celltype.idx-1] ### 1 col before actual data
    tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
    spm <- c(spm,tmp.spm)
  }
  spm[is.na(spm)] <- 0 ### Replace NA with 0
  celltype_spm <- cbind(celltype_spm,spm)
  names(celltype_spm)[ncol(celltype_spm)] <- paste("CellType_SPM",celltype,sep = "_")
}
# Merge
to_merge_data <- cbind(
  celltype_spm)
### Rename
names(to_merge_data) <- paste("GSE112013_Human",names(to_merge_data),sep = "_")
#write.csv(cbind(data,to_merge_data),"./GSE112013_Human.csv")
sn_data <- cbind(sn_data,to_merge_data)

# 2.10 (Hom Somat) GSE142585 -----
data <- sn_data[,c(1,grep("GSE142585_Human",names(sn_data)))]
# Cell-Type SPM
celltype_list <- gsub("GSE142585_Human_","",names(data))[-1]
celltype_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (celltype in celltype_list) {
  celltype.idx <- which(names(data) == paste("GSE142585_Human",celltype,sep = "_"))
  spm <- c()
  for (i in 1:nrow(data)) {
    v1 <- as.numeric(data[i,2:ncol(data)]) ### From col 2 to the end
    v2 <- rep(0,length(v1))
    v2[celltype.idx-1] <- v1[celltype.idx-1] ### 1 col before actual data
    tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
    spm <- c(spm,tmp.spm)
  }
  spm[is.na(spm)] <- 0 ### Replace NA with 0
  celltype_spm <- cbind(celltype_spm,spm)
  names(celltype_spm)[ncol(celltype_spm)] <- paste("CellType_SPM",celltype,sep = "_")
}
# Merge
to_merge_data <- cbind(
  celltype_spm)
### Rename
names(to_merge_data) <- paste("GSE142585_Human",names(to_merge_data),sep = "_")
#write.csv(cbind(data,to_merge_data),"./GSE142585_Human.csv")
sn_data <- cbind(sn_data,to_merge_data)

# 2.11 (Hom Gonad) GSE86146 -----
data <- sn_data[,c(1,grep("GSE86146_Human",names(sn_data)))]
# Cell-Type SPM
celltype_list <- gsub("GSE86146_Human_","",names(data))[-1]
celltype_spm <- as.data.frame(array(dim = c(nrow(data),0)))
for (celltype in celltype_list) {
  celltype.idx <- which(names(data) == paste("GSE86146_Human",celltype,sep = "_"))
  spm <- c()
  for (i in 1:nrow(data)) {
    v1 <- as.numeric(data[i,2:ncol(data)]) ### From col 2 to the end
    v2 <- rep(0,length(v1))
    v2[celltype.idx-1] <- v1[celltype.idx-1] ### 1 col before actual data
    tmp.spm <- (v1 %*% v2)/(modulus(v1)*modulus(v2))
    spm <- c(spm,tmp.spm)
  }
  spm[is.na(spm)] <- 0 ### Replace NA with 0
  celltype_spm <- cbind(celltype_spm,spm)
  names(celltype_spm)[ncol(celltype_spm)] <- paste("CellType_SPM",celltype,sep = "_")
}
# Merge
to_merge_data <- cbind(
  celltype_spm)
### Rename
names(to_merge_data) <- paste("GSE86146_Human",names(to_merge_data),sep = "_")
#write.csv(cbind(data,to_merge_data),"./GSE86146_Human.csv")
sn_data <- cbind(sn_data,to_merge_data)

# Display -----
dim(sn_data)
table(sn_data$Phenotype.Male)




#####
##### Save -----------------------------------------------------------
#####





save(sn_data,file = "All_Merged, Imputated, Calculated Features.RData")
