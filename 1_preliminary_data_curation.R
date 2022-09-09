source("./0_presets.R")

#####
##### Genelist -----------------------------------------------------------------
#####





noa <- read.csv("../../data/genelist_NOA.csv") %>% .[,c(1:3)]
normal <- read.csv("../../data/genelist_Normal.csv") %>% .[,c(1:3)]
# Remake Genelist
gl <- rbind(noa,normal)
gl$Phenotype.Male <- gsub("NOA","Abnormal",gl$Phenotype.Male) ### Rename Abnormal
gl$Phenotype.Male <- gsub("Oligozoospermia","Abnormal",gl$Phenotype.Male) ### Rename Abnormal

#write.csv(gl,"../../data/genelist_All.csv",row.names = FALSE)
table(gl$Phenotype.Male)





#####
##### Genomic Features ---------------------------------------------------------
#####





# 2.1 RVIS Score
RVIS <- read.csv("../../data/pLI & lethal.csv")
names(RVIS)[1] <- "humanGene"
# From Human Symbol to Mouse Symbol
genoinfo <- human2mouse(RVIS$humanGene)
# Merge with Symbol Annotation
RVIS <- merge(RVIS,genoinfo,all.x = TRUE) %>%
  .[!is.na(.$mouseGene),] %>% 
  .[,c(15,9)]
names(RVIS)[1] <- "GeneName"
RVIS[RVIS == "."] <- NA
RVIS$RVIS_ExAc_score <- as.numeric(RVIS$RVIS_ExAc_score)
# Tackle Dplications
dp_list <- RVIS[duplicated(RVIS$GeneName),] %>% 
  .[!duplicated(.$GeneName),]
for (i in dp_list$GeneName) {
  tp_min_RVIS <- min(RVIS[RVIS$GeneName == i & !RVIS$RVIS_ExAc_score == ".",2])
  RVIS[RVIS$GeneName == i,2] <- tp_min_RVIS
}
RVIS <- RVIS[!duplicated(RVIS$GeneName),]

# 2.2 Shet Score
Shet <- read.csv("../../data/shet.csv")
names(Shet)[1] <- "humanGene"
# From Human Symbol to Mouse Symbol
genoinfo <- human2mouse(Shet$humanGene)
# Merge with Symbol Annotation (leave duplication to tackle when merge with gl)
Shet <- merge(Shet,genoinfo,all.x = TRUE) %>%
  .[!is.na(.$mouseGene),] %>% 
  .[,c(19,2)]
names(Shet)[1] <- "GeneName"
Shet[Shet == "."] <- NA
Shet[,2] <- as.numeric(Shet[,2])
# Tackle Dplications
dp_list <- Shet[duplicated(Shet$GeneName),] %>% 
  .[!duplicated(.$GeneName),]
for (i in dp_list$GeneName) {
  tp_max_Shet <- max(Shet[Shet$GeneName == i,2])
  Shet[Shet$GeneName == i,2] <- tp_max_Shet
}
Shet <- Shet[!duplicated(Shet$GeneName),]

# 2.3 Constraint
constraint <- read.table("../../data/forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt",header = TRUE)
names(constraint)[2] <- "humanGene"
### From Human Symbol to Mouse Symbol
genoinfo <- human2mouse(constraint$humanGene)
### Merge with Symbol Annotation
constraint <- merge(constraint,genoinfo,all.x = TRUE) %>%
  .[!is.na(.$mouseGene),] %>% 
  .[,c(24,17:19,23,20)]
names(constraint)[1] <- "GeneName"
### Tackle Dplications
dp_list <- constraint[duplicated(constraint$GeneName),] %>% 
  .[!duplicated(.$GeneName),]
for (i in dp_list$GeneName) {
  tp <- constraint[constraint$GeneName == i,-5]
  ##### Calculate Shannon’s Entropy
  tp_H <- apply(tp[,2:ncol(tp)],1,function(x) {
    x <- as.numeric(x) %>% abs()
    if (sum(x) == 0) {
      H = 0
    } else {
      Prob <- x/sum(x)
      H <- H(Prob)
    }
  })
  if (!sum(tp_H == max(tp_H)) == nrow(tp)) {
    tp_max_H <- constraint[rownames(constraint) == names(tp_H[tp_H == max(tp_H)][1]),]
  } else {
    tp_max_H <- constraint[rownames(constraint) == names(tp_H[1]),]
  }
  constraint[constraint$GeneName == i,] <- tp_max_H
}
constraint <- constraint[!duplicated(constraint$GeneName),]





#####
##### mRNA Expression Features -------------------------------------------------
#####





# 3.1 ZMBH Mouse Testis
mouse.data <- mouse.testis
# Rename
names(mouse.data)[2:ncol(mouse.data)] <- paste("ZMBH_Mouse",names(mouse.data)[2:ncol(mouse.data)],sep = "_")
# Statistics
gsub("Testis.","",names(mouse.testis))[-1] ### Sample Information
nrow(mouse.data) ### Num of Genes Detected
sum(gl$Ensembl.ID %in% mouse.data$Ensembl.ID)/nrow(gl) ### Coverage of Genelist

# 3.2 ZMBH Human Testis
human.data <- human.testis
# From Human Ensembl to Human Symbol
genoinfo1 <- bitr(human.data$Ensembl.ID,fromType = 'ENSEMBL',toType = c('SYMBOL'),OrgDb = 'org.Hs.eg.db')
names(genoinfo1) <- c("Ensembl.ID","humanGene")
# From Human Symbol to Mouse Symbol
genoinfo2 <- human2mouse(genoinfo1$humanGene)
genoinfo <- merge(genoinfo1,genoinfo2,all.x = TRUE) %>% .[,2:3]
# Merge with Symbol Annotation
human.data <- merge(human.data,genoinfo,all.x = TRUE) %>%
  .[!is.na(.$mouseGene),] %>% 
  .[,c(41,2:40)]
names(human.data)[1] <- "GeneName"
# Tackle Dplications
dp_list <- human.data[duplicated(human.data$GeneName),] %>% 
  .[!duplicated(.$GeneName),]
for (i in dp_list$GeneName) {
  tp <- human.data[human.data$GeneName == i,]
  ### Calculate Shannon’s Entropy
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
    tp_max_H <- human.data[rownames(human.data) == names(tp_H[tp_H == max(tp_H)][1]),]
  } else {
    tp_max_H <- human.data[rownames(human.data) == names(tp_H[1]),]
  }
  human.data[human.data$GeneName == i,] <- tp_max_H
}
human.data <- human.data[!duplicated(human.data$GeneName),]
# Rename
names(human.data)[2:ncol(human.data)] <- paste("ZMBH_Human",names(human.data)[2:ncol(human.data)],sep = "_")
# Statistics
gsub("Testis.","",names(human.testis))[-1] ### Sample Information
nrow(human.testis) ### Num of Genes Detected
nrow(human.data) ### Num of Genes Translated
sum(gl$GeneName %in% human.data$GeneName)/nrow(gl) ### Coverage of Genelist

# 3.3 GTEx Human Testis
gtex.human.data <- gtex.all.median[,c(2,28)]
names(gtex.human.data)[1] <- "humanGene"
# From Human Symbol to Mouse Symbol
genoinfo <- human2mouse(gtex.human.data$humanGene)
# Merge with Symbol Annotation
gtex.human.data <- merge(gtex.human.data,genoinfo,all.x = TRUE) %>%
  .[!is.na(.$mouseGene),] %>% 
  .[,c(3,2)]
names(gtex.human.data)[1] <- "GeneName"
# Tackle Dplications
dp_list <- gtex.human.data[duplicated(gtex.human.data$GeneName),] %>% 
  .[!duplicated(.$GeneName),]
for (i in dp_list$GeneName) {
  tp <- gtex.human.data[gtex.human.data$GeneName == i,]
  ### Calculate Mean
  tp_mean <- cbind(
    GeneName = i,
    GTEx_Human_Testis = mean(tp$GTEx_Human_Testis))
  gtex.human.data[gtex.human.data$GeneName == i,2] <- mean(tp$GTEx_Human_Testis)
}
gtex.human.data <- gtex.human.data[!duplicated(gtex.human.data$GeneName),]
# Statistics
nrow(gtex.all.median) ### Num of Genes Detected
nrow(gtex.human.data) ### Num of Genes Translated
sum(gl$GeneName %in% gtex.human.data$GeneName)/nrow(gl) ### Coverage of Genelist





#####
##### Single Cell Features -----------------------------------------------------
#####





# 4.1 (Mus Spermatogenesis) GSE107644
# Statistics
gsub("GSE107644_Mouse_","",names(gse107644))[-1] ### Sample Information
nrow(gse107644) ### Num of Genes Detected
sum(gl$GeneName %in% gse107644$GeneName)/nrow(gl) ### Coverage of Genelist

# 4.2 (Mus Somat) GSE112393
# Statistics
gsub("GSE112393_Mouse_","",names(gse112393))[-1] ### Sample Information
nrow(gse112393) ### Num of Genes Detected
sum(gl$GeneName %in% gse112393$GeneName)/nrow(gl) ### Coverage of Genelist

# 4.3 (Mus Gonad) GSE136220
gse136220$Ensembl.ID <- substr(gse136220$Ensembl.ID,1,18)
# Statistics
gsub("GSE136220_Mouse_","",names(gse136220))[-1] ### Sample Information
nrow(gse136220) ### Num of Genes Detected
sum(gl$Ensembl.ID %in% gse136220$Ensembl.ID)/nrow(gl) ### Coverage of Genelist

# 4.4 (Mus Gonad) GSE148032
# Statistics
gsub("GSE148032_Mouse_","",names(gse148032))[-1] ### Sample Information
nrow(gse148032) ### Num of Genes Detected
sum(gl$GeneName %in% gse148032$GeneName)/nrow(gl) ### Coverage of Genelist

# 4.5 (Mus E10.5-E16.5 Gonad) GSE97519
# Statistics
gsub("GSE97519_Mouse_","",names(gse97519))[-1] ### Sample Information
nrow(gse97519) ### Num of Genes Detected
sum(gl$GeneName %in% gse97519$GeneName)/nrow(gl) ### Coverage of Genelist

# 4.6 (Mus Adult SPG) GSE108974
# Statistics
gsub("GSE108974_Mouse_","",names(gse108974))[-1] ### Sample Information
nrow(gse108974) ### Num of Genes Detected
sum(gl$GeneName %in% gse108974$GeneName)/nrow(gl) ### Coverage of Genelist

# 4.7 (Mus P6 SPG) GSE108970
# Statistics
gsub("GSE108970_Mouse_","",names(gse108970))[-1] ### Sample Information
nrow(gse108970) ### Num of Genes Detected
sum(gl$GeneName %in% gse108970$GeneName)/nrow(gl) ### Coverage of Genelist

# 4.8 (Mus P3&7 SPG) GSE82174
# Statistics
gsub("GSE82174_Mouse_","",names(gse82174))[-1] ### Sample Information
nrow(gse82174) ### Num of Genes Detected
sum(gl$GeneName %in% gse82174$GeneName)/nrow(gl) ### Coverage of Genelist

# 4.9 (Mus P5.5 SCS) GSE107711
# Statistics
gsub("GSE107711_Mouse_","",names(gse107711))[-1] ### Sample Information
nrow(gse107711) ### Num of Genes Detected
sum(gl$GeneName %in% gse107711$GeneName)/nrow(gl) ### Coverage of Genelist

# 4.10 (Hom Spermatogenesis) GSE106487
human.scs <- gse106487
names(human.scs)[1] <- "humanGene"
# From Human Symbol to Mouse Symbol
genoinfo <- human2mouse(human.scs$humanGene)
# Merge with Symbol Annotation
human.scs <- merge(human.scs,genoinfo,all.x = TRUE) %>%
  .[!is.na(.$mouseGene),] %>% 
  .[,c(18,2:17)]
names(human.scs)[1] <- "GeneName"
# Tackle Dplications
dp_list <- human.scs[duplicated(human.scs$GeneName),] %>% 
  .[!duplicated(.$GeneName),]
for (i in dp_list$GeneName) {
  tp <- human.scs[human.scs$GeneName == i,]
  ### Calculate Shannon’s Entropy
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
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[tp_H == max(tp_H)][1]),]
  } else {
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[1]),]
  }
  human.scs[human.scs$GeneName == i,] <- tp_max_H
}
human.scs <- human.scs[!duplicated(human.scs$GeneName),]
# Statistics
gsub("GSE106487_Human_","",names(human.scs))[-1] ### Sample Information
nrow(gse106487) ### Num of Genes Detected
nrow(human.scs) ### Num of Genes Translated
sum(gl$GeneName %in% human.scs$GeneName)/nrow(gl) ### Coverage of Genelist
gse106487 <- human.scs

# 4.11 (Hom Spermatogenesis) Lab
human.scs <- lab.SCSeq
names(human.scs)[1] <- "humanGene"
# From Human Symbol to Mouse Symbol
genoinfo <- human2mouse(human.scs$humanGene)
# Merge with Symbol Annotation
human.scs <- merge(human.scs,genoinfo,all.x = TRUE) %>%
  .[!is.na(.$mouseGene),] %>% 
  .[,c(11,2:10)]
names(human.scs)[1] <- "GeneName"
# Tackle Dplications
dp_list <- human.scs[duplicated(human.scs$GeneName),] %>% 
  .[!duplicated(.$GeneName),]
for (i in dp_list$GeneName) {
  tp <- human.scs[human.scs$GeneName == i,]
  ### Calculate Shannon’s Entropy
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
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[tp_H == max(tp_H)][1]),]
  } else {
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[1]),]
  }
  human.scs[human.scs$GeneName == i,] <- tp_max_H
}
human.scs <- human.scs[!duplicated(human.scs$GeneName),]
# Statistics
gsub("Lab_Human_","",names(human.scs))[-1] ### Sample Information
nrow(lab.SCSeq) ### Num of Genes Detected
nrow(human.scs) ### Num of Genes Translated
sum(gl$GeneName %in% human.scs$GeneName)/nrow(gl) ### Coverage of Genelist
lab.SCSeq <- human.scs

# 4.12 (Hom Puberty) GSE134144
human.scs <- gse134144
names(human.scs)[1] <- "humanGene"
# From Human Symbol to Mouse Symbol
genoinfo <- human2mouse(human.scs$humanGene)
# Merge with Symbol Annotation
human.scs <- merge(human.scs,genoinfo,all.x = TRUE) %>%
  .[!is.na(.$mouseGene),] %>% 
  .[,c(32,2:31)]
names(human.scs)[1] <- "GeneName"
# Tackle Dplications
dp_list <- human.scs[duplicated(human.scs$GeneName),] %>% 
  .[!duplicated(.$GeneName),]
for (i in dp_list$GeneName) {
  tp <- human.scs[human.scs$GeneName == i,]
  ### Calculate Shannon’s Entropy
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
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[tp_H == max(tp_H)][1]),]
  } else {
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[1]),]
  }
  human.scs[human.scs$GeneName == i,] <- tp_max_H
}
human.scs <- human.scs[!duplicated(human.scs$GeneName),]
# Statistics
gsub("GSE134144_Human_","",names(human.scs))[-1] ### Sample Information
nrow(gse134144) ### Num of Genes Detected
nrow(human.scs) ### Num of Genes Translated
sum(gl$GeneName %in% human.scs$GeneName)/nrow(gl) ### Coverage of Genelist
gse134144 <- human.scs

# 4.13 (Hom Somat) GSE124263
human.scs <- gse124263
# From Human Ensembl to Human Symbol
genoinfo1 <- bitr(human.scs$Ensembl.ID,fromType = 'ENSEMBL',toType = c('SYMBOL'),OrgDb = 'org.Hs.eg.db')
names(genoinfo1) <- c("Ensembl.ID","humanGene")
# From Human Symbol to Mouse Symbol
genoinfo2 <- human2mouse(genoinfo1$humanGene)
genoinfo <- merge(genoinfo1,genoinfo2,all.x = TRUE) %>% .[,2:3]
# Merge with Symbol Annotation
human.scs <- merge(human.scs,genoinfo,all.x = TRUE) %>%
  .[!is.na(.$mouseGene),] %>% 
  .[,c(20,2:19)]
names(human.scs)[1] <- "GeneName"
# Tackle Dplications
dp_list <- human.scs[duplicated(human.scs$GeneName),] %>% 
  .[!duplicated(.$GeneName),]
for (i in dp_list$GeneName) {
  tp <- human.scs[human.scs$GeneName == i,]
  ### Calculate Shannon’s Entropy
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
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[tp_H == max(tp_H)][1]),]
  } else {
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[1]),]
  }
  human.scs[human.scs$GeneName == i,] <- tp_max_H
}
human.scs <- human.scs[!duplicated(human.scs$GeneName),]
# Statistics
gsub("GSE124263_Human_","",names(human.scs))[-1] ### Sample Information
nrow(gse124263) ### Num of Genes Detected
nrow(human.scs) ### Num of Genes Translated
sum(gl$GeneName %in% human.scs$GeneName)/nrow(gl) ### Coverage of Genelist
gse124263 <- human.scs

# 4.14 (Hom Somat) GSE112013
human.scs <- gse112013
names(human.scs)[1] <- "humanGene"
# From Human Symbol to Mouse Symbol
genoinfo <- human2mouse(human.scs$humanGene)
# Merge with Symbol Annotation
human.scs <- merge(human.scs,genoinfo,all.x = TRUE) %>%
  .[!is.na(.$mouseGene),] %>% 
  .[,c(15,2:14)]
names(human.scs)[1] <- "GeneName"
# Tackle Dplications
dp_list <- human.scs[duplicated(human.scs$GeneName),] %>% 
  .[!duplicated(.$GeneName),]
for (i in dp_list$GeneName) {
  tp <- human.scs[human.scs$GeneName == i,]
  ### Calculate Shannon’s Entropy
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
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[tp_H == max(tp_H)][1]),]
  } else {
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[1]),]
  }
  human.scs[human.scs$GeneName == i,] <- tp_max_H
}
human.scs <- human.scs[!duplicated(human.scs$GeneName),]
# Statistics
gsub("GSE112013_Human_","",names(human.scs))[-1] ### Sample Information
nrow(gse112013) ### Num of Genes Detected
nrow(human.scs) ### Num of Genes Translated
sum(gl$GeneName %in% human.scs$GeneName)/nrow(gl) ### Coverage of Genelist
gse112013 <- human.scs

# 4.15 (Hom Somat) GSE142585
human.scs <- gse142585
names(human.scs)[1] <- "humanGene"
# From Human Symbol to Mouse Symbol
genoinfo <- human2mouse(human.scs$humanGene)
# Merge with Symbol Annotation
human.scs <- merge(human.scs,genoinfo,all.x = TRUE) %>%
  .[!is.na(.$mouseGene),] %>% 
  .[,c(13,2:12)]
names(human.scs)[1] <- "GeneName"
# Tackle Dplications
dp_list <- human.scs[duplicated(human.scs$GeneName),] %>% 
  .[!duplicated(.$GeneName),]
for (i in dp_list$GeneName) {
  tp <- human.scs[human.scs$GeneName == i,]
  ### Calculate Shannon’s Entropy
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
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[tp_H == max(tp_H)][1]),]
  } else {
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[1]),]
  }
  human.scs[human.scs$GeneName == i,] <- tp_max_H
}
human.scs <- human.scs[!duplicated(human.scs$GeneName),]
# Statistics
gsub("GSE142585_Human_","",names(human.scs))[-1] ### Sample Information
nrow(gse142585) ### Num of Genes Detected
nrow(human.scs) ### Num of Genes Translated
sum(gl$GeneName %in% human.scs$GeneName)/nrow(gl) ### Coverage of Genelist
gse142585 <- human.scs

# 4.16 (Hom Gonad) GSE86146
human.scs <- gse86146
names(human.scs)[1] <- "humanGene"
# From Human Symbol to Mouse Symbol
genoinfo <- human2mouse(human.scs$humanGene)
# Merge with Symbol Annotation
human.scs <- merge(human.scs,genoinfo,all.x = TRUE) %>%
  .[!is.na(.$mouseGene),] %>% 
  .[,c(80,2:79)]
names(human.scs)[1] <- "GeneName"
# Tackle Dplications
dp_list <- human.scs[duplicated(human.scs$GeneName),] %>% 
  .[!duplicated(.$GeneName),]
for (i in dp_list$GeneName) {
  tp <- human.scs[human.scs$GeneName == i,]
  ### Calculate Shannon’s Entropy
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
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[tp_H == max(tp_H)][1]),]
  } else {
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[1]),]
  }
  human.scs[human.scs$GeneName == i,] <- tp_max_H
}
human.scs <- human.scs[!duplicated(human.scs$GeneName),]
# Statistics
gsub("GSE86146_Human_","",names(human.scs))[-1] ### Sample Information
nrow(gse86146) ### Num of Genes Detected
nrow(human.scs) ### Num of Genes Translated
sum(gl$GeneName %in% human.scs$GeneName)/nrow(gl) ### Coverage of Genelist
gse86146 <- human.scs

# 4.17 (Hom Adult SPG) GSE108977
human.scs <- gse108977
names(human.scs)[1] <- "humanGene"
# From Human Symbol to Mouse Symbol
genoinfo <- human2mouse(human.scs$humanGene)
# Merge with Symbol Annotation
human.scs <- merge(human.scs,genoinfo,all.x = TRUE) %>%
  .[!is.na(.$mouseGene),] %>% 
  .[,c(11,2:10)]
names(human.scs)[1] <- "GeneName"
# Tackle Dplications
dp_list <- human.scs[duplicated(human.scs$GeneName),] %>% 
  .[!duplicated(.$GeneName),]
for (i in dp_list$GeneName) {
  tp <- human.scs[human.scs$GeneName == i,]
  ### Calculate Shannon’s Entropy
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
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[tp_H == max(tp_H)][1]),]
  } else {
    tp_max_H <- human.scs[rownames(human.scs) == names(tp_H[1]),]
  }
  human.scs[human.scs$GeneName == i,] <- tp_max_H
}
human.scs <- human.scs[!duplicated(human.scs$GeneName),]
# Statistics
gsub("GSE108977_Human_","",names(human.scs))[-1] ### Sample Information
nrow(gse108977) ### Num of Genes Detected
nrow(human.scs) ### Num of Genes Translated
sum(gl$GeneName %in% human.scs$GeneName)/nrow(gl) ### Coverage of Genelist
gse108977 <- human.scs





#####
##### Save ---------------------------------------------------------------------
#####





save(
  noa,normal,gl,
  RVIS,Shet,constraint,
  mouse.data,human.data,gtex.human.data,
  gse107644,gse112393,gse136220,gse148032,gse97519,gse108974,gse108970,gse82174,gse107711,
  gse106487,lab.SCSeq,gse134144,gse124263,gse112013,gse142585,gse86146,gse108977,
  file = "All_Unmerged, Unimputated, Uncalculated Features.RData")