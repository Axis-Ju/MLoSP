source("./0_presets.R")
#load("./All_Unmerged, Unimputated, Uncalculated Features.RData")
#load("./All_Merged, Imputated, Calculated Features.RData")
load("./OthN_Unmerged, Unimputated, Uncalculated Features.RData")
load("./OthN_Merged, Imputated, Calculated Features.RData")

names(sn_data) <- gsub("\\s","_",names(sn_data)) ### Replace space with _ to avoid str2lang error
names(sn_data) <- gsub("_\\&_","_",names(sn_data)) ### Replace & with _ to avoid str2lang error
names(sn_data) <- gsub("#","_",names(sn_data)) ### Replace # with _ to avoid formula.character error
names(sn_data) <- gsub("-","_",names(sn_data)) ### Replace - with _ to avoid eval(predvars, data, env) error
# ##### Temp
# sn_data$Ensembl.ID <- substr(sn_data$Ensembl.ID,1,18)
# sn_data$Phenotype.Male <- NULL
# sn_data <- merge(gl,sn_data)
# sn_data <- sn_data[!duplicated(sn_data$Ensembl.ID),]
# ##### Temp E
# sn_data <- sn_data[!is.na(sn_data$Phenotype.Male),]

# 01.Feature Selection -----
# Screening Based on Difference



### C1.A Recoding Outliers
sn_data.m <- data.frame(apply(sn_data[,-(1:3)],2,outHigh))
sn_data.m <- data.frame(apply(sn_data.m,2,outLow))
sn_data.m <- data.frame(sn_data[,1:3],sn_data.m)
### C1.B No Recoding Outliers
#sn_data.m <- sn_data



### Calculating and Visualization
plist <- list()
wilcoxon_list <- as.data.frame(array(dim = c(0,2)))
names(wilcoxon_list) <- c("Feature","PValue")
for (i in 4:ncol(sn_data.m)) {
  tp <- reshape2::melt(sn_data.m[,c(1:3,i)],id.vars = c("Ensembl.ID","GeneName","Phenotype.Male"))
  ##### Wilcoxon Rank Sum Test
  noa <- tp[tp$Phenotype.Male == "Abnormal",5]
  normal <- tp[tp$Phenotype.Male == "Normal",5]
  w <- wilcox.test(noa,normal)
  wilcoxon_list <- rbind(
    wilcoxon_list,
    data.frame(
      Features = names(sn_data.m)[i],
      PValue = w$p.value))
  ##### Plotting
  g <- ggplot(tp,aes(x = Phenotype.Male,y = value,fill = Phenotype.Male)) +
    geom_boxplot2(width = 0.8, width.errorbar = 0.5) +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = NA,colour = "Black",size = 1,linetype = 1),
      legend.position="none",
      plot.title = element_text(size = 7)) +
    labs(title = paste(names(sn_data.m)[i],w$p.value,sep = " = "))
  plist[[i]] <- g
}
#g <- ggarrange(plotlist = plist,ncol = 26,nrow = 26)
#ggsave(filename = "../../results/05.Model Construction/20220516 Feature Comparison.pdf",plot = g,width = 100,height = 100,limitsize = FALSE)
wilcoxon_list$QValue <- p.adjust(wilcoxon_list$PValue,method = "BH")
#write.csv(wilcoxon_list,"../../results/05.Model Construction/20220516 Wilcoxon Test of Features.csv",row.names = FALSE)
### Feature Selection
features_selection <- wilcoxon_list
sn_data <- sn_data.m[,names(sn_data.m) %in% c(
  "Ensembl.ID",
  "Phenotype.Male",
  wilcoxon_list[wilcoxon_list$QValue < 0.01,1])];dim(sn_data)
features_selection$Diff <- features_selection$Features %in% names(sn_data)
# Formation for Learning
sn_class <- factor(sn_data$Phenotype.Male,ordered = TRUE)
sn_rownames <- sn_data$Ensembl.ID
sn_data <- sn_data[,4:ncol(sn_data)]
row.names(sn_data) <- sn_rownames
# Zero- and Near Zero-Variance Predictors
zerovar <- nearZeroVar(sn_data)
### Feature Selection
if (!length(zerovar) == 0) {
  sn_data <- sn_data[,-zerovar];dim(sn_data)
}
features_selection$ZeroVar <- features_selection$Features %in% names(sn_data)
# Identifying Correlated Predictors
descrCorr = cor(sn_data)
highCorr = findCorrelation(descrCorr, 0.75)
### Feature Selection
sn_data <- sn_data[,-highCorr];dim(sn_data)
features_selection$HighCorr <- features_selection$Features %in% names(sn_data)
# Linear Dependencies
comboInfo <- findLinearCombos(sn_data)
### Feature Selection
if(!is.null(comboInfo$remove)){
  sn_data <- sn_data[,!comboInfo$remove];dim(sn_data)
}
features_selection$Linear <- features_selection$Features %in% names(sn_data)
# Save Features
features_selection$Employed <- features_selection$Linear
employed_feature <- features_selection[features_selection$Employed,"Features"]


# 02.Centering and Scaling -----
preProcValues <- preProcess(sn_data, method = c("center", "scale"))

# 03.Create Data Partition V1 -----
# # Select Abnormal for Validation
# pub_year <- read.csv("../../data/genelist_All.csv")
# selected_data_list <- which(rownames(sn_data) %in% pub_year[(pub_year$Published.Year >= 2020),2])
# selected_data <- sn_data[selected_data_list,]
# selected_data_class <- sn_class[selected_data_list];table(selected_data_class)
# selected_data_rowname <- rownames(selected_data)
# # Remove Selected Abnormal from Train and Test Set
# sn_data <- sn_data[-selected_data_list,]
# sn_class <- sn_class[-selected_data_list]
# # Create Train and Test Set
# set.seed(777)
# inTrain = createDataPartition(sn_class, p = 4/5, list = FALSE)
# trainx = sn_data[inTrain,];dim(trainx)
# testx = sn_data[-inTrain,];dim(testx)
# trainy = sn_class[inTrain];table(trainy)
# testy = sn_class[-inTrain];table(testy)
# # preProcess
# trainx <- predict(preProcValues, trainx)
# testx <- predict(preProcValues, testx)
# selected_data_tc <- predict(preProcValues, selected_data)

# 03.Create Data Partition V2 -----
# Create Train and Test Set
pub_year <- read.csv("../../data/genelist_All.csv")
inTrain = which(rownames(sn_data) %in% pub_year[(pub_year$Published.Year < 2020),2])
trainx = sn_data[inTrain,];dim(trainx)
testx = sn_data[-inTrain,];dim(testx)
trainy = sn_class[inTrain];table(trainy)
testy = sn_class[-inTrain];table(testy)
# preProcess
trainx <- predict(preProcValues, trainx)
testx <- predict(preProcValues, testx)

# 04.Save -----
save(
  features_selection,employed_feature,file = "Employed Features.RData")
save(
  sn_data,sn_class,preProcValues,
  trainx,trainy,testx,testy,file = "Training, Test, OutsideTest.RData")
