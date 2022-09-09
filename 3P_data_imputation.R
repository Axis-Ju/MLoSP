source("./0_presets.R")
load("./All_Merged, Unimputated, Uncalculated Features.RData")
load("./Merged, Unimputated, Uncalculated Features (Prediction).RData")

# 01.Screen -----
### Missing Pattern
to_judge <- pred_data[,c(
  1,
  4:10,
  11,38,77,
  78,98,120,134,161,170,174,179,183,
  187,203,212,242,260,273,284,362)]
missing_pattern <- md.pattern(to_judge,rotate.names = TRUE)
### Screening of Features Missing in over 80% Genes
pred_data <- pred_data[,colSums(is.na(sn_data)) <= 0.2*nrow(sn_data)]
### Screening of Genes Missing in over 20% Datasets
pred_data <- pred_data[rowSums(is.na(to_judge)) <= 0.8*ncol(to_judge),]
### Display
dim(pred_data)

# 02.Multiple Imputation of Missing Value -----
names(pred_data) <- gsub("\\s","_",names(pred_data)) ### Replace space with _ to avoid str2lang error
names(pred_data) <- gsub("_\\&_","_",names(pred_data)) ### Replace & with _ to avoid str2lang error
names(pred_data) <- gsub("#","_",names(pred_data)) ### Replace # with _ to avoid formula.character error
names(pred_data) <- gsub("-","_",names(pred_data)) ### Replace - with _ to avoid eval(predvars, data, env) error
to_imputate_data <- pred_data[,-c(1:3)]

imputation_data <- mice(
  to_imputate_data,
  m = 1,
  maxit = 5,
  seed = 777,
  method = "pmm", ### "rf","cart","norm.boot"
  remove.collinear = FALSE,
  remove.constant = FALSE)

load("./imputation_data_pmm (Pred).RData")

### Completion
imputed_data <- complete(imputation_data)
pred_data[,4:ncol(pred_data)] <- imputed_data
### Display
dim(pred_data)
table(pred_data$Phenotype.Male)

# 03.Save -----
save(pred_data,file = "Merged, Imputated, Uncalculated Features (Prediction).RData")
