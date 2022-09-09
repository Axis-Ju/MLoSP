source("./0_presets.R")
load("./All_Merged, Unimputated, Uncalculated Features.RData")

# 01.Screen -----
### Missing Pattern
to_judge <- sn_data[,c(
  1,
  4:10,
  11,38,77,
  78,98,120,134,161,170,174,179,183,
  187,203,212,242,260,273,284,362)]
missing_pattern <- md.pattern(to_judge,rotate.names = TRUE)
### Screening of Features Missing in over 80% Genes
sn_data <- sn_data[,colSums(is.na(sn_data)) <= 0.2*nrow(sn_data)]
### Screening of Genes Missing in over 20% Datasets
sn_data <- sn_data[rowSums(is.na(to_judge)) <= 0.8*ncol(to_judge),]
### Display
dim(sn_data)
table(sn_data$Phenotype.Male)

# 02.Multiple Imputation of Missing Value -----
names(sn_data) <- gsub("\\s","_",names(sn_data)) ### Replace space with _ to avoid str2lang error
names(sn_data) <- gsub("_\\&_","_",names(sn_data)) ### Replace & with _ to avoid str2lang error
names(sn_data) <- gsub("#","_",names(sn_data)) ### Replace # with _ to avoid formula.character error
names(sn_data) <- gsub("-","_",names(sn_data)) ### Replace - with _ to avoid eval(predvars, data, env) error
to_imputate_data <- sn_data[,-c(1:3)]
to_imputate_data$Response <- sn_data$Phenotype.Male %>% 
  gsub("Normal",0,.) %>% 
  gsub("Abnormal",1,.) %>% 
  as.numeric()

imputation_data <- mice(
  to_imputate_data,
  m = 3,
  maxit = 5,
  seed = 777,
  method = "pmm", ### "rf","cart","norm.boot"
  remove.collinear = FALSE,
  remove.constant = FALSE)

load("./imputation_data_pmm.RData")

### Visualization of Imputation
#summary(imputation_data)
#stripplot(imputation_data,col = c("grey",mdc(2)),pch = c(1,20))
### Fitting to Diagnose
imputation_data_fit <- lm.mids(Response ~ .,imputation_data)
imputation_data_pooled <- pool(imputation_data_fit)
### Diagnosis of Pooled Imputation
res.1 <- imputation_data_pooled$glanced;res.1
res.2 <- cbind(imputation_data_pooled$pooled,summary(imputation_data_pooled))
pool.r.squared(imputation_data_fit)
### Completion
imputed_data <- complete(imputation_data,action = 1)
sn_data[,4:ncol(sn_data)] <- imputed_data[,-which(names(imputed_data) == "Response")]

# 03.Save -----
save(sn_data,file = "All_Merged, Imputated, Uncalculated Features.RData")
