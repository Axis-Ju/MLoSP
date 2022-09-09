source("./0_presets.R")
load("./Employed Features.RData")
load("./Training, Test, OutsideTest.RData")
load("./20220703 OthN_RF.RData")
load("./20220703 OthN_XGBoost.RData")
load("./20220703 OthN_SVM.RData")
load("./Merged, Imputated, Calculated Features (Prediction).RData")

# 01.Extract Probability -----
pred_data <- pred_data[,c(1:2,which(colnames(pred_data) %in% employed_feature))]
### Merge with Train Data
pred_data_nrow <- nrow(pred_data)
pred_data_rowname <- c(pred_data$Ensembl.ID,rownames(sn_data))
pred_data <- predict(preProcValues,pred_data[,-c(1:2)])
pred_data <- rbind(
  pred_data,
  trainx,
  testx) %>% as.matrix()

### Extract Probability
model_to_predict <- list(rfFit,xgbFit,svmRadialFit) ##### Model to Test
prediction_exp <- predict(model_to_predict,newdata = pred_data, type = "prob")
prediction_exp <- Reduce(cbind,prediction_exp)
names(prediction_exp) <- paste(c("RF","RF","XGBoost","XGBoost","SVM","SVM"),rep(c("Abnormal","Normal"),3),sep = "_")
prediction_exp <- data.frame(
  Ensembl.ID = pred_data_rowname,
  Obs = c(
    rep("Others",sum(pred_data_nrow)),
    paste("Train",as.character(trainy),sep = "_"),
    paste("Test",as.character(testy),sep = "_")),
  prediction_exp)

# 02.Ensembl Model -----
prediction_exp$RF_Pred <- ifelse(prediction_exp$RF_Abnormal > best_thresh_rfFit,"Abnormal","Normal") %>% as.factor()
prediction_exp$XGBoost_Pred <- ifelse(prediction_exp$XGBoost_Abnormal > best_thresh_xgbFit,"Abnormal","Normal") %>% as.factor()
prediction_exp$SVM_Pred <- ifelse(prediction_exp$SVM_Abnormal > best_thresh_svmRadialFit,"Abnormal","Normal") %>% as.factor()

best_thresh_softVoting <- 0.5
prediction_exp$Soft_Voting <- ifelse(rowMeans(data.frame(prediction_exp$RF_Abnormal,prediction_exp$XGBoost_Abnormal,prediction_exp$SVM_Abnormal) %>% as.matrix()) > best_thresh_softVoting,"Abnormal","Normal") %>% as.factor()

prediction_exp$Hard_Voting<- apply(data.frame(prediction_exp$RF_Pred,prediction_exp$XGBoost_Pred,prediction_exp$SVM_Pred),1,function(v){
  n_Abnormal <- sum(v == "Abnormal")
  n_Normal <- sum(v == "Normal")
  ifelse(n_Abnormal == max(n_Abnormal,n_Normal),"Abnormal","Normal") %>% as.factor()
})

prediction_exp$Intersect <- apply(data.frame(prediction_exp$RF_Pred,prediction_exp$XGBoost_Pred,prediction_exp$SVM_Pred),1,function(v){
  n_Abnormal <- sum(v == "Abnormal")
  ifelse(n_Abnormal == 3,"Abnormal","Normal") %>% as.factor()
})

# 03.Evaluation in Train -----
prediction_exp_train <- prediction_exp[grep("Train",prediction_exp$Obs),]
prediction_exp_train_stat <- as.data.frame(array(dim = c(19,0)))
for (i in 9:14) {
  ### Confusion Matrix
  pred <- prediction_exp_train[,i]
  obs <- gsub("Train_","",prediction_exp_train$Obs) %>% as.factor()
  tmp_CM <- confusionMatrix(pred,obs)
  
  ### Train ROC
  if (i %in% 9:12) {
    obs <- ifelse(prediction_exp_train$Obs == 'Train_Abnormal',1,0)
    if (i == 12) {
      score <- rowMeans(data.frame(prediction_exp_train$RF_Abnormal,prediction_exp_train$XGBoost_Abnormal,prediction_exp_train$SVM_Abnormal) %>% as.matrix())
    } else {
      score <- prediction_exp_train[,3+2*(i-9)]
    }
    tmp_pred <- prediction(score,obs)
    tmp_auc <- performance(tmp_pred, "auc")
    tmp_auc_train <- tmp_auc@y.values[[1]]
  } else {
    tmp_auc_train <- "Can't be Calculated"
  }
  
  ### Merge
  tmp_stat <- as.data.frame(c(tmp_CM$overall,tmp_CM$byClass,AUC = tmp_auc_train))
  prediction_exp_train_stat <- cbind(prediction_exp_train_stat,tmp_stat)
  names(prediction_exp_train_stat)[ncol(prediction_exp_train_stat)] <- names(prediction_exp_train)[i]
}
prediction_exp_train_stat <- cbind(
  Statistics = rownames(prediction_exp_train_stat),
  prediction_exp_train_stat)

# 04.Evaluation in Test -----
prediction_exp_test <- prediction_exp[grep("Test",prediction_exp$Obs),]
prediction_exp_test_stat <- as.data.frame(array(dim = c(19,0)))
for (i in 9:14) {
  ### Confusion Matrix
  pred <- prediction_exp_test[,i]
  obs <- gsub("Test_","",prediction_exp_test$Obs) %>% as.factor()
  tmp_CM <- confusionMatrix(pred,obs)
  
  ### Test ROC
  if (i %in% 9:12) {
    obs <- ifelse(prediction_exp_test$Obs == 'Test_Abnormal',1,0)
    if (i == 12) {
      score <- rowMeans(data.frame(prediction_exp_test$RF_Abnormal,prediction_exp_test$XGBoost_Abnormal,prediction_exp_test$SVM_Abnormal) %>% as.matrix())
    } else {
      score <- prediction_exp_test[,3+2*(i-9)]
    }
    tmp_pred <- prediction(score,obs)
    tmp_auc <- performance(tmp_pred, "auc")
    tmp_auc_test <- tmp_auc@y.values[[1]]
  } else {
    tmp_auc_test <- "Can't be Calculated"
  }
  
  ### Merge
  tmp_stat <- as.data.frame(c(tmp_CM$overall,tmp_CM$byClass,AUC = tmp_auc_test))
  prediction_exp_test_stat <- cbind(prediction_exp_test_stat,tmp_stat)
  names(prediction_exp_test_stat)[ncol(prediction_exp_test_stat)] <- names(prediction_exp_test)[i]
}
prediction_exp_test_stat <- cbind(
  Statistics = rownames(prediction_exp_test_stat),
  prediction_exp_test_stat)

# 05.Save -----
sheets <- list(
  "Prediction Result" = prediction_exp,
  "Statistics of Train Set" = prediction_exp_train_stat,
  "Statistics of Test Set" = prediction_exp_test_stat
)
#write_xlsx(sheets,"./20220701 Prediction.xlsx",)

