source("./0_presets.R")
load("./Employed Features.RData")
load("./Training, Test, OutsideTest.RData")

# 01.Parameter Tuning -----
# Set Train Control
set.seed(777)
fitControl = trainControl(
  method = "boot", ### boot, LOOCV
  number = 30,
  #repeats = 10,
  classProbs = TRUE,
  returnData = TRUE,
  verboseIter = TRUE,
  summaryFunction = twoClassSummary ### twoClassSummary, mnLogLoss, prSummary
)
# n.trees
gbm.grid <- expand.grid(
  #n.trees = seq(50,500,50),
  n.trees = seq(5,100,1),
  #n.trees = seq(15,30,1),
  interaction.depth = 1,
  shrinkage = 0.1,
  n.minobsinnode = 10)
set.seed(777)
gbmFit = train(
  trainx,
  trainy,
  method = "gbm",
  metric = "ROC",
  tuneGrid = gbm.grid,
  trControl = fitControl,
  verbose = FALSE);gbmFit
best_ntree <- gbmFit$bestTune$n.trees;best_ntree
#best_ntree <- 77
plot(gbmFit)
# interaction.depth
gbm.grid <- expand.grid(
  n.trees = best_ntree,
  interaction.depth = seq(1,50,1),
  shrinkage = 0.1,
  n.minobsinnode = 10)
set.seed(777)
gbmFit = train(
  trainx,
  trainy,
  method = "gbm",
  metric = "ROC",
  tuneGrid = gbm.grid,
  trControl = fitControl,
  verbose = FALSE);gbmFit
best_interaction_depth <- gbmFit$bestTune$interaction.depth;best_interaction_depth
#best_interaction_depth <- 3
plot(gbmFit)
# shrinkage
gbm.grid <- expand.grid(
  n.trees = best_ntree,
  interaction.depth = best_interaction_depth,
  #shrinkage = seq(0,1,0.1),
  #shrinkage = seq(0,0.3,0.01),
  #shrinkage = seq(0,0.07,0.0001),
  n.minobsinnode = 10)
set.seed(777)
gbmFit = train(
  trainx,
  trainy,
  method = "gbm",
  metric = "ROC",
  tuneGrid = gbm.grid,
  trControl = fitControl,
  verbose = FALSE);gbmFit
best_shrinkage <- gbmFit$bestTune$shrinkage;best_shrinkage
#best_shrinkage <- 0.04
plot(gbmFit)
# n.minobsinnode
gbm.grid <- expand.grid(
  n.trees = best_ntree,
  interaction.depth = best_interaction_depth,
  shrinkage = best_shrinkage,
  n.minobsinnode = seq(1,100,1)
  #n.minobsinnode = seq(1,30,1)
)
set.seed(777)
gbmFit = train(
  trainx,
  trainy,
  method = "gbm",
  metric = "ROC",
  tuneGrid = gbm.grid,
  trControl = fitControl,
  verbose = FALSE);gbmFit
best_n_minobsinnode <- gbmFit$bestTune$n.minobsinnode;best_n_minobsinnode
#best_n_minobsinnode <- 37
plot(gbmFit)

# 02.Best Model -----
# Set Train Control
set.seed(777)
fitControl = trainControl(
  method = "boot", ### boot, LOOCV
  #number = 5,
  #repeats = 10,
  classProbs = TRUE,
  returnData = TRUE,
  verboseIter = TRUE,
  summaryFunction = twoClassSummary ### twoClassSummary, mnLogLoss, prSummary
)
set.seed(777)
gbm.grid <- expand.grid(
  n.trees = best_ntree,
  interaction.depth = best_interaction_depth,
  shrinkage = best_shrinkage,
  n.minobsinnode = best_n_minobsinnode)
gbmFit = train(
  trainx,
  trainy,
  method = "gbm",
  metric = "ROC",
  tuneGrid = gbm.grid,
  trControl = fitControl,
  verbose = FALSE);gbmFit
# Feature Importance
dotPlot(varImp(gbmFit))
# Performance on Test Set
### Extract Prediction
probValues = extractProb(list(gbmFit),testX = testx, testY = testy)
testProbs = subset(probValues, dataType == "Test")
### Confusion Matrix
CM <- confusionMatrix(testProbs$pred, testProbs$obs);CM
fourfoldplot(CM$table, color = c("cyan", "pink"),
  conf.level = 0, margin = 1, main = "Confusion Matrix")
### Select Threshold
prediction_prob <- predict(gbmFit, newdata = testx, type = "prob")
roc_curve <- roc(testy, prediction_prob[,1])
best_thresh_gbmFit <- coords(roc = roc_curve,x = "best",input = "threshold",transpose = F, best.method = "youden") %>% .[1,1];best_thresh_gbmFit
### Extract Prediction with Best Threshold
probValues = extractProb(list(gbmFit),testX = testx, testY = testy)
probValues$pred <- ifelse(probValues$Abnormal > best_thresh_gbmFit,"Abnormal","Normal") %>% as.factor()
testProbs = subset(probValues, dataType == "Test")
### Confusion Matrix with Best Threshold
CM <- confusionMatrix(testProbs$pred, testProbs$obs);CM
fourfoldplot(CM$table, color = c("cyan", "pink"),
  conf.level = 0, margin = 1, main = "Confusion Matrix")
### ROC Plot with Best Threshold
testProbs$lable <- ifelse(testProbs$obs=='Abnormal',yes=1,0) ##### abnormal = true
testPreds <- prediction(testProbs$Abnormal,testProbs$lable)
testPerf <- performance(testPreds, measure="tpr", x.measure="fpr" )
auc <- performance(testPreds, "auc")
auc@y.values[[1]]
### Plotting ROC Curve
tibble::tibble(x = testPerf@x.values[[1]], y = testPerf@y.values[[1]]) %>%
  ggplot(aes(x = x, y = y)) + 
  geom_line() +
  geom_point() +
  theme_bw()+
  labs(x = "True positive rate", y = "False positive rate",title = "Stochastic Gradient Boosting") +
  scale_color_brewer(palette = "Reds") + 
  theme(legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    axis.title = element_text(face = "bold")) +
  geom_text(aes(x=1, y=0, label=paste0("AUC: ",auc@y.values[[1]] %>% round(2))),hjust="right",vjust="bottom")
# Performance on Selected Data
probValues = extractProb(list(gbmFit),testX = selected_data, testY = selected_data_class)
probValues$pred <- ifelse(probValues$Abnormal > best_thresh_gbmFit,"Abnormal","Normal") %>% as.factor()
testProbs = subset(probValues, dataType == "Test")
### Confusion Matrix with Best Threshold
CM <- confusionMatrix(testProbs$pred, testProbs$obs);CM
fourfoldplot(CM$table, color = c("cyan", "pink"),
  conf.level = 0, margin = 1, main = "Confusion Matrix")
### ROC Plot with Best Threshold
testProbs$lable <- ifelse(testProbs$obs=='Abnormal',yes=1,0) ##### abnormal = true
testPreds <- prediction(testProbs$Abnormal,testProbs$lable)
testPerf <- performance(testPreds, measure="tpr", x.measure="fpr")
auc <- performance(testPreds, "auc")
auc@y.values[[1]]
### Plotting ROC Curve
tibble::tibble(x = testPerf@x.values[[1]], y = testPerf@y.values[[1]]) %>%
  ggplot(aes(x = x, y = y)) + 
  geom_line() +
  geom_point() +
  theme_bw()+
  labs(x = "True positive rate", y = "False positive rate",title = "Stochastic Gradient Boosting") +
  scale_color_brewer(palette = "Reds") + 
  theme(legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    axis.title = element_text(face = "bold")) +
  geom_text(aes(x=1, y=0, label=paste0("AUC: ",auc@y.values[[1]] %>% round(2))),hjust="right",vjust="bottom")

# 03.Predictions -----
# Predicting with Train Data (Any Celltype Exp)
load("Merged, Imputated, Calculated Features (Prediction).RData")
pred_data <- pred_data[,c(1:2,which(colnames(pred_data) %in% employed_feature))]
### Merge with Train Data
pred_data_rowname <- c(pred_data$Ensembl.ID,rownames(sn_data))
pred_data <- rbind(pred_data[,-c(1:2)],sn_data) %>% as.matrix()
pred_data <- predict(preProcValues, pred_data)
### Extract Probability
model_to_predict <- gbmFit ##### Model to Test
prediction_exp_withTrain <- predict(model_to_predict,newdata = pred_data, type = "prob")
prediction_exp_withTrain <- data.frame(
  Ensembl.ID = pred_data_rowname,
  Obs = c(rep("Others",sum(max_exp_all)),as.character(sn_class)),
  prediction_exp_withTrain)
#write.csv(prediction_exp_withTrain,"Prediction of Selected Abnormal + Lab Normal + All Genes.csv",row.names = FALSE)
### Train&Test (Abnormal) in Y
x <- data.frame(Ensembl.ID = rownames(sn_data))[sn_class == "Abnormal",] ##### Abnormal
y <- prediction_exp_withTrain[order(prediction_exp_withTrain$Abnormal,decreasing = TRUE),1] ##### All
distribution <- sort(match(x,y));distribution
### Train&Test (Normal) in Y
x <- data.frame(Ensembl.ID = rownames(sn_data))[sn_class == "Normal",] ##### Normal
y <- prediction_exp_withTrain[order(prediction_exp_withTrain$Abnormal,decreasing = TRUE),1] ##### All
distribution <- sort(match(x,y));distribution

# Predicting with Selected Data (Any Celltype Exp)
load("Merged, Imputated, Calculated Features (Prediction).RData")
pred_data <- pred_data[,c(1:2,which(colnames(pred_data) %in% employed_feature))]
### Merge with Selected Data
pred_data_rowname <- c(pred_data$Ensembl.ID,selected_data_rowname)
pred_data <- rbind(pred_data[,-c(1:2)],selected_data) %>% as.matrix()
pred_data <- predict(preProcValues, pred_data)
### Extract Probability
model_to_predict <- gbmFit ##### Model to Test
prediction_exp_withSelected <- predict(model_to_predict,newdata = pred_data, type = "prob")
prediction_exp_withSelected <- data.frame(
  Ensembl.ID = pred_data_rowname,
  Obs = c(rep("Others",sum(max_exp_all)),as.character(selected_data_class)),
  prediction_exp_withSelected)
#write.csv(prediction_exp_withSelected,"Prediction of Selected Abnormal + Lab Normal + All Genes.csv",row.names = FALSE)
### Selected (Abnormal) in Y
x <- data.frame(Ensembl.ID = rownames(selected_data))[selected_data_class == "Abnormal",1] ##### Abnormal
y <- prediction_exp_withSelected[order(prediction_exp_withSelected$Abnormal,decreasing = TRUE),1] ##### All
distribution <- sort(match(x,y));distribution
### Selected (Normal) in Y
x <- data.frame(Ensembl.ID = rownames(selected_data))[selected_data_class == "Normal",1] ##### Normal
y <- prediction_exp_withSelected[order(prediction_exp_withSelected$Abnormal,decreasing = TRUE),1] ##### All
distribution <- sort(match(x,y));distribution

# 04.Final Model -----
imbal_all <- data.frame(
  Class = factor(ifelse(sn_class == "Abnormal",TRUE,FALSE)),
  sn_data)
levels(imbal_all$Class) = c("Normal","Abnormal")
set.seed(777)
fitControl = trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  returnData = TRUE,
  savePredictions = "final",
  summaryFunction = twoClassSummary,
  index = createResample(imbal_all$Class,5),
  classProbs = TRUE,
  sampling = "down")
set.seed(777)
gbm.grid <- expand.grid(
  n.trees = best_ntree,
  interaction.depth = best_interaction_depth,
  shrinkage = best_shrinkage,
  n.minobsinnode = best_n_minobsinnode)
gbmFit = train(
  Class ~ .,
  imbal_all,
  method = "gbm",
  metric = "ROC",
  tuneGrid = gbm.grid,
  trControl = fitControl,
  verbose = FALSE);gbmFit

# 05.Save -----
save(gbmFit,best_thresh_gbmFit,file = "GBM.RData")
