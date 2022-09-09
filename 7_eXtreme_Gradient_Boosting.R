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
# nrounds
xgb.grid <- expand.grid(
  #nrounds = c(100,300,500,700,1000,1500), ### Boosting Iterations
  nrounds = seq(10,150,10), ### Boosting Iterations
  max_depth = 5, ### Max Tree Depth
  eta = 0.3, ### Shrinkage
  gamma = 0, ### Minimum Loss Reduction
  subsample = 1, ### Subsample Percentage
  colsample_bytree = 1, ### Subsample Ratio of Columns
  min_child_weight = 1 ### Minimum Sum of Instance Weight
)
set.seed(777)
xgbFit = train(
  trainx,
  trainy,
  method = "xgbTree",
  metric = "ROC",
  tuneGrid = xgb.grid,
  trControl = fitControl,
  verbose = FALSE);xgbFit
best_nrounds <- xgbFit$bestTune$nrounds;best_nrounds
#best_nrounds <- 90
plot(xgbFit)
# max_depth
xgb.grid <- expand.grid(
  nrounds = best_nrounds, ### Boosting Iterations
  max_depth = seq(1,15,1), ### Max Tree Depth
  eta = 0.3, ### Shrinkage
  gamma = 0, ### Minimum Loss Reduction
  subsample = 1, ### Subsample Percentage
  colsample_bytree = 1, ### Subsample Ratio of Columns
  min_child_weight = 1 ### Minimum Sum of Instance Weight
)
set.seed(777)
xgbFit = train(
  trainx,
  trainy,
  method = "xgbTree",
  metric = "ROC",
  tuneGrid = xgb.grid,
  trControl = fitControl,
  verbose = FALSE);xgbFit
best_max_depth <- xgbFit$bestTune$max_depth;best_max_depth
#best_max_depth <- 6
plot(xgbFit)
# eta
xgb.grid <- expand.grid(
  nrounds = best_nrounds, ### Boosting Iterations
  max_depth = best_max_depth, ### Max Tree Depth
  eta = c(0.0001,0.001,seq(0.01,0.5,0.01)), ### Shrinkage
  gamma = 0, ### Minimum Loss Reduction
  subsample = 1, ### Subsample Percentage
  colsample_bytree = 1, ### Subsample Ratio of Columns
  min_child_weight = 1 ### Minimum Sum of Instance Weight
)
set.seed(777)
xgbFit = train(
  trainx,
  trainy,
  method = "xgbTree",
  metric = "ROC",
  tuneGrid = xgb.grid,
  trControl = fitControl,
  verbose = FALSE);xgbFit
best_eta <- xgbFit$bestTune$eta;best_eta
#best_eta <- 0.26
plot(xgbFit)
# gamma
xgb.grid <- expand.grid(
  nrounds = best_nrounds, ### Boosting Iterations
  max_depth = best_max_depth, ### Max Tree Depth
  eta = best_eta, ### Shrinkage
  gamma = seq(0,1,0.1), ### Minimum Loss Reduction
  subsample = 1, ### Subsample Percentage
  colsample_bytree = 1, ### Subsample Ratio of Columns
  min_child_weight = 1 ### Minimum Sum of Instance Weight
)
set.seed(777)
xgbFit = train(
  trainx,
  trainy,
  method = "xgbTree",
  metric = "ROC",
  tuneGrid = xgb.grid,
  trControl = fitControl,
  verbose = FALSE);xgbFit
best_gamma <- xgbFit$bestTune$gamma;best_gamma
#best_gamma <- 0
plot(xgbFit)
# subsample
xgb.grid <- expand.grid(
  nrounds = best_nrounds, ### Boosting Iterations
  max_depth = best_max_depth, ### Max Tree Depth
  eta = best_eta, ### Shrinkage
  gamma = best_gamma, ### Minimum Loss Reduction
  subsample = seq(0,1,0.1), ### Subsample Percentage
  colsample_bytree = 1, ### Subsample Ratio of Columns
  min_child_weight = 1 ### Minimum Sum of Instance Weight
)
set.seed(777)
xgbFit = train(
  trainx,
  trainy,
  method = "xgbTree",
  metric = "ROC",
  tuneGrid = xgb.grid,
  trControl = fitControl,
  verbose = FALSE);xgbFit
best_subsample <- xgbFit$bestTune$subsample;best_subsample
#best_subsample <- 1
plot(xgbFit)
# colsample_bytree
xgb.grid <- expand.grid(
  nrounds = best_nrounds, ### Boosting Iterations
  max_depth = best_max_depth, ### Max Tree Depth
  eta = best_eta, ### Shrinkage
  gamma = best_gamma, ### Minimum Loss Reduction
  subsample = best_subsample, ### Subsample Percentage
  colsample_bytree = seq(0,1,0.1), ### Subsample Ratio of Columns
  min_child_weight = 1 ### Minimum Sum of Instance Weight
)
set.seed(777)
xgbFit = train(
  trainx,
  trainy,
  method = "xgbTree",
  metric = "ROC",
  tuneGrid = xgb.grid,
  trControl = fitControl,
  verbose = FALSE);xgbFit
best_colsample_bytree <- xgbFit$bestTune$colsample_bytree;best_colsample_bytree
#best_colsample_bytree <- 0.2
plot(xgbFit)
# min_child_weight
xgb.grid <- expand.grid(
  nrounds = best_nrounds, ### Boosting Iterations
  max_depth = best_max_depth, ### Max Tree Depth
  eta = best_eta, ### Shrinkage
  gamma = best_gamma, ### Minimum Loss Reduction
  subsample = best_subsample, ### Subsample Percentage
  colsample_bytree = best_colsample_bytree, ### Subsample Ratio of Columns
  min_child_weight = c(seq(0,1,0.1),3,5,7,9) ### Minimum Sum of Instance Weight
)
set.seed(777)
xgbFit = train(
  trainx,
  trainy,
  method = "xgbTree",
  metric = "ROC",
  tuneGrid = xgb.grid,
  trControl = fitControl,
  verbose = FALSE);xgbFit
best_min_child_weight <- xgbFit$bestTune$min_child_weight;best_min_child_weight
#best_min_child_weight <- 0.7
plot(xgbFit)

# 02.Best Model -----
# Set Train Control
set.seed(777)
fitControl = trainControl(
  method = "boot", ### boot, LOOCV
  number = 30,
  #repeats = 10,
  classProbs = TRUE,
  returnData = TRUE,
  verboseIter = TRUE,
  savePredictions = TRUE,
  summaryFunction = twoClassSummary ### twoClassSummary, mnLogLoss, prSummary
)
set.seed(777)
xgb.grid <- expand.grid(
  nrounds = best_nrounds, ### Boosting Iterations
  max_depth = best_max_depth, ### Max Tree Depth
  eta = best_eta, ### Shrinkage
  gamma = best_gamma, ### Minimum Loss Reduction
  subsample = best_subsample, ### Subsample Percentage
  colsample_bytree = best_colsample_bytree, ### Subsample Ratio of Columns
  min_child_weight = best_min_child_weight ### Minimum Sum of Instance Weight
)
xgbFit = train(
  trainx,
  trainy,
  method = "xgbTree",
  metric = "ROC",
  tuneGrid = xgb.grid,
  trControl = fitControl,
  verbose = FALSE);xgbFit
# Feature Importance
dotPlot(varImp(xgbFit))
# Performance on Test Set
### Extract Prediction
probValues = extractProb(list(xgbFit),testX = testx, testY = testy)
testProbs = subset(probValues, dataType == "Test")
### Confusion Matrix
CM <- confusionMatrix(testProbs$pred, testProbs$obs);CM
fourfoldplot(CM$table, color = c("cyan", "pink"),
  conf.level = 0, margin = 1, main = "Confusion Matrix")
### Select Threshold
prediction_prob <- predict(xgbFit, newdata = testx, type = "prob")
roc_curve <- roc(testy, prediction_prob[,1])
best_thresh_xgbFit <- coords(roc = roc_curve,x = "best",input = "threshold",transpose = F, best.method = "youden") %>% .[1,1];best_thresh_xgbFit
### Extract Prediction with Best Threshold
probValues = extractProb(list(xgbFit),testX = testx, testY = testy)
probValues$pred <- ifelse(probValues$Abnormal > best_thresh_xgbFit,"Abnormal","Normal") %>% as.factor()
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

# 03.Final Model -----
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
xgb.grid <- expand.grid(
  n.trees = best_ntree,
  interaction.depth = best_interaction_depth,
  shrinkage = best_shrinkage,
  n.minobsinnode = best_n_minobsinnode)
xgbFit = train(
  Class ~ .,
  imbal_all,
  method = "xgb",
  metric = "ROC",
  tuneGrid = xgb.grid,
  trControl = fitControl,
  verbose = FALSE);xgbFit

# 04.Save -----
save(xgbFit,best_thresh_xgbFit,file = "20220630 XGBoost.RData")
save(xgbFit,best_thresh_xgbFit,file = "20220703 OthN_XGBoost.RData")
