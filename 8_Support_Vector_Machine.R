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
svm.grid <- expand.grid(
  #C = c(seq(0.1,3,0.1),seq(3,10,1)),
  C = c(seq(1.5,3.0,0.01)),
  sigma = 0.03379925)
set.seed(777)
svmRadialFit = train(
  trainx,
  trainy, 
  method = "svmRadial",
  tuneGrid = svm.grid,
  trControl = fitControl,
  metric = "ROC",
  verbose = FALSE);svmRadialFit
best_C <- svmRadialFit$bestTune$C;best_C
#best_C <- 2.1
plot(svmRadialFit)
# sigma
svm.grid <- expand.grid(
  #sigma = seq(0,0.1,0.001),
  sigma = seq(0.0001,0.002,0.0001),
  #sigma = seq(0.0075,0.01,0.00001),
  #sigma = seq(0.0082,0.0085,0.000001),
  C = best_C)
set.seed(777)
svmRadialFit = train(
  trainx,
  trainy,
  method = "svmRadial",
  metric = "ROC",
  tuneGrid = svm.grid,
  trControl = fitControl,
  verbose = FALSE);svmRadialFit
best_sigma <- svmRadialFit$bestTune$sigma;best_sigma
#best_sigma <- 0.0007
plot(svmRadialFit)

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
svmRadial.grid <- expand.grid(
  C = best_C,
  sigma = best_sigma)
svmRadialFit = train(
  trainx,
  trainy,
  method = "svmRadial",
  metric = "ROC",
  tuneGrid = svmRadial.grid,
  trControl = fitControl,
  verbose = FALSE);svmRadialFit
# Feature Importance
dotPlot(varImp(svmRadialFit))
# Performance on Test Set
### Extract Prediction
probValues = extractProb(list(svmRadialFit),testX = testx, testY = testy)
testProbs = subset(probValues, dataType == "Test")
### Confusion Matrix
CM <- confusionMatrix(testProbs$pred, testProbs$obs);CM
fourfoldplot(CM$table, color = c("cyan", "pink"),
  conf.level = 0, margin = 1, main = "Confusion Matrix")
### Select Threshold
prediction_prob <- predict(svmRadialFit, newdata = testx, type = "prob")
roc_curve <- roc(testy, prediction_prob[,1])
best_thresh_svmRadialFit <- coords(roc = roc_curve,x = "best",input = "threshold",transpose = F, best.method = "youden") %>% .[1,1];best_thresh_svmRadialFit
### Extract Prediction with Best Threshold
probValues = extractProb(list(svmRadialFit),testX = testx, testY = testy)
probValues$pred <- ifelse(probValues$Abnormal > best_thresh_svmRadialFit,"Abnormal","Normal") %>% as.factor()
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
  labs(x = "True positive rate", y = "False positive rate",title = "Support Vector Machine") +
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
svmRadial.grid <- expand.grid(
  n.trees = best_ntree,
  interaction.depth = best_interaction_depth,
  shrinkage = best_shrinkage,
  n.minobsinnode = best_n_minobsinnode)
svmRadialFit = train(
  Class ~ .,
  imbal_all,
  method = "svmRadial",
  metric = "ROC",
  tuneGrid = svmRadial.grid,
  trControl = fitControl,
  verbose = FALSE);svmRadialFit

# 04.Save -----
save(svmRadialFit,best_thresh_svmRadialFit,file = "20220630 SVM.RData")
save(svmRadialFit,best_thresh_svmRadialFit,file = "20220703 OthN_SVM.RData")
