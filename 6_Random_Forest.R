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
# mtry
rf.grid <- expand.grid(
  mtry = seq(2,sqrt(ncol(sn_data))+20,1)
  )
set.seed(777)
rfFit = train(
  trainx,
  trainy,
  method = "rf",
  metric = "ROC", ### Accuracy, Kappa
  tuneGrid = rf.grid,
  trControl = fitControl);rfFit
best_mtry <- rfFit$bestTune$mtry;best_mtry
#best_mtry <- 7
plot(rfFit)
# maxnodes
rf.grid <- expand.grid(mtry = best_mtry)
maxnodes_tuning <- list()
for (i in seq(50,ncol(sn_data),5)) {
  print(i)
  set.seed(777)
  rfFit = train(
    trainx,
    trainy,
    method = "rf",
    maxnodes = i,
    metric = "ROC",
    tuneGrid = rf.grid,
    trControl = fitControl)
  current_iteration <- toString(i)
  maxnodes_tuning[[current_iteration]] <- rfFit
}
result_maxnodes_tuning <- resamples(maxnodes_tuning)
summary(result_maxnodes_tuning)
bwplot(result_maxnodes_tuning)
best_maxnodes <- 75
# nodesize
rf.grid <- expand.grid(mtry = best_mtry)
nodesize_tuning <- list()
for (i in seq(1,10,1)) {
  print(i)
  set.seed(777)
  rfFit = train(
    trainx,
    trainy,
    method = "rf",
    maxnodes = best_maxnodes,
    nodesize = i,
    metric = "ROC",
    tuneGrid = rf.grid,
    trControl = fitControl)
  current_iteration <- toString(i)
  nodesize_tuning[[current_iteration]] <- rfFit
}
result_nodesize_tuning <- resamples(nodesize_tuning)
summary(result_nodesize_tuning)
bwplot(result_nodesize_tuning)
best_nodesize <- 1
# ntree
rf.grid <- expand.grid(mtry = best_mtry)
ntree_tuning <- list()
for (i in c(100,250,300,350,400,450,500,550,600,800,1000,2000,3000)) {
  print(i)
  set.seed(777)
  rfFit = train(
    trainx,
    trainy,
    method = "rf",
    maxnodes = best_maxnodes,
    #nodesize = best_nodesize,
    ntree = i,
    metric = "ROC",
    tuneGrid = rf.grid,
    trControl = fitControl)
  current_iteration <- toString(i)
  ntree_tuning[[current_iteration]] <- rfFit
}
result_ntree_tuning <- resamples(ntree_tuning)
summary(result_ntree_tuning)
bwplot(result_ntree_tuning)
best_ntree <- 100

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
rf.grid <- expand.grid(mtry = best_mtry)
set.seed(777)
rfFit = train(
  trainx,
  trainy,
  method = "rf",
  #maxnodes = best_maxnodes,
  #nodesize = best_nodesize,
  #ntree = best_ntree,
  metric = "ROC",
  tuneGrid = rf.grid,
  trControl = fitControl);rfFit
# Feature Importance
dotPlot(varImp(rfFit))
# Performance on Test Set
### Extract Prediction
probValues = extractProb(list(rfFit),testX = testx, testY = testy)
testProbs = subset(probValues, dataType == "Test")
### Confusion Matrix
CM <- confusionMatrix(testProbs$pred, testProbs$obs);CM
fourfoldplot(CM$table, color = c("cyan", "pink"),
  conf.level = 0, margin = 1, main = "Confusion Matrix")
### Select Threshold
prediction_prob <- predict(rfFit, newdata = testx, type = "prob")
roc_curve <- roc(testy, prediction_prob[,1])
best_thresh_rfFit <- coords(roc = roc_curve,x = "best",input = "threshold",transpose = F, best.method = "youden") %>% .[1,1];best_thresh_rfFit
### Extract Prediction with Best Threshold
probValues = extractProb(list(rfFit),testX = testx, testY = testy)
probValues$pred <- ifelse(probValues$Abnormal > best_thresh_rfFit,"Abnormal","Normal") %>% as.factor()
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
  labs(x = "True positive rate", y = "False positive rate",title = "Random Forest") +
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
rf.grid <- expand.grid(mtry = best_mtry)
rfFit = train(
  Class ~ .,
  imbal_all,
  method = "rf",
  #maxnodes = best_maxnodes,
  #nodesize = best_nodesize,
  #ntree = best_ntree,
  metric = "ROC",
  tuneGrid = rf.grid,
  trControl = fitControl);rfFit

# 04.Save -----
save(rfFit,best_thresh_rfFit,file = "20220630 RF.RData")
save(rfFit,best_thresh_rfFit,file = "20220703 OthN_RF.RData")
