predict_bagging <- function(model, train, test, optimizecutoff = FALSE){
  if (optimizecutoff){
    prob <- predict(model, newdata = train, type = "prob")
    cutoff <- InformationValue::optimalCutoff(actuals = as.numeric(train$answer) - 1, 
                                              predictedScores = prob[, 2], 
                                              optimiseFor = "Both")
    prob <- predict(model, newdata = test, type = "prob")
    pred <- factor(ifelse(prob[, 2] >= cutoff, levels(train$answer)[2], 
                  levels(train$answer)[1]), levels = levels(train$answer))
    prob <- prob[, 2]
  }
  else{
    cutoff <- 0.5
    prob <- predict(model, newdata = test, type = "prob")
    pred <- factor(ifelse(prob[, 2] >= cutoff, levels(train$answer)[2], 
                          levels(train$answer)[1]), levels = levels(train$answer))
    prob <- prob[, 2]
  }
  return(list(pred = pred, prob = prob, cutoff = cutoff))
}


bagging_basic <- function(train, test, optimizecutoff = FALSE, ...) {
  model <- ipred::bagging(answer ~ ., data = train)
  predict_bagging(model, train, test, optimizecutoff)
}
