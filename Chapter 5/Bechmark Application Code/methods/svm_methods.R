predict_svm <- function(model, train, test, optimizecutoff){
  if (optimizecutoff){
    prob <- attr(predict(model, newdata = train, probability = T), "probabilities")
    prob <- prob[, attr(prob, "dimnames")[[2]] == levels(train$answer)[2]]
    cutoff <- InformationValue::optimalCutoff(actuals = as.numeric(train$answer) - 1, 
                                              predictedScores = prob, 
                                              optimiseFor = "Both")
    prob <- attr(predict(model, newdata = test, probability = T), "probabilities")
    prob <- prob[, attr(prob, "dimnames")[[2]] == levels(train$answer)[2]]
    pred <- factor(ifelse(prob >= cutoff, levels(train$answer)[2], 
                  levels(train$answer)[1]), levels = levels(train$answer))
  }
  else{
    cutoff <- 0.5
    prob <- attr(predict(model, newdata = test, probability = T), "probabilities")
    prob <- prob[, attr(prob, "dimnames")[[2]] == levels(train$answer)[2]]
    pred <- factor(ifelse(prob >= cutoff, levels(train$answer)[2], 
                          levels(train$answer)[1]), levels = levels(train$answer))
  }
  return(list(pred = pred, prob = prob, cutoff = cutoff))
}


svm_basic <- function(train, test, optimizecutoff = FALSE, ...) {
  model <- e1071::svm(answer ~ ., data = train, probability = TRUE)
  predict_svm(model, train, test, optimizecutoff)
}
