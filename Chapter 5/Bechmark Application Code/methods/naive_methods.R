naive_bayes_1 <- function(train, test) {
  model <- e1071::naiveBayes(answer ~ ., data = train)
  predict(model, newdata = test)
}



