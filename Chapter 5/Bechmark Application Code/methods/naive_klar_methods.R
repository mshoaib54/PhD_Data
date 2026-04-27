naive_bayes_klar <- function(train, test) {
  model <- klaR::NaiveBayes(answer ~ ., data = train)
  predict(model, newdata = test)
}
