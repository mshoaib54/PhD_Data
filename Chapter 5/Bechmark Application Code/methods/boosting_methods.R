boosting_basic <- function(train, test) {
  model <- adabag::boosting(answer ~ ., data = train, mfinal = 20)
  predict(model, train)$class
}
