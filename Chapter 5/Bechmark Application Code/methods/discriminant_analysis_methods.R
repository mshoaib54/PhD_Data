regularized_da <- function(train, test, optimizecutoff = FALSE, ...) {
  model <- klaR::rda(answer ~ ., data = train)
  predict(model, test, aslist = FALSE)
}
