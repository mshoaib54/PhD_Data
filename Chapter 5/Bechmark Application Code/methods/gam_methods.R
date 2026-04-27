gam_basic <- function(train, test, optimizecutoff = FALSE, ...) {
  suppressWarnings(model <- gam::gam(answer ~ ., data = train, family = binomial))
  suppressWarnings(prob <- predict(model, test, type = "response"))
  factor(ifelse(prob >= 0.5, levels(train$answer)[1], 
                levels(train$answer)[2]), levels = levels(train$answer))
}
