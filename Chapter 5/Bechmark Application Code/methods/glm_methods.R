glm_binomial <- function(train, test, optimizecutoff = FALSE){
  suppressWarnings(model <- glm(answer ~ ., data = train, family = binomial))
  suppressWarnings(prob <- predict(model, test, type = "response"))
  factor(ifelse(prob[, 2] >= cutoff, levels(train$answer)[2], 
                        levels(train$answer)[1]), levels = levels(train$answer))
}
