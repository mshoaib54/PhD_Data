nnet_1 <- function(train, test){
  model <- nnet(answer ~ ., data = train, size = 20, 
       decay = 0.01, maxit = 5000, linout = FALSE, MaxNWts = 5000, trace = FALSE)
  predict(model, newdata = test)
}