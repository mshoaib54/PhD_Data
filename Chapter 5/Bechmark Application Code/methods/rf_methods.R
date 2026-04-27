rf_1 <- function(train, test){
  model <- randomForest::randomForest(answer ~ ., data = train,
                                      nodesize = 1,
                                      ntree = 100,
                                      mtry = round(sqrt(NCOL(train))) + 1)
  predict(model, newdata = test)
}

rf_2 <- function(train, test){
  model <- randomForest::randomForest(answer ~ ., data = train,
                                      nodesize = 1, ntree = 200,
                                      mtry = round(sqrt(NCOL(train))))
  predict(model, newdata = test)
}