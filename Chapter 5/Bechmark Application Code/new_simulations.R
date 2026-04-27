library("stagedtrees")
source("methods.R")
source("statistics.R")

nreps <- 100
n <- 1000
ps <- 2:15
k <- 2
stats <- c("time", "accuracy", "f1")

torun <- c("st_naive_totvar_cmi", 
             "st_naive_jefferys_cmi",
             "st_naive_kaniadakis_cmi",
             "st_naive_jensen_shannon_cmi",
             "st_naive_hellinger_cmi",
           "rf_1")

TABLE <- array(
  data = NA,
  dim = c(
    length(stats),
    length(torun),
    length(ps),
    1,
    nreps
  ),
  dimnames = list(
    stat = stats,
    classifier = torun,
    p = ps,
    k = k,
    rep = 1:nreps
  )
)

for (p in ps){
  tree <- lapply(1:p, function(i) paste0(1:k))
  names(tree) <- paste0("X",1:p)
  tree <- c(list(answer = c("0","1")), tree)
  pb <- txtProgressBar(min = 1, max = nreps, style = 3)
  message(p)
  for (rep in 1:nreps){
    setTxtProgressBar(pb, rep)
    todo <- TRUE
    while (todo) {
      model <- random_sevt(x = tree, q = 0.5)
      train <- sample_from(model, size =  n)
      test <- sample_from(model, size =  n)
      if (length(unique(train$answer)) == 2){
        todo <- FALSE
      }
    }
    for (c_name in torun){
      #message(c_name)
      c_fun <- get(c_name)
      time_ <- system.time(predict <- c_fun(train, test))[3]
      res <- list(time = time_, predict = predict)
      for (stat in stats){
        stat_fun <- get(stat)
        TABLE[stat, c_name, paste0(p), paste0(k), rep] <- stat_fun(res, test$answer)
      }
    }
  }
  close(pb)
}


saveRDS(TABLE, "TABLE_NEW_SIMULATION.rds")


TABLE <- array(
  data = NA,
  dim = c(
    length(stats),
    length(torun),
    length(ps),
    1,
    nreps
  ),
  dimnames = list(
    stat = stats,
    classifier = torun,
    p = ps,
    k = k,
    rep = 1:nreps
  )
)

  for (p in ps){
    pb <- txtProgressBar(min = 1, max = nreps, style = 3)
    message(p)
    for (rep in 1:nreps){
      setTxtProgressBar(pb, rep)
      data <- generate_linear_dataset(p, n = 2 * n, eps = 0, gamma = 0, 
                                      alpha = runif(p, min = -1, max = 1))
      colnames(data)[1] <- "answer"
      train <- data[1:n,]
      test <- data[-c(1:n),]
      for (c_name in torun){
        #message(c_name)
        c_fun <- get(c_name)
        time_ <- system.time(predict <- c_fun(train, test))[3]
        res <- list(time = time_, predict = predict)
        for (stat in stats){
          stat_fun <- get(stat)
          TABLE[stat, c_name, paste0(p), paste0(k), rep] <- stat_fun(res, test$answer)
        }
      }
    }
    close(pb)
  }

saveRDS(TABLE, "TABLE_NEW_SIMULATION_LINEAR.rds")

TABLE <- array(
  data = NA,
  dim = c(
    length(stats),
    length(torun),
    length(ps),
    1,
    nreps
  ),
  dimnames = list(
    stat = stats,
    classifier = torun,
    p = ps,
    k = k,
    rep = 1:nreps
  )
)

  for (p in ps){
    pb <- txtProgressBar(min = 1, max = nreps, style = 3)
    message(p)
    for (rep in 1:nreps){
      setTxtProgressBar(pb, rep)
      data <- generate_xor_dataset(p, n = 2 * n, eps = 1.2)
      colnames(data)[1] <- "answer"
      train <- data[1:n,]
      test <- data[-c(1:n),]
      for (c_name in torun){
        #message(c_name)
        c_fun <- get(c_name)
        time_ <- system.time(predict <- c_fun(train, test))[3]
        res <- list(time = time_, predict = predict)
        for (stat in stats){
          stat_fun <- get(stat)
          TABLE[stat, c_name, paste0(p), paste0(k), rep] <- stat_fun(res, test$answer)
        }
      }
    }
    close(pb)
  }

saveRDS(TABLE, "TABLE_NEW_SIMULATION_XOR.rds")
