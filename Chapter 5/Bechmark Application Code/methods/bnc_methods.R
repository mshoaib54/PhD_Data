bnc_nb <- function(train, test){
  model <- bnclassify::nb(class = "answer", dataset = train)
  model <- bnclassify::lp(model, dataset = train, smooth = 1)
  predict(model, test)
}

bnc_tan_cl <- function(train, test){
  model <- bnclassify::tan_cl(class = "answer", dataset = train)
  model <- bnclassify::lp(model, dataset = train, smooth = 1)
  predict(model, test)
}

bnc_tan_hc <- function(train, test){
  model <- bnclassify::tan_hc(class = "answer", dataset = train, k = 5)
  model <- bnclassify::lp(model, dataset = train, smooth = 1)
  predict(model, test)
}

bnc_fssj <- function(train, test){
  model <- bnclassify::fssj(class = "answer", dataset = train, 
                            k = 5, smooth = 1)
  model <- bnclassify::lp(model, dataset = train, smooth = 1)
  predict(model, test)
}

bnc_bsej <- function(train, test){
  model <- bnclassify::bsej(class = "answer", dataset = train, 
                            k = 5, smooth = 1)
  model <- bnclassify::lp(model, dataset = train, smooth = 1)
  predict(model, test)
}


bnc_3db <- function(train, test){
  kdbk <- min(3, ncol(train)-1)
  model <- bnclassify::kdb(class = "answer", dataset = train, 
                           k = 5, smooth = 1, kdbk = kdbk)
  model <- bnclassify::lp(model, dataset = train, smooth = 1)
  predict(model, test)
}

bnc_5db <- function(train, test){
  kdbk <- min(5, ncol(train)-1)
  model <- bnclassify::kdb(class = "answer", dataset = train, 
                           k = 5, smooth = 1, kdbk = kdbk)
  model <- bnclassify::lp(model, dataset = train, smooth = 1)
  predict(model, test)
}