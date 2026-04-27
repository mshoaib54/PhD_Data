args <- commandArgs(trailingOnly = TRUE)
datasets <- factor(read.table("binary_datasets_names_final.tsv")[, 1])
source("methods.R")
library("cli")

if (length(args) > 0){
  if (base::endsWith(args[1], ".tsv")){ 
    datasets <- read.table(args[1])[,1]    
  }else{ ## it is the name of a dataset
    datasets <- args[1]
  }
}
if (length(args) > 1){
  torun <- unlist(sapply(args[-1], function(x){
    if (base::endsWith(x, "_")){
      classifiers[base::startsWith(classifiers, x)]
    }else{
      x
    }
  }))
}else{
  torun <- classifiers
}

nreps <- 10

dir.create("results/", showWarnings = FALSE)

set.seed(2020)
for (d in datasets){
  message(d)
  timestamp()
  res_path <- paste0("results/",d,"/") 
  dir.create(res_path, showWarnings = FALSE)
  data <- readRDS(paste0("datasets/", d, ".rds"))
  split_path <- paste0("splits/", d, "/")
  cli::cli_progress_bar(total = nreps * length(torun), type = "iterator")
  for (r in 1:nreps){ 
    id_test <- readRDS(paste0(split_path, r, "_id_test.rds"))
    train <- data[-id_test,]
    test <- data[id_test,]
    for (c_name in torun){
      c_fun <- get(c_name)
      time <- system.time(predict <- c_fun(train, test))[3]
      filename <- paste0(res_path, c_name, "_", r, ".rds")
      saveRDS(list(time = time, predict = predict), file = filename)
      cli::cli_progress_update(status = c_name)
    }
  }
}



