rm(list = ls())
all_names <- read.table("datasets_names.tsv")[, 1]

nreps <- 10
k <- 5

splits_path <- "splits/"
dir.create(splits_path, showWarnings = FALSE)
for (nam in all_names){
  path <- paste0(splits_path, nam, "/")
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  data <- readRDS(paste0("datasets/", nam, ".rds" ))
  N <- NROW(data)
  ## replicability dataset-wise
  set.seed(2020)
  for (r in 1:nreps){
    zero_count <- 1
    while(zero_count > 0) {
      zero_count <- 0
      id_test <- sample(1:N, N / k, replace = FALSE)
      marginal_counts <- rep(list(list()), NCOL(data))
      for(t in 1:NCOL(data)) {
        marginal_counts[[t]] <- table(data[-id_test, t])
        if(any(marginal_counts[[t]] == 0)) { 
          zero_count <- zero_count + 1 
        }
      }
    }
    saveRDS(id_test, file = paste0(path, r, "_id_test.rds"))
  }
}

