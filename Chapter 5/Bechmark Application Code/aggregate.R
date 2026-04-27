datasets <- read.table("binary_datasets_names_final.tsv")[, 1]

classifiers <- c("st_naive_totvar_cmi", 
                 "st_naive_jefferys_cmi",
                 "st_naive_kaniadakis_cmi",
                 "st_naive_jensen_shannon_cmi",
                 "st_naive_hellinger_cmi")




source("statistics.R")
nreps <- 10

TABLE <- array(
  data = NA,
  dim = c(
    length(statistics),
    length(datasets),
    length(classifiers),
    nreps
  ),
  dimnames = list(
    stat = statistics,
    data = datasets,
    classifier = classifiers,
    rep = 1:nreps
  )
)

for (d in datasets) {
  res_path <- paste0("results/", d, "/")
  data <- readRDS(paste0("datasets/", d, ".rds"))
  split_path <- paste0("splits/", d, "/")
  for (r in 1:nreps) {  
    id_test <- readRDS(paste0(split_path, r, "_id_test.rds"))
    true <- data$answer[id_test]
    for (c_name in classifiers) {
      filename <- paste0(res_path, c_name, "_", r, ".rds" )
      if (file.exists(filename)){
        res <- readRDS(filename)
        for (stat in statistics){
          stat_fun <- get(stat)
          TABLE[stat, d, c_name, r] <- stat_fun(res, true)
        }
      }
    }
  }
}

AVG <- apply(TABLE, c(1,2,3), mean, na.rm = TRUE)

saveRDS(TABLE, "TABLE.rds")

saveRDS(AVG, "AVG.rds")


