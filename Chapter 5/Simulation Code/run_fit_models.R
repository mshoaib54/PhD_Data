
library("stagedtrees")
source("fit_models.R")
library("progress")



sample_sizes <- c(2^7,2^9,2^11,2^13)

# parameter values for fit models
k_values <- c(2,NA)
methods <- c("complete", "average",
             "ward.D2", "mcquitty")
divergences <- c('totvar', 'hellinger',
                 'jensenshannon', 'kaniadakis', 'jeffreys')

#methods <- c("complete")
#divergences <- c('totvar', 'kaniadakis')


#folders <- c("simulated_data/p5/","simulated_data/p7/", "simulated_data/p9/",
#            "simulated_data/p11/","simulated_data/p13/","simulated_data/p15/")

#folders <- c("simulated_data/p5/","simulated_data/p7/", "simulated_data/p9/")
folders <-c("simulated_data/p11/", "simulated_data/p13/", "simulater_data/p15/")
# Collect .rds files from all folders
all_files <- unlist(lapply(folders, function(path) {
  list.files(path,
             recursive = TRUE,
             pattern = "*.rds$",
             full.names = TRUE)
}))

dir.create("fitted_models", showWarnings = FALSE)
pb <- progress_bar$new(format = "fitting models [:bar] :percent :eta",
                       total = length(all_files))
for (file in all_files){
  results <- fit_models(file, sample_sizes = sample_sizes,
                        k_values = k_values,
                        methods = methods,
                        divergences = divergences,
                        fit_bhc = FALSE
  )
  dir.create(dirname(file.path("fitted_models", file)),
             showWarnings = FALSE, recursive = TRUE)
  saveRDS(results, file = file.path("fitted_models", file))
  pb$tick()
}
