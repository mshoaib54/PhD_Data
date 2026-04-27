source("generate_data.R")

# define paths
data_dir = "simulated_data"

# Define parameter values for data generation
p_values <- c(5, 7, 9, 11, 13, 15)
q_values <- c(0.5, 0.9)
ktrue_values <- c(2)
sample_sizes <- 2**c(7, 9, 11, 13)
reps <- 100


### generate and save data ####
dir.create(data_dir, showWarnings = FALSE)


all_parameters <- rbind(
  expand.grid(p = p_values, ktrue = ktrue_values, q = NA),
  expand.grid(p = p_values, ktrue = NA, q = q_values)
)

pb <- progress_bar$new(format = "make data [:bar] :percent :eta",
                       total = reps * nrow(all_parameters))
for (i in seq_len(nrow(all_parameters))){
  p <- all_parameters$p[i]
  ktrue <- all_parameters$ktrue[i]
  q <- all_parameters$q[i]
  dir_out <- file.path(data_dir, paste0("p", p), paste0("q", q),
                       paste0("k", ktrue))
  dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
  for (rep in seq_len(reps)){
    gt <- generate_data(p = p, q = q, ktrue = ktrue,
                        sample_size = max(sample_sizes),
                        rfun = rexp)
    saveRDS(object = gt, file = file.path(dir_out, paste0("rep", rep, ".rds")))
    pb$tick()
  }
}
