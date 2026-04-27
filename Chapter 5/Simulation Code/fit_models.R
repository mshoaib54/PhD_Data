source("divergences.R")
library("stagedtrees")

fit_models <- function(data_path,
                       sample_sizes = c(100, 1000),
                       k_values = c(NA, 2),
                       methods = c("complete", "average",
                                   "ward.D2", "mcquitty"),
                       divergences = c("euclidean"),
                       lambda = 1e-6,
                       fit_bhc = TRUE) {

  dataAll <- readRDS(data_path)$data
  all_results <- list()
  for (size in sample_sizes){
  stopifnot(size <= nrow(dataAll))
  data <- head(dataAll, size)
  time_full <- system.time(fm1 <- full(data, lambda = lambda))["user.self"]

  # Filter distance methods if a subset is provided
  if (is.null(divergences)) {
    divergences <- names(divergence_measures)
  }

  results <- list()
   # Loop through methods
  for (method in methods){
    # Loop through each k value
    for (k_val in k_values) {
      # Loop through each selected divergence measure, fit model, and measure time
      for (divergence in divergences) {
        distance_function <- divergence_measures[[divergence]]
        model_time <- system.time(
        model <- tryCatch(
            stages_hclust(fm1, distance = distance_function,
                          k = k_val, max_k = Inf, first_max = TRUE, method = method),
            error = function(e) {
              message(paste("Error in", divergence, " k =", k_val,
                            " method =",method, ":", e$message))
              return(NA)
            }
          )
        , gcFirst = TRUE)["user.self"]

        # Store model, time, and k value
        res <- list(
          model = model,
          time = model_time,  # Store only user time
          k = k_val,  # Store k value
          alg = "hclust",
          divergence = divergence,
          method = method,
          data_path = data_path,
          sample_size = size
        )

        results[[paste("hclust", method, divergence, k_val, sep = "_")]] <- res
      }
    }
  }

  # Special case for BHC
  if (ncol(data) <= 7 & fit_bhc) {
    bhc_time <- system.time(bhc_model <- stages_bhc(fm1))["user.self"]
  } else {
    bhc_time <- NA
    bhc_model <- NA
  }
    results[["bhc"]] <- list(
      model = bhc_model,
      alg = "bhc",
      time = bhc_time,
      data_path = data_path,
      sample_size = size,
      k = NA,
      divergence = NA,
      method = NA
    )


  results[["full"]] <- list(
    model = fm1,
    alg = "full",
    time = time_full,
    data_path = data_path,
    sample_size = size,
    k = NA,
    divergence = NA,
    method = NA
  )

  all_results[[paste0("size", size)]] <- results
  }
  return(all_results)
}
