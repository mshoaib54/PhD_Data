# Load necessary libraries
library("stagedtrees")
library("progress")

# Function to create random staged event tree with k=2
# stages at each level
random_sevt_k <- function(tree, k = 2, rfun = rexp){
  model <- sevt(tree)
  model$stages <- lapply(model$stages, FUN = function(stages){
    paste0(sample(1:k, size = length(stages), replace = TRUE))
  })
  model <- random_sevt(model, q = 0, rfun = rfun)
  return(model)
}


# Function to generate and sample data from a staged tree model
generate_data <- function(p, q = NA, ktrue = NA,
                                 sample_size = 5000,
                                 rfun = rexp) {
  pl <- rep(list(c("0", "1")), p)
  names(pl) <- paste0("X", 1:p)

  true_model <- if (is.na(ktrue)) {
    random_sevt(pl, q = q, rfun = rfun)
  } else {
    random_sevt_k(pl, k = ktrue, rfun = rfun)
  }

  data <- sample_from(true_model, size = sample_size)

  return(list(
    true_model = true_model,
    data = data,
    p = p,
    q = q,
    ktrue = ktrue,
    sample_size = sample_size
  ))
}

