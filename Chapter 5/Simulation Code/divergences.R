library("Rcpp")
sourceCpp("divergences.cpp")

## symmetrized kaniadakis
kaniadakis <- function(data) {
  D <- Ckaniadakis(data)
  D <- 0.5 * (D + t(D))
  D <- D[lower.tri(D)]
  attr(D, which = "Size") <- nrow(data)
  attr(D, "Labels") <- rownames(data)
  class(D) <- "dist"
  return(D)
}



### aka symmetrized Kullback-Leibler
jeffreys <- function(data){
  D <- Cjeffreys(data)
  D <- D[lower.tri(D)]
  attr(D, which = "Size") <- nrow(data)
  attr(D, "Labels") <- rownames(data)
  class(D) <- "dist"

  return(D)
}

## jensenshannon
jensenshannon <- function(data){
  D <- Cjensensshannon(data)
  D <- D[lower.tri(D)]
  attr(D, which = "Size") <- nrow(data)
  attr(D, "Labels") <- rownames(data)
  class(D) <- "dist"

  return(D)
}

# symmetrized Renyi Divergence
renyi <- function(data, alpha = 0.1) {
  if (alpha <= 0 || alpha == 1) {
    stop("Parameter alpha must be greater than 0 and not equal to 1.")
  }
  D <- Creny(data, alpha)
  D <- 0.5 * (D + t(D))
  D <- D[lower.tri(D)]
  attr(D, which = "Size") <- nrow(data)
  attr(D, "Labels") <- rownames(data)
  class(D) <- "dist"

  return(D)
}

# Define the divergence measures to be used
divergence_measures <- list(
  jeffreys = jeffreys,
  kaniadakis = kaniadakis,
  renyi = renyi,
  hellinger = "hellinger", 
  totvar = "totvar",
  euclidean = "euclidean",
  jensenshannon = jensenshannon
)

