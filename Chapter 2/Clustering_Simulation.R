
# Exporting the libraries
#####

library(aricode)  # for Adjusted Rand Index
library(entropy)  # for Normalized Mutual Information and V-measure
library(Compositional)
library(robCompositions)
library(transport)
library(DirichletReg)
library(compositions)
library(RColorBrewer)
library(coda.base)
library(ggplot2)
library(reshape2)
library(plotly)
library(ggtern)
library(triptych)
library(zCompositions)
library(chemometrics)
library("Ternary")
library(nexus)
library(ContaminatedMixt)
library(fpc)
library(cluster)
library(grid)
library(flashClust)
library(mclust)
library(fpc)      # for Davies-Bouldin Index and Dunn Index
library(cluster)  # for silhouette calculations
library(aricode)  # for Adjusted Rand Index
library(entropy)  # for Normalized Mutual Information and V-measure
############
# Several Divergence and metric functions
#############
# Distances & divergences
# -------------------------

# Additive log-ratio distance
d_alr <- function(x, y, denom_index = NULL) {
  x <- as.numeric(x); y <- as.numeric(y)
  n <- length(x)
  if (length(y) != n) stop("x and y must have same dimension.")
  if (any(x <= 0) || any(y <= 0)) stop("Additive log-ratio distance requires strictly positive components.")
  if (is.null(denom_index)) denom_index <- n
  idx <- setdiff(1:n, denom_index)
  v <- log(x[idx] / x[denom_index]) - log(y[idx] / y[denom_index])
  sqrt(sum(v^2))
}

# Aitchison distance (scaled clr form)
# d_A(x,y) = sqrt( (1/n) sum_i ( clr_i(x) - clr_i(y) )^2 )
d_A <- function(x, y) {
  z <- clr(x) - clr(y)
  sqrt(mean(z^2))
}

# Pairwise (unscaled) variant (rarely used directly)
d_A_pairwise_unscaled <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  n <- length(x)
  if (any(x <= 0) || any(y <= 0)) stop("Requires strictly positive.")
  Lx <- outer(log(x), log(x), `-`)
  Ly <- outer(log(y), log(y), `-`)
  sqrt(sum((Lx - Ly)^2))
}

# Angular distance
d_angular <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (any(x < 0) || any(y < 0)) stop("Angular distance requires nonnegative components.")
  nx <- sqrt(sum(x^2)); ny <- sqrt(sum(y^2))
  if (nx == 0 || ny == 0) stop("Zero vector not allowed.")
  u <- sum((x/nx) * (y/ny))
  safe_acos(u)
}

# Bhattacharyya coefficient (after closure to 1)
BC <- function(x, y) {
  x <- close_vec(x, 1); y <- close_vec(y, 1)
  sum(sqrt(x * y))
}

# Bhattacharyya distance (standard; not always metric)
d_bhattacharyya <- function(x, y) {
  u <- BC(x, y)
  -log(pmax(u, 1e-15))
}

# Bhattacharyya angular form (metric on simplex)
d_bhattacharyya_angular <- function(x, y) {
  u <- BC(x, y)
  safe_acos(u)
}

# Fidelity distance
d_fidelity <- function(x, y) {
  u <- BC(x, y)
  sqrt(pmax(0, 1 - u))
}

# Fisher distance (geodesic on sqrt map)
d_fisher <- function(x, y) {
  u <- BC(x, y)
  2 * safe_acos(u)
}

# Hellinger distance
d_hellinger <- function(x, y) {
  x <- close_vec(x, 1); y <- close_vec(y, 1)
  sqrt(0.5 * sum((sqrt(x) - sqrt(y))^2))
}

# Compositional KL (symmetric/compositional form)
d_compositional_kl <- function(x, y) {
  x <- close_vec(x, 1); y <- close_vec(y, 1)
  if (any(x <= 0) || any(y <= 0)) stop("Compositional KL requires strictly positive components.")
  n <- length(x)
  a <- sum(x / y)
  b <- sum(y / x)
  sqrt((n/2) * log((a*b) / (n^2)))
}

# City-block (L1)
d_cityblock <- function(x, y) sum(abs(x - y))

# Total variation
d_total_variation <- function(x, y) 0.5 * sum(abs(x - y))

# Triangular discrimination
d_triangular_discrimination <- function(x, y) {
  num <- (x - y)^2
  den <- x + y
  if (any(den == 0)) stop("Triangular discrimination undefined when x_i + y_i = 0.")
  sum(num / den)
}

# Rényi divergence
d_renyi <- function(x, y, alpha = 0.5) {
  if (alpha <= 0 || alpha == 1) stop("alpha must be > 0 and alpha != 1.")
  x <- close_vec(x, 1); y <- close_vec(y, 1)
  if (any(x <= 0) || any(y <= 0)) stop("Rényi divergence requires strictly positive components.")
  (1/(alpha - 1)) * log(sum((x^alpha) * (y^(1 - alpha))))
}

# Zero-adjusted Aitchison distance (pseudocount + closure)
d_zero_adjusted_aitchison <- function(x, y, eps = 1e-6) {
  x2 <- add_eps_close(x, eps = eps, kappa = 1)
  y2 <- add_eps_close(y, eps = eps, kappa = 1)
  d_A(x2, y2)
}

# Mahalanobis (crude) with pseudoinverse covariance
d_mahalanobis_crude <- function(x, mu, Sigma) {
  x <- as.numeric(x); mu <- as.numeric(mu)
  if (length(x) != length(mu)) stop("x and mu must match.")
  Sinv <- MASS::ginv(Sigma)
  as.numeric(sqrt(t(x - mu) %*% Sinv %*% (x - mu)))
}

# Mahalanobis in clr space
d_mahalanobis_clr <- function(x, y, Sigma_clr) {
  zx <- clr(x); zy <- clr(y)
  Sinv <- MASS::ginv(Sigma_clr)
  as.numeric(sqrt(t(zx - zy) %*% Sinv %*% (zx - zy)))
}


# -------------------------
# ilr distance (Euclidean distance in ilr space)
# -------------------------

# (n-1) x n orthonormal contrast matrix (rows orthonormal, row-sums = 0)
# We use a normalized Helmert sub-matrix; any orthonormal basis gives the same ilr distance.
psi_helmert <- function(n) {
  H <- stats::contr.helmert(n)     # n x (n-1)
  Psi <- t(H)                      # (n-1) x n
  Psi <- Psi / sqrt(rowSums(Psi^2))
  Psi
}

d_ilr <- function(x, y, Psi = NULL) {
  x <- as.numeric(x); y <- as.numeric(y)
  n <- length(x)
  if (length(y) != n) stop("x and y must have the same dimension.")
  if (any(x <= 0) || any(y <= 0)) stop("ilr distance requires strictly positive components.")
  if (is.null(Psi)) Psi <- psi_helmert(n)
  zx <- ilr(x, Psi)
  zy <- ilr(y, Psi)
  sqrt(sum((zx - zy)^2))
}

# -------------------------
# Kaniadakis divergence (escort-based; as in your definitions)
# -------------------------

A <- function(p) 2 * (p^2) / (1 + p^2)                 # Nudit’s function

ptilde <- function(p) {                               # escort distribution
  ap <- A(p)
  ap / sum(ap)
}

Ek <- function(p, u) sum(u * ptilde(p))               # escort expectation eE_p[u]

Kl <- function(p, eps = 1e-12) {                      # Kaniadakis log (your definition)
  p <- pmax(as.numeric(p), eps)
  0.5 * (p - 1/p)
}

Kexp <- function(x) {                                 # inverse of Kl: kexp
  x <- as.numeric(x)
  x + sqrt(x^2 + 1)
}

Kdiv <- function(p, q) {                              # Kaniadakis divergence
  Ek(p, Kl(p)) - Ek(p, Kl(q))
}

d_kaniadakis <- function(x, y) {
  # Typically used on the probability simplex (closed to 1)
  x <- close_vec(x, 1)
  y <- close_vec(y, 1)
  Kdiv(x, y)
}

#####
# Function to calculate Calinski index


# Define metric calculation functions
compute_wcss <- function(data, labels) {
  wcss <- function(data, centers, cluster) {
    sum(sapply(unique(cluster), function(k) {
      sum(rowSums((data[cluster == k, ] - centers[k, ])^2))
    }))
  }
  
  # Compute cluster centers
  centers <- aggregate(data, by = list(cluster = labels), FUN = mean)
  centers <- centers[, -1]  # remove the cluster column
  
  # Calculate WCSS
  wcss_value <- wcss(data, centers, labels)
  return(wcss_value)
}

compute_ch_index <- function(data, labels) {
  ch_index <- calinhara(data, labels)
  return(ch_index)
}

compute_avg_silhouette_width <- function(data, labels) {
  sil_values <- silhouette(labels, dist(data))[, 3]
  # Exclude non-positive silhouette values before computing the geometric mean
  sil_values <- sil_values[sil_values > 0]
  if (length(sil_values) == 0) {
    return(NA)
  }
  geometric_mean_sil_width <- exp(mean(log(sil_values)))
  return(geometric_mean_sil_width)
}

compute_davies_bouldin <- function(data, labels) {
  clusters <- unique(labels)
  n_clusters <- length(clusters)
  centroids <- sapply(clusters, function(k) {
    colMeans(data[labels == k, , drop = FALSE])
  }, simplify = "matrix")
  
  s <- sapply(clusters, function(k) {
    cluster_data <- data[labels == k, , drop = FALSE]
    mean(sqrt(rowSums((cluster_data - centroids[, k])^2)))
  })
  
  d <- matrix(NA, n_clusters, n_clusters)
  for (i in 1:n_clusters) {
    for (j in 1:n_clusters) {
      if (i != j) {
        d[i, j] <- sqrt(sum((centroids[, i] - centroids[, j])^2))
      }
    }
  }
  
  r <- matrix(NA, n_clusters, n_clusters)
  for (i in 1:n_clusters) {
    for (j in 1:n_clusters) {
      if (i != j) {
        r[i, j] <- (s[i] + s[j]) / d[i, j]
      }
    }
  }
  
  dbi <- mean(apply(r, 1, max, na.rm = TRUE))
  return(dbi)
}

compute_dunn_index <- function(data, labels) {
  clusters <- unique(labels)
  n_clusters <- length(clusters)
  
  # Compute the minimum inter-cluster distance
  min_inter_cluster_dist <- Inf
  for (i in 1:(n_clusters - 1)) {
    for (j in (i + 1):n_clusters) {
      cluster_i <- data[labels == clusters[i], , drop = FALSE]
      cluster_j <- data[labels == clusters[j], , drop = FALSE]
      dist_ij <- as.matrix(dist(rbind(cluster_i, cluster_j)))
      inter_cluster_dist <- min(dist_ij[1:nrow(cluster_i), (nrow(cluster_i) + 1):(nrow(cluster_i) + nrow(cluster_j))])
      if (inter_cluster_dist < min_inter_cluster_dist) {
        min_inter_cluster_dist <- inter_cluster_dist
      }
    }
  }
  
  # Compute the maximum intra-cluster
  # Compute the maximum intra-cluster distance
  max_intra_cluster_dist <- 0
  for (k in clusters) {
    cluster_k <- data[labels == k, , drop = FALSE]
    if (nrow(cluster_k) > 1) {  # Ensure there is more than one point in the cluster
      intra_cluster_dists <- dist(cluster_k)
      max_dist_in_cluster <- max(intra_cluster_dists)
      if (max_dist_in_cluster > max_intra_cluster_dist) {
        max_intra_cluster_dist <- max_dist_in_cluster
      }
    }
  }
  
  # Compute Dunn Index
  dunn_index <- min_inter_cluster_dist / max_intra_cluster_dist
  return(dunn_index)
}


#####
# Dummy Example
#####
# Example compositions x and y
x <- matrix(c(0.3, 0.5, 0.2), nrow = 1, ncol = 3)
y <- matrix(c(0.5, 0.3, 0.2), nrow = 1, ncol = 3)
z <- matrix(c(0.6,0.2,0.2), nrow = 1, ncol = 3)
perturb<- function(x,z) x*z/sum(x*z)
x_1 <- perturb(x,z)
y_1 <- perturb(y,z)
# Checking the on Sub_Compositions

#####
# Checking on the subcomposition
#####
sx_12 <- matrix(x[c(1, 2)], nrow = 1, ncol = 2)
sy_23 <- matrix(y[c(1, 2)], nrow = 1, ncol = 2)
s_12<-t(apply(sx_12,1, function(x) x/sum(x)))
s_23 <-t(apply(sy_23,1, function(x) x/sum(x)))

#####
# Examples
######
# Aitchison Distance

aDist(x, y) 
aDist(x_1, y_1)
aDist(s_12,s_23)
# CKL Divergence
d_compositional_kl(x, y)
d_compositional_kl(x_1, y_1)
d_compositional_kl(s_12, s_23)
# KD Divergence
d_kaniadakis(x,y)
d_kaniadakis(x_1,y_1)
d_kaniadakis(s_12,s_23)
# Angular Distance
d_ang(x,y)
d_ang(x_1,y_1)
d_ang(s_12,s_23)
# Hellinger Distance
d_hellinger(x,y)
d_hellinger(x_1,y_1)
d_hellinger(s_12,s_23)
# J Divergence
d_fisher(x,y)
d_fisher(x_1,y_1)
d_fisher(s_12,s_23)


kmeans_compositional <- function(data, numClusters,
                                 distFunc = NULL,   # string key
                                 distFun  = NULL,   # function(x,y)
                                 maxIterations = 100, tolerance = 1e-6, seed = 1) {
  
  dataMatrix <- as.matrix(data)
  
  availableDistFuncs <- list(
    "alr"                       = d_alr,
    "aitchison"                 = d_A,
    "aitchison_zero_adjusted"   = d_zero_adjusted_aitchison,
    "aitchison_pairwise_raw"    = d_A_pairwise_unscaled,
    "angular"                   = d_angular,
    "bhattacharyya"             = d_bhattacharyya,
    "bhattacharyya_angular"     = d_bhattacharyya_angular,
    "fidelity"                  = d_fidelity,
    "fisher"                    = d_fisher,
    "hellinger"                 = d_hellinger,
    "compositional_kl"          = d_compositional_kl,
    "cityblock"                 = d_cityblock,
    "total_variation"           = d_total_variation,
    "triangular_discrimination" = d_triangular_discrimination,
    "renyi_alpha_0_5"           = function(x, y) d_renyi(x, y, alpha = 0.5),
    "ilr"                       = d_ilr,
    "kaniadakis"                = d_kaniadakis
  )
  
  # ---- choose distance function ----
  if (!is.null(distFun)) {
    if (!is.function(distFun)) stop("distFun must be a function(x,y).")
    selectedDistFunc <- distFun
  } else {
    if (is.null(distFunc)) stop("Provide either distFunc (string) or distFun (function).")
    if (!distFunc %in% names(availableDistFuncs)) {
      stop("Unknown distFunc. Available:\n  ", paste(names(availableDistFuncs), collapse = "\n  "))
    }
    selectedDistFunc <- availableDistFuncs[[distFunc]]
  }
  
  # ---- init ----
  set.seed(seed)
  n <- nrow(dataMatrix)
  centroids <- dataMatrix[sample(n, numClusters), , drop = FALSE]
  
  for (iter in 1:maxIterations) {
    
    distsToCentroids <- sapply(1:numClusters, function(k) {
      apply(dataMatrix, 1, function(x) selectedDistFunc(x, centroids[k, ]))
    })
    
    clusterAssignments <- max.col(-distsToCentroids)
    
    newCentroids <- t(sapply(1:numClusters, function(k) {
      Xk <- dataMatrix[clusterAssignments == k, , drop = FALSE]
      if (nrow(Xk) > 0) colMeans(Xk, na.rm = TRUE) else centroids[k, ]
    }))
    
    if (sqrt(sum((centroids - newCentroids)^2)) < tolerance) break
    centroids <- newCentroids
  }
  
  # ---- metrics (your functions must exist) ----
  wcss <- compute_wcss(dataMatrix, clusterAssignments)
  ch   <- compute_ch_index(dataMatrix, clusterAssignments)
  sil  <- compute_avg_silhouette_width(dataMatrix, clusterAssignments)
  dbi  <- compute_davies_bouldin(dataMatrix, clusterAssignments)
  dunn <- compute_dunn_index(dataMatrix, clusterAssignments)
  
  list(
    clusterAssignments = clusterAssignments,
    centroids = centroids,
    withinClusterSumOfSquares = wcss,
    calinskiHarabaszIndex = ch,
    averageSilhouetteWidth = sil,
    daviesBouldinIndex = dbi,
    dunnIndex = dunn,
    iterations = iter
  )
}


#####
# Simulation Study
#####

## =========================
## 4) Simulation Study + clustering (ALL measures)
## =========================
set.seed(123)
N <- 200  # Number of samples
D <- 3    # Number of components
df <- rdirichlet(N, rep(1, D))  # Generate compositional data
#log_ratios <- log(data[, c("X1", "X2", "X3")]) - colMeans(log(data[, c("X1", "X2", "X3")]))
# Perform K-means clustering for compositional data
K <- 3


## --- list ALL pairwise measures that work directly on compositions ---
## Note: some are distances, some are divergences; kmeans just needs nonnegative dissimilarity.
## Replace your old availableDistFuncs with this FULL catalog
## (uses the exact function names you defined above)



## Run clustering for each measure
cluster_results <- purrr::imap(dist_list, function(f, name) {
  kmeans_compositional(df, numClusters = K, distFun = f, seed = 123)
})

## Attach cluster labels to data with systematic names
data_df <- as.data.frame(df)   # df is your 200x3 data.frame from Dirichlet

for (nm in names(cluster_results)) {
  colname <- paste0("cluster_", gsub("[^A-Za-z0-9]+", "_", nm))
  cl <- cluster_results[[nm]]$clusterAssignments
  
  stopifnot(nrow(data_df) == length(cl))  # sanity check
  
  data_df[[colname]] <- factor(cl)
}

## =========================
## 5) Facet-wrapped ternary plot for ALL clusterings
## =========================
methods_tbl <- tibble::tibble(
  Method = names(cluster_results),
  cluster_col = paste0("cluster_", gsub("[^A-Za-z0-9]+", "_", names(cluster_results)))
)
data_long <- purrr::pmap_dfr(methods_tbl, function(Method, cluster_col) {
  data_df %>%
    transmute(
      V1 = V1, V2 = V2, V3 = V3,
      Cluster = factor(.data[[cluster_col]]),
      Method = Method
    )
})


library(dplyr)
library(purrr)
library(stringr)
library(ggtern)
library(ggplot2)
library(RColorBrewer)

## -----------------------------
## 1) Ensure UNIQUE method names
## -----------------------------
orig_names <- names(all_results)
names(all_results) <- make.unique(orig_names)   # avoids duplicates like Mahalanobis_clr, Renyi...

## Pretty labels for facets (same mapping everywhere)
pretty_method <- function(x) {
  x %>%
    gsub("_", " ", .) %>%
    gsub("\\s+", " ", .) %>%
    trimws()
}

method_map <- tibble(
  MethodKey   = names(all_results),
  MethodPretty = pretty_method(names(all_results))
)

## -----------------------------
## 2) Attach cluster columns to data_df
## -----------------------------
data_df <- as.data.frame(df)
colnames(data_df) <- c("V1","V2","V3")

for (nm in names(all_results)) {
  colname <- paste0("cluster__", nm)   # keep key, no gsub collisions
  data_df[[colname]] <- factor(all_results[[nm]]$clusterAssignments)
}

## -----------------------------
## 3) Build LONG data for points (consistent MethodPretty)
## -----------------------------
data_long <- method_map %>%
  mutate(cluster_col = paste0("cluster__", MethodKey)) %>%
  pmap_dfr(function(MethodKey, MethodPretty, cluster_col) {
    data_df %>%
      transmute(
        V1 = V1, V2 = V2, V3 = V3,
        Cluster = .data[[cluster_col]],
        Method  = MethodPretty
      )
  })

## -----------------------------
## 4) Build LONG data for centroids (same MethodPretty)
## -----------------------------
centroids_long <- method_map %>%
  pmap_dfr(function(MethodKey, MethodPretty) {
    C <- as.data.frame(all_results[[MethodKey]]$centroids)
    colnames(C) <- c("V1","V2","V3")
    C$Cluster <- factor(seq_len(nrow(C)))   # centroid id
    C$Method  <- MethodPretty
    C
  })

## -----------------------------
## 5) Dark palette (clusters)
## -----------------------------
K <- length(levels(data_long$Cluster))
pal_dark <- if (K <= 8) brewer.pal(K, "Dark2") else grDevices::hcl.colors(K, "Dark 3")
names(pal_dark) <- levels(data_long$Cluster)

## -----------------------------
## 6) Plot (cleaner theme + wrapped strip titles)
## -----------------------------
p_facet_clean <- ggtern(data_long, aes(V1, V2, V3, color = Cluster)) +
  geom_point(alpha = 0.80, size = 1.4) +
  geom_point(
    data = centroids_long,
    mapping = aes(V1, V2, V3),
    inherit.aes = FALSE,
    shape = 4, size = 3.2, stroke = 1.2,
    color = "black"
  ) +
  facet_wrap(~ Method, nrow = 3, labeller = label_wrap_gen(width = 16)) +
  scale_color_manual(values = pal_dark, name = "Cluster") +
  theme_rgbw() +
  theme(
    # make facets readable
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    
    # reduce ternary clutter
    tern.axis.text = element_text(size = 6, face = "bold"),
    tern.axis.title = element_blank()
  ) +
  labs(title = "")

print(p_facet_clean)

ggsave("Ternary_clusters_ALL_methods_with_centroids_clean.pdf",
       p_facet_clean, device = cairo_pdf, width = 18, height = 10)




## =========================
## Run compositional k-means with ALL your measures
## (uses kmeans_compositional() that accepts distFun OR distFunc)
## =========================

## ---- make sure df is a data.frame with columns V1,V2,V3 ----
# df <- as.data.frame(MCMCpack::rdirichlet(N, rep(1, D)))
# colnames(df) <- c("V1","V2","V3")

K <- 3

## ============================================================
## 1) Measures that are directly usable as pairwise dissimilarities
##    (function(x,y) only)
## ============================================================
result_alr        <- kmeans_compositional(df, numClusters = K, distFun = d_alr, seed = 123)
result_A          <- kmeans_compositional(df, numClusters = K, distFun = d_A, seed = 123)
result_A0         <- kmeans_compositional(df, numClusters = K, distFun = d_zero_adjusted_aitchison, seed = 123)
#result_A_raw      <- kmeans_compositional(df, numClusters = K, distFun = d_A_pairwise_unscaled, seed = 123)

result_ang        <- kmeans_compositional(df, numClusters = K, distFun = d_angular, seed = 123)

result_bhat       <- kmeans_compositional(df, numClusters = K, distFun = d_bhattacharyya, seed = 123)
result_bhat_ang   <- kmeans_compositional(df, numClusters = K, distFun = d_bhattacharyya_angular, seed = 123)

result_fid        <- kmeans_compositional(df, numClusters = K, distFun = d_fidelity, seed = 123)
result_fisher     <- kmeans_compositional(df, numClusters = K, distFun = d_fisher, seed = 123)

result_hel        <- kmeans_compositional(df, numClusters = K, distFun = d_hellinger, seed = 123)
result_ckl        <- kmeans_compositional(df, numClusters = K, distFun = d_compositional_kl, seed = 123)

result_l1         <- kmeans_compositional(df, numClusters = K, distFun = d_cityblock, seed = 123)
result_tv         <- kmeans_compositional(df, numClusters = K, distFun = d_total_variation, seed = 123)
result_tri        <- kmeans_compositional(df, numClusters = K, distFun = d_triangular_discrimination, seed = 123)

result_renyi      <- kmeans_compositional(df, numClusters = K,
                                          distFun = function(x,y) d_renyi(x,y, alpha = 0.5),
                                          seed = 123)

result_ilr        <- kmeans_compositional(df, numClusters = K, distFun = d_ilr, seed = 123)

result_kaniadakis <- kmeans_compositional(df, numClusters = K, distFun = d_kaniadakis, seed = 123)

## ============================================================
## 2) Mahalanobis measures (need sample-dependent Sigma / mu)
##    We wrap them as functions(x,y) so they can be used in kmeans.
## ============================================================
X <- as.matrix(df)

## crude Mahalanobis depends on a reference mean mu and covariance Sigma
mu_crude    <- colMeans(X)
Sigma_crude <- cov(X)
dist_maha_crude_wrapped <- function(x, y) d_mahalanobis_crude(x, mu = mu_crude, Sigma = Sigma_crude)

## clr Mahalanobis depends on Sigma_clr
# clr() must exist in your environment; if not, replace clr(x) with clr_vec(x)
X_clr    <- t(apply(X, 1, clr))          # uses your clr() function
Sigma_clr <- cov(X_clr)
dist_maha_clr_wrapped <- function(x, y) d_mahalanobis_clr(x, y, Sigma_clr = Sigma_clr)

result_maha_crude <- kmeans_compositional(df, numClusters = K, distFun = dist_maha_crude_wrapped, seed = 123)
result_maha_clr   <- kmeans_compositional(df, numClusters = K, distFun = dist_maha_clr_wrapped, seed = 123)

## ============================================================
## 3) Collect everything in one object (optional, convenient)
## ============================================================
all_results <- list(
  ALR = result_alr,
  Aitchison = result_A,
  Aitchison_zero_adjusted = result_A0,
  #Aitchison_pairwise_raw = result_A_raw,
  Angular = result_ang,
  Bhattacharyya = result_bhat,
  Bhattacharyya_angular = result_bhat_ang,
  Fidelity = result_fid,
  Fisher = result_fisher,
  Hellinger = result_hel,
  CKL = result_ckl,
  Cityblock_L1 = result_l1,
  Total_variation = result_tv,
  Triangular_discrimination = result_tri,
  Renyi_alpha_0_5 = result_renyi,
  ILR = result_ilr,
  Kaniadakis = result_kaniadakis,
  Mahalanobis_crude = result_maha_crude,
  Mahalanobis_clr = result_maha_clr
)

## =========================
## Extract ALL cluster labels (one variable per method)
## =========================

cluster_ALR                     <- all_results$ALR$clusterAssignments
cluster_Aitchison               <- all_results$Aitchison$clusterAssignments
cluster_Aitchison_zero_adjusted <- all_results$Aitchison_zero_adjusted$clusterAssignments
#cluster_Aitchison_pairwise_raw  <- all_results$Aitchison_pairwise_raw$clusterAssignments

cluster_Angular                 <- all_results$Angular$clusterAssignments

cluster_Bhattacharyya           <- all_results$Bhattacharyya$clusterAssignments
cluster_Bhattacharyya_angular   <- all_results$Bhattacharyya_angular$clusterAssignments
cluster_Fidelity                <- all_results$Fidelity$clusterAssignments
cluster_Fisher                  <- all_results$Fisher$clusterAssignments

cluster_Hellinger               <- all_results$Hellinger$clusterAssignments
cluster_CKL                     <- all_results$CKL$clusterAssignments

cluster_Cityblock_L1            <- all_results$Cityblock_L1$clusterAssignments
cluster_Total_variation         <- all_results$Total_variation$clusterAssignments
cluster_Triangular_discrimination <- all_results$Triangular_discrimination$clusterAssignments

cluster_Renyi_alpha_0_5         <- all_results$Renyi_alpha_0_5$clusterAssignments

cluster_ILR                     <- all_results$ILR$clusterAssignments

cluster_Kaniadakis              <- all_results$Kaniadakis$clusterAssignments

cluster_Mahalanobis_crude       <- all_results$Mahalanobis_crude$clusterAssignments
cluster_Mahalanobis_clr         <- all_results$Mahalanobis_clr$clusterAssignments
# Centriods 
centroids <- all_results$centroids
# Convert data and centroids to data frames for plotting
data_df <- as.data.frame(df)
centroids_list <- lapply(all_results, `[[`, "centroids")


## =========================
## Build centroids_df for ALL methods + build data_df with ALL clusters
## =========================



# ---- 1) Data for plotting (keep as data.frame) ----
data_df <- as.data.frame(df)   # df has V1,V2,V3

# ---- 2) Add ALL cluster assignments (from all_results) ----
# creates columns like cluster_Aitchison, cluster_Kaniadakis, ...
for (nm in names(all_results)) {
  colname <- paste0("cluster_", gsub("[^A-Za-z0-9]+", "_", nm))
  data_df[[colname]] <- factor(all_results[[nm]]$clusterAssignments)
}

# ---- 3) Collect ALL centroids into one long data frame ----
centroids_df <- purrr::imap_dfr(all_results, function(res, method) {
  C <- as.data.frame(res$centroids)
  colnames(C) <- c("V1","V2","V3")
  C$cluster <- factor(seq_len(nrow(C)))
  C$Method  <- method
  C
})

# ---- 4) Example: ternary plot of raw data (no clustering) ----
p0 <- ggtern(data = data_df, aes(V1, V2, V3)) +
  geom_point(alpha = 0.5, size = 2) +
  theme_rgbw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))

print(p0)
ggsave("Ternary_clusters_ALL_methods_facetwrap.pdf",
       p0, device = cairo_pdf, width = 10, height = 10)
# ---- 5) Facet ternary plot: points colored by cluster for EACH method ----
# reshape to long format (Method-wise clustering)
cluster_cols <- grep("^cluster_", names(data_df), value = TRUE)

data_long <- data_df %>%
  pivot_longer(cols = all_of(cluster_cols),
               names_to = "Method",
               values_to = "Cluster") %>%
  mutate(Method = gsub("^cluster_", "", Method),
         Method = gsub("_", " ", Method))

p_clusters <- ggtern(data_long, aes(V1, V2, V3, color = Cluster)) +
  geom_point(alpha = 0.85, size = 1.6) +
  facet_wrap(~ Method, nrow = 3) +
  theme_rgbw() +
  labs(color = "Cluster") +
  theme(strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom")

print(p_clusters)

## Centroids in a different color (black), while points keep cluster colors

# Fixed palette for cluster points
cluster_levels <- levels(data_long$Cluster)
pal <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A",
         "4" = "#984EA3", "5" = "#FF7F00")[cluster_levels]

p_clusters_centroids <- ggtern(data_long, aes(V1, V2, V3, color = Cluster)) +
  geom_point(alpha = 0.85, size = 1.6) +
  # centroids: constant color (black) + larger X markers
  geom_point(
    data = centroids_df,
    mapping = ggtern::aes(x = V1, y = V2, z = V3),
    inherit.aes = FALSE,
    shape = 4, size = 3.4, stroke = 1.3,
    color = "black"
  ) +
  facet_wrap(~ Method, nrow = 3) +
  scale_color_manual(values = pal, name = "Cluster") +
  theme_rgbw() +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom"
  )

print(p_clusters_centroids)
## ============================================================
## Plot ternary clusters for ALL methods (auto-facet, no manual repeats)
## Assumes:
##   - data_df has columns V1,V2,V3 and cluster_* columns for every method
##   - all_results exists (or at least the cluster_* columns exist)
## ============================================================



## ---- 1) Make sure you are using a DATA FRAME (not matrix) ----
data_df <- as.data.frame(df)  # df must have V1,V2,V3

## ---- 2) Add ALL cluster columns from all_results (if not already added) ----
for (nm in names(all_results)) {
  colname <- paste0("cluster_", gsub("[^A-Za-z0-9]+", "_", nm))
  data_df[[colname]] <- factor(all_results[[nm]]$clusterAssignments)
}

## ---- 3) Convert wide -> long (one panel per method) ----
cluster_cols <- grep("^cluster_", names(data_df), value = TRUE)

data_long <- data_df %>%
  pivot_longer(cols = all_of(cluster_cols),
               names_to = "Method",
               values_to = "Cluster") %>%
  mutate(
    Method = gsub("^cluster_", "", Method),
    Method = gsub("_", " ", Method),
    Cluster = factor(Cluster)
  )

## ---- 4) Fixed cluster palette (works even if K > 3) ----
K <- nlevels(data_long$Cluster)
pal <- RColorBrewer::brewer.pal(max(3, min(8, K)), "Set1")
pal <- setNames(rep(pal, length.out = K), levels(data_long$Cluster))

## ---- 5) ONE plot: all methods with facet_wrap ----
p_all <- ggtern(data_long, aes(V1, V2, V3, color = Cluster)) +
  geom_point(alpha = 0.85, size = 1.8) +
  facet_wrap(~ Method, nrow = 3) +
  scale_color_manual(values = pal, name = "Cluster") +
  theme_rgbw() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(title = "Ternary clustering results across all distances/divergences")

print(p_all)

## Optional: save
# ggsave("Ternary_All_Methods_Facet.pdf", p_all, width = 14, height = 10)



## ============================================================
## PCA on CLR data + cluster overlays for ALL methods (facet plot)
## Assumes:
##   - df is your N x 3 data.frame with columns V1,V2,V3
##   - all_results is the list of clustering outputs (each has $clusterAssignments)
## ============================================================

## -------------------------
## 1) CLR transform + PCA
## -------------------------
clr_transformation <- function(x, eps = 1e-15) {
  x <- pmax(as.numeric(x), eps)
  gm <- exp(mean(log(x)))
  log(x / gm)
}

data_clr <- t(apply(as.matrix(df), 1, clr_transformation))
colnames(data_clr) <- colnames(df)

pca_res <- prcomp(data_clr, center = FALSE, scale. = FALSE)

scores <- as.data.frame(pca_res$x[, 1:2])
colnames(scores) <- c("PC1", "PC2")

loadings <- as.data.frame(pca_res$rotation[, 1:2])
loadings$variable <- rownames(loadings)

## scale arrows nicely
max_scale <- max(abs(scores$PC1), abs(scores$PC2)) * 0.9
arrows_data <- loadings %>%
  transmute(
    x = 0, y = 0,
    xend = PC1 * max_scale,
    yend = PC2 * max_scale,
    variable = variable
  )

## -------------------------
## 2) Build one long table of PCA scores + clusters for ALL methods
## -------------------------
methods_tbl <- tibble(
  Method = names(all_results),
  cluster_col = names(all_results)  # will map by lookup below
)

scores_long <- purrr::map_dfr(names(all_results), function(m) {
  tibble(
    PC1 = scores$PC1,
    PC2 = scores$PC2,
    Cluster = factor(all_results[[m]]$clusterAssignments),
    Method = m
  )
})

## Optional: clean method labels
scores_long <- scores_long %>%
  mutate(Method = gsub("[^A-Za-z0-9]+", " ", Method))

## -------------------------
## 3) One facet PCA plot for ALL methods
## -------------------------
p_pca_all <- ggplot(scores_long, aes(PC1, PC2, color = Cluster)) +
  geom_point(alpha = 0.85, size = 1.8) +
  geom_segment(
    data = arrows_data,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(0.18, "inches")),
    alpha = 0.25,
    color = "gray40"
  ) +
  geom_text(
    data = arrows_data,
    aes(x = xend, y = yend, label = variable),
    inherit.aes = FALSE,
    size = 3.2,
    color = "black",
    fontface = "bold",
    hjust = 1.15,
    vjust = 1.15
  ) +
  facet_wrap(~ Method, nrow = 3) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  ) +
  labs(
    title = "PCA on CLR-transformed compositions with clustering overlays",
    x = "PC1 (CLR-PCA)",
    y = "PC2 (CLR-PCA)",
    color = "Cluster"
  )

print(p_pca_all)


## -------------------------
## 1) Load and prepare real data
## -------------------------
data(alimentation)

# Use first 9 columns as you did
new_data <- alimentation[1:9]

# Closure (kappa=1)
closed_data <- clo(new_data)

# Keep three parts for ternary visualization: RM, WM, E
# (Make sure these exist in your dataset; if not, change names accordingly.)
closed_3 <- closed_data[, c("RM","WM","E")]
closed_3 <- as.data.frame(closed_3)

# Optional row names / labels
rowNames <- alimentation$Country

# Sanity checks
stopifnot(all(rowSums(closed_3) > 0))
any_na  <- any(is.na(closed_3))
any_inf <- any(is.infinite(as.matrix(closed_3)))
print(list(any_na = any_na, any_inf = any_inf))

K <- 3

## -------------------------
## 2) Define the list of measures to run (ALL)
## -------------------------
dist_list <- list(
  "Aitchison"                 = d_A,
  "Aitchison (zero-adjust)"   = d_zero_adjusted_aitchison,
  "C-KL"                      = d_compositional_kl,
  "Kaniadakis"                = d_kaniadakis,
  "Hellinger"                 = d_hellinger,
  "Fisher"                    = d_fisher,
  "Bhattacharyya"             = d_bhattacharyya,
  "Bhattacharyya (angular)"   = d_bhattacharyya_angular,
  "Fidelity"                  = d_fidelity,
  "Angular"                   = d_angular,
  "City-block (L1)"           = d_cityblock,
  "Total variation"           = d_total_variation,
  "Triangular discrimination" = d_triangular_discrimination,
  "Renyi (alpha=0.5)"         = function(x,y) d_renyi(x,y,alpha=0.5),
  "ilr distance"              = d_ilr,
  "ALR distance"              = function(x,y) d_alr(x,y,denom_index = 3)  # denom = E (3rd)
)

## Mahalanobis measures (sample-dependent) on the SAME dataset
X <- as.matrix(closed_3)
mu_crude    <- colMeans(X)
Sigma_crude <- cov(X)
dist_maha_crude_wrapped <- function(x,y) d_mahalanobis_crude(x, mu = mu_crude, Sigma = Sigma_crude)

X_clr    <- t(apply(X, 1, clr))
Sigma_clr <- cov(X_clr)
dist_maha_clr_wrapped <- function(x,y) d_mahalanobis_clr(x, y, Sigma_clr = Sigma_clr)

dist_list[["Mahalanobis (crude)"]] <- dist_maha_crude_wrapped
dist_list[["Mahalanobis (clr)"]]   <- dist_maha_clr_wrapped

## -------------------------
## 3) Run clustering for ALL methods
## -------------------------
al_results <- purrr::imap(dist_list, function(f, name) {
  kmeans_compositional(closed_3, numClusters = K, distFun = f, seed = 123)
})

## -------------------------
## 4) Build data frame with ALL cluster columns
## -------------------------
al_df <- closed_3
for (nm in names(al_results)) {
  colname <- paste0("cluster_", gsub("[^A-Za-z0-9]+", "_", nm))
  al_df[[colname]] <- factor(al_results[[nm]]$clusterAssignments)
}

## -------------------------
## 5) Ternary facet plot (ALL methods at once)
## -------------------------
cluster_cols <- grep("^cluster_", names(al_df), value = TRUE)

al_long <- al_df %>%
  pivot_longer(cols = all_of(cluster_cols),
               names_to = "Method",
               values_to = "Cluster") %>%
  mutate(
    Method = gsub("^cluster_", "", Method),
    Method = gsub("_", " ", Method),
    Cluster = factor(Cluster)
  )

# fixed palette for clusters (points)
K_here <- nlevels(al_long$Cluster)
pal <- brewer.pal(max(3, min(8, K_here)), "Set1")
pal <- setNames(rep(pal, length.out = K_here), levels(al_long$Cluster))

p_ternary_all <- ggtern(al_long, aes(RM, WM, E, color = Cluster)) +
  geom_point(alpha = 0.85, size = 2.8) +
  facet_wrap(~ Method, nrow = 3) +
  scale_color_manual(values = pal, name = "Cluster") +
  theme_rgbw() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(title = "Alimentation data: ternary clustering across all measures")

print(p_ternary_all)

## Optional save
# ggsave("Alimentation_Ternary_AllMethods.pdf", p_ternary_all, width = 14, height = 10)

## -------------------------
## 1) Build a long centroids table for ALL methods (ternary space)
## -------------------------
centroids_long <- purrr::imap_dfr(al_results, function(res, method) {
  C <- as.data.frame(res$centroids)
  colnames(C) <- c("RM","WM","E")
  C$Cluster <- factor(seq_len(nrow(C)))          # centroid id (1..K)
  C$Method  <- gsub("[^A-Za-z0-9]+", " ", method)
  C$lab <- sprintf("(%.2f, %.2f, %.2f)", C$RM, C$WM, C$E)
  C
})

## -------------------------
## 2) Ternary facet plot + centroids (black X) + centroid labels
## -------------------------


p_ternary_all_centroids <- ggtern(al_long, aes(RM, WM, E, color = Cluster)) +
  geom_point(alpha = 0.85, size = 2.8) +
  
  # centroids as black X (keep z for ggtern)
  geom_point(
    data = centroids_long,
    mapping = ggtern::aes(x = RM, y = WM, z = E),
    inherit.aes = FALSE,
    shape = 4, size = 3.6, stroke = 1.3,
    color = "black"
  ) +
  
  facet_wrap(~ Method, nrow = 3) +
  scale_color_manual(values = pal, name = "Cluster") +
  theme_rgbw() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(title = "")

suppressWarnings(print(p_ternary_all_centroids))
## Optional save
 ggsave("Alimentation_Ternary_AllMethods.pdf", p_ternary_all_centroids, width = 14, height = 10)
## -------------------------
## 6) CLR-PCA + facet overlays (ALL methods)
## -------------------------
clr_transformation <- function(x, eps = 1e-15) {
  x <- pmax(as.numeric(x), eps)
  gm <- exp(mean(log(x)))
  log(x/gm)
}

data_clr <- t(apply(as.matrix(al_df[,c("RM","WM","E")]), 1, clr_transformation))
colnames(data_clr) <- c("RM","WM","E")

pca_res <- prcomp(data_clr, center = FALSE, scale. = FALSE)

scores <- as.data.frame(pca_res$x[,1:2])
colnames(scores) <- c("PC1","PC2")

loadings <- as.data.frame(pca_res$rotation[,1:2])
loadings$variable <- rownames(loadings)

max_scale <- max(abs(scores$PC1), abs(scores$PC2)) * 0.9
arrows_data <- loadings %>%
  transmute(
    x=0, y=0,
    xend = PC1 * max_scale,
    yend = PC2 * max_scale,
    variable = variable
  )

scores_long <- purrr::map_dfr(names(al_results), function(m) {
  tibble(
    PC1 = scores$PC1,
    PC2 = scores$PC2,
    Cluster = factor(al_results[[m]]$clusterAssignments),
    Method = m
  )
}) %>%
  mutate(Method = gsub("[^A-Za-z0-9]+"," ", Method))

p_pca_all <- ggplot(scores_long, aes(PC1, PC2, color = Cluster)) +
  geom_point(alpha = 0.85, size = 2.2) +
  geom_segment(
    data = arrows_data,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(0.18, "inches")),
    alpha = 0.25, color = "gray40"
  ) +
  geom_text(
    data = arrows_data,
    aes(x = xend, y = yend, label = variable),
    inherit.aes = FALSE,
    size = 3.2, color = "black", fontface = "bold",
    hjust = 1.15, vjust = 1.15
  ) +
  facet_wrap(~ Method, nrow = 3) +
  scale_color_manual(values = pal, name = "Cluster") +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "Alimentation data: CLR-PCA with clustering overlays (all measures)",
    x = "PC1 (CLR-PCA)", y = "PC2 (CLR-PCA)"
  )

print(p_pca_all)
