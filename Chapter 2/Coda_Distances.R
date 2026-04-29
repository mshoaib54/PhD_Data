# ============================================================
# CoDA Distances & Divergences on the Probability Simplex (R)
# (Robust version: NO data.table ':=' syntax)
# ------------------------------------------------------------
# Provides:
#   - distances/divergences (as discussed)
#   - pairwise matrices + bubble heatmaps
#   - correlation matrix across measures
#   - ternary plots (n=3) of distance-to-reference for ALL measures
#     with per-panel mean shown, and both facet_wrap and patchwork wrap
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtern)
  library(reshape2)
  library(MASS)
  # optional (only used for wrap layout)
  # install.packages("patchwork")
  library(patchwork)
})

# -------------------------
# Utilities
# -------------------------

# Closure to constant kappa (default 1 => probability simplex)
close_to <- function(P, kappa = 1) {
  P <- as.matrix(P)
  rs <- rowSums(P)
  if (any(!is.finite(rs)) || any(rs <= 0)) stop("Row sum <= 0 found; cannot close.")
  P * (kappa / rs)
}

close_vec <- function(x, kappa = 1) {
  x <- as.numeric(x)
  s <- sum(x)
  if (!is.finite(s) || s <= 0) stop("Invalid sum; cannot close.")
  x * (kappa / s)
}

gm <- function(x) exp(mean(log(x)))

safe_log  <- function(x, eps = 1e-12) log(pmax(x, eps))
safe_acos <- function(u) acos(pmin(1, pmax(-1, u)))

add_eps_close <- function(x, eps = 1e-6, kappa = 1) {
  x <- as.numeric(x)
  close_vec(x + eps, kappa = kappa)
}

# -------------------------
# Log-ratio transforms
# -------------------------

alr <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (any(x <= 0)) stop("alr requires strictly positive components.")
  log(x[1:(n-1)] / x[n])
}

# clr(x) = log(x) - (1/n) sum_j log(x_j) * 1_n   (row-centered log vector)
clr <- function(x) {
  x <- as.numeric(x)
  if (any(x <= 0)) stop("clr requires strictly positive components.")
  lx <- log(x)
  lx - mean(lx)
}

# ilr(x) = log(x) Psi^T, where Psi is contrast matrix (n-1 x n)
ilr <- function(x, Psi) {
  x <- as.numeric(x)
  if (any(x <= 0)) stop("ilr requires strictly positive components.")
  as.numeric(log(x) %*% t(Psi))
}

# -------------------------
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

# -------------------------
# Pairwise distance matrices
# -------------------------

pairwise_matrix <- function(P, fun, ..., diag_zero = TRUE) {
  P <- as.matrix(P)
  m <- nrow(P)
  M <- matrix(NA_real_, m, m)
  for (i in 1:m) {
    for (j in 1:m) {
      if (diag_zero && i == j) {
        M[i, j] <- 0
      } else {
        out <- try(fun(P[i, ], P[j, ], ...), silent = TRUE)
        M[i, j] <- if (inherits(out, "try-error") || !is.finite(out)) NA_real_ else out
      }
    }
  }
  rownames(M) <- rownames(P)
  colnames(M) <- rownames(P)
  M
}

# -------------------------
# Measure-correlation matrix
# -------------------------

measure_corr_matrix <- function(P, measures_named) {
  P <- as.matrix(P)
  m <- nrow(P)

  get_vec <- function(fun) {
    vals <- c()
    for (i in 1:(m-1)) for (j in (i+1):m) {
      out <- try(fun(P[i,], P[j,]), silent = TRUE)
      vals <- c(vals, if (inherits(out, "try-error") || !is.finite(out)) NA_real_ else out)
    }
    vals
  }

  Dall <- sapply(measures_named, get_vec)
  cor(Dall, use = "complete.obs")
}

# -------------------------
# Heatmap (bubble) plot
# -------------------------

plot_heat_bubble <- function(M, fill_name = "Value", title = "") {
  df_long <- reshape2::melt(M)
  colnames(df_long) <- c("Row","Col","Val")

  ggplot(df_long, aes(x = factor(Col), y = factor(Row))) +
    geom_point(aes(size = Val, fill = Val),
               shape = 21, color = "black", stroke = 0.8, alpha = 1, na.rm = TRUE) +
    geom_label(aes(label = ifelse(is.na(Val), "NA", sprintf("%.3f", Val))),
               size = 3.0, fontface = "bold",
               label.size = 0.15, fill = "white", alpha = 0.90,
               vjust = -0.8, na.rm = TRUE) +
    scale_fill_gradient(low = "lightblue", high = "navy", name = fill_name) +
    scale_size(range = c(3, 14), name = fill_name) +
    coord_fixed() +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
      axis.text.y = element_text(face = "bold", size = 11),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 12),
      plot.title  = element_text(face = "bold", hjust = 0.5, size = 13)
    ) +
    labs(title = title, x = "", y = "")
}

# -------------------------
# Ternary (n=3) helpers
# -------------------------

# Distance-to-reference in long format (no data.table)
distance_to_ref_long <- function(P, x0, measures, measure_labels = NULL) {
  P <- as.matrix(P)
  if (ncol(P) != 3) stop("Ternary plot requires 3-part compositions.")
  colnames(P) <- c("x1","x2","x3")
  x0 <- close_vec(x0, 1)

  df <- as.data.frame(P)

  # compute each measure as a column
  for (nm in names(measures)) {
    f <- measures[[nm]]
    df[[nm]] <- apply(P, 1, function(x) {
      out <- try(f(x, x0), silent = TRUE)
      if (inherits(out, "try-error") || !is.finite(out)) NA_real_ else out
    })
  }

  long <- reshape2::melt(
    df,
    id.vars = c("x1","x2","x3"),
    variable.name = "measure",
    value.name = "value"
  )

  # replace NaN/Inf with NA
  long$value[!is.finite(long$value)] <- NA_real_

  # attach readable labels
  if (!is.null(measure_labels)) {
    # keep only labels we know about
    keep_levels <- intersect(names(measure_labels), unique(long$measure))
    long <- long[long$measure %in% keep_levels, , drop = FALSE]
    long$measure_label <- factor(long$measure, levels = keep_levels, labels = unname(measure_labels[keep_levels]))
  } else {
    long$measure_label <- long$measure
  }

  long
}

# Mean per measure (to annotate inside each ternary panel)
measure_means <- function(long_df, label_digits = 4, where = c(1/3, 1/3, 1/3)) {
  # aggregate over the plotted label (so it matches facets)
  tmp <- aggregate(value ~ measure_label, data = long_df, FUN = function(v) mean(v, na.rm = TRUE))
  names(tmp)[2] <- "mean_value"

  tmp$x1 <- where[1]; tmp$x2 <- where[2]; tmp$x3 <- where[3]
  tmp$lab <- sprintf(paste0("mean = %.", label_digits, "f"), tmp$mean_value)
  tmp
}

plot_ternary_facet_all <- function(long_df, means_df,
                                   title = "Distances to reference composition x0") {
  ggtern(long_df, aes(x = x1, y = x2, z = x3)) +
    geom_point(aes(color = value), size = 1.8, alpha = 0.9, na.rm = TRUE) +
    geom_label(
      data = means_df,
      aes(x = x1, y = x2, z = x3, label = lab),
      inherit.aes = FALSE,
      size = 3.1, fontface = "bold",
      label.size = 0.15, fill = "white", alpha = 0.9
    ) +
    facet_wrap(~ measure_label, scales = "free") +
    theme_bw() +
    labs(title = title, color = "")
}

plot_one_ternary <- function(df_one, mean_one, title = NULL) {
  ggtern(df_one, aes(x = x1, y = x2, z = x3)) +
    geom_point(aes(color = value), size = 2.0, alpha = 0.9, na.rm = TRUE) +
    geom_label(
      data = mean_one,
      aes(x = x1, y = x2, z = x3, label = lab),
      inherit.aes = FALSE,
      size = 3.1, fontface = "bold",
      label.size = 0.15, fill = "white", alpha = 0.9
    ) +
    theme_bw() +
    labs(title = title, color = "")
}

plot_ternary_wrap_all <- function(long_df, means_df, ncol = 3) {
  measures_list <- levels(long_df$measure_label)
  if (is.null(measures_list)) measures_list <- sort(unique(as.character(long_df$measure_label)))

  plots <- lapply(measures_list, function(m) {
    df_m <- subset(long_df, measure_label == m)
    mean_m <- subset(means_df, measure_label == m)
    plot_one_ternary(df_m, mean_m, title = m)
  })

  wrap_plots(plots, ncol = ncol)
}

# ============================================================
# EXAMPLE RUN (edit/replace with your own dataset)
# ============================================================

set.seed(1)
P_raw <- matrix(rexp(60, rate = 1), ncol = 3)
rownames(P_raw) <- paste0("obs", 1:nrow(P_raw))
P <- close_to(P_raw, 1)  # probability simplex

# Measures (internal IDs)
measures <- list(
  d_A       = function(x,y) d_A(x,y),
  d_ilr     = function(x,y) d_ilr(x,y),
  d_alr     = function(x,y) d_alr(x,y),
  d_angular = function(x,y) d_angular(x,y),
  d_bcos    = function(x,y) d_bhattacharyya_angular(x,y),
  d_b       = function(x,y) d_bhattacharyya(x,y),
  d_fidelity= function(x,y) d_fidelity(x,y),
  d_fisher  = function(x,y) d_fisher(x,y),
  d_hellinger=function(x,y) d_hellinger(x,y),
  d_ckl     = function(x,y) d_compositional_kl(x,y),
  d_kaniadakis = function(x,y) d_kaniadakis(x,y),
  d_cityblock=function(x,y) d_cityblock(x,y),
  d_tv      = function(x,y) d_total_variation(x,y),
  d_td      = function(x,y) d_triangular_discrimination(x,y),
  d_renyi05 = function(x,y) d_renyi(x,y, alpha = 0.5),
  d_za      = function(x,y) d_zero_adjusted_aitchison(x,y, eps = 1e-6)
)

# Full names (used in plots)
measure_labels <- c(
  d_A        = "Aitchison",
  d_ilr      = "Isometric log-ratio (ilr) ",
  d_alr      = "Additive log-ratio",
  d_angular  = "Angular ",
  d_bcos     = "Bhattacharyya angular",
  d_b        = "Bhattacharyya standard",
  d_fidelity = "Fidelity ",
  d_fisher   = "Fisher ",
  d_hellinger= "Hellinger ",
  d_ckl      = "Compositional KL",
  d_kaniadakis = "Kaniadakis divergence",
  d_cityblock= "City-block (L1) ",
  d_tv       = "Total variation",
  d_td       = "Triangular discrimination",
  d_renyi05  = "Rényi divergence",
  d_za       = "Zero-adjusted Aitchison"
)

# 1) Correlation matrix of measures (based on upper-triangular pairwise values)
Cmat <- measure_corr_matrix(P, measures)
colnames(Cmat) <- unname(measure_labels[colnames(Cmat)])
rownames(Cmat) <- unname(measure_labels[rownames(Cmat)])
print(Cmat)
p2 <- plot_heat_bubble(Cmat, fill_name = "Correlation", title = "Correlation of measures")

ggsave("CWPO.pdf", plot = p2, device = cairo_pdf, width = 12, height = 10)
# 2) Pairwise distance matrix for one measure (Aitchison)
M_A <- pairwise_matrix(P, measures$d_A)
print(plot_heat_bubble(M_A, fill_name = "Aitchison distance", title = "Pairwise Aitchison distance"))

# 3) Ternary plots of distance-to-reference for ALL measures (n=3)
x0 <- c(1/3, 1/3, 1/3)


# 3) Now your commands should work:
long_df  <- distance_to_ref_long(P, x0, measures, measure_labels = measure_labels)
means_df <- measure_means(long_df, label_digits = 4, where = c(1/3, 1/3, 1/3))

p_facet <- plot_ternary_facet_all(
  long_df, means_df,
  title = "All distances/divergences to x0 (mean shown per panel)"
)
print(p_facet)
ggsave("CWPOT.pdf", plot = p_facet, device = cairo_pdf, width = 12, height = 10)
p_wrap <- plot_ternary_wrap_all(long_df, means_df, ncol = 3)
print(p_wrap)

dir.create("plots", showWarnings = FALSE)



# ============================================================
# End of script
# ============================================================

setwd("/Users/muhammadshoaib/Desktop/PhD Data/Thesis/Thesis_data/")

st_data <- load("statisticianmebudget.RData")
timebud <- get(st_data[1])

# ============================================================
# REAL DATA RUN: Statistician time-budget data (20 days)
# ============================================================

# 0) Optional: Cairo PDF device (recommended)
# install.packages("Cairo")
 library(Cairo)

# 1) Enter (or read) the real dataset -----------------------------------------



# Name rows by day for nicer heatmaps
rownames(timebud) <- paste0("Day", timebud$Day)

# 2) Build 6-part composition on the probability simplex ----------------------
# Parts: (T, C, A, R, O, S)
P6_raw <- as.matrix(timebud[, c("T","C","A","R","O","S")])
P6 <- close_to(P6_raw, 1)                 # close each day to sum=1
rownames(P6) <- rownames(timebud)

# 3) Correlation matrix across measures (on full 6-part compositions) ---------
# (Uses your measure_corr_matrix and measures list)
Cmat_real <- measure_corr_matrix(P6, measures)
colnames(Cmat_real) <- unname(measure_labels[colnames(Cmat_real)])
rownames(Cmat_real) <- unname(measure_labels[rownames(Cmat_real)])

p_corr_real <- plot_heat_bubble(
  Cmat_real,
  fill_name = "Correlation",
  title = "Time-budget data: correlation across measures"
)
print(p_corr_real)

# Save correlation plot
ggsave("TimeBudget_Correlation.pdf", plot = p_corr_real,
       device = cairo_pdf, width = 12, height = 10)

# 4) Pairwise matrix + heatmap for selected measures --------------------------
# Example: Aitchison, Compositional KL, Zero-adjusted Aitchison, Additive log-ratio, Kaniadakis

M_A   <- pairwise_matrix(P6, measures$d_A)
M_CKL <- pairwise_matrix(P6, measures$d_ckl)
M_ZA  <- pairwise_matrix(P6, measures$d_za)
M_ALR <- pairwise_matrix(P6, measures$d_alr)
M_KAN <- pairwise_matrix(P6, measures$d_kaniadakis)

pA   <- plot_heat_bubble(M_A,   fill_name="Aitchison",        title="Pairwise Aitchison (20 days)")
pCKL <- plot_heat_bubble(M_CKL, fill_name="Compositional KL", title="Pairwise Compositional KL (20 days)")
pZA  <- plot_heat_bubble(M_ZA,  fill_name="Zero-adjusted A.", title="Pairwise Zero-adjusted Aitchison (20 days)")
pALR <- plot_heat_bubble(M_ALR, fill_name="Additive log-ratio",title="Pairwise Additive log-ratio (20 days)")
pKAN <- plot_heat_bubble(M_KAN, fill_name="Kaniadakis",       title="Pairwise Kaniadakis divergence (20 days)")

print(pA);   ggsave("TimeBudget_Pairwise_Aitchison.pdf", plot=pA,   device=cairo_pdf, width=12, height=7)
print(pCKL); ggsave("TimeBudget_Pairwise_CKL.pdf",      plot=pCKL, device=cairo_pdf, width=12, height=7)
print(pZA);  ggsave("TimeBudget_Pairwise_ZAitchison.pdf",plot=pZA, device=cairo_pdf, width=12, height=7)
print(pALR); ggsave("TimeBudget_Pairwise_ALR.pdf",      plot=pALR, device=cairo_pdf, width=12, height=7)
print(pKAN); ggsave("TimeBudget_Pairwise_Kaniadakis.pdf",plot=pKAN,device=cairo_pdf, width=12, height=7)

# 5) Ternary visualization: (Work, O, S) -------------------------------------
# Build a 3-part composition: Work = T+C+A+R, plus O, plus S
work_raw <- rowSums(timebud[, c("T","C","A","R")])
P3_raw <- cbind(Work = work_raw, O = timebud$O, S = timebud$S)
P3 <- close_to(P3_raw, 1)
colnames(P3) <- c("x1","x2","x3")  # required by your ternary helpers (x1,x2,x3)
rownames(P3) <- rownames(timebud)

# Reference point on ternary simplex (barycenter)
x0 <- c(1/3, 1/3, 1/3)

# Long format + means (mean printed in each facet panel)
long_df  <- distance_to_ref_long(P3, x0, measures, measure_labels = measure_labels)
means_df <- measure_means(long_df, label_digits = 4, where = c(1/3, 1/3, 1/3))

# IMPORTANT: avoid the "Ignoring unknown aesthetics: z" warning
# by using T/L/R aesthetics for labels (ggtern-native aesthetics).
plot_ternary_facet_all_TLR <- function(long_df, means_df,
                                       title = "Distances/divergences to x0 on (Work, O, S)") {
  ggtern(long_df, aes(T = x1, L = x2, R = x3)) +
    geom_point(aes(color = value), size = 1.9, alpha = 0.9, na.rm = TRUE) +
    ggtern::geom_text(   # ggtern geom, understands ternary aesthetics
      data = means_df,
      aes(T = x1, L = x2, R = x3, label = lab),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 3.6
    ) +
    facet_wrap(~ measure_label, scales = "free") +
    theme_bw() +
    labs(title = title, color = "")
}

p_facet <- plot_ternary_facet_all_TLR(
  long_df, means_df,
  title = "Time-budget data: distances/divergences to x0 on (Work, O, S)\n(mean shown per panel)"
)
print(p_facet)

ggsave("TimeBudget_Ternary_Facets_WorkO_S.pdf", plot = p_facet,
       device = cairo_pdf, width = 14, height = 10)

# (Optional) Patchwork wrap instead of facets
p_wrap <- plot_ternary_wrap_all(long_df, means_df, ncol = 3)
print(p_wrap)
ggsave("TimeBudget_Ternary_Wrap_WorkO_S.pdf", plot = p_wrap,
       device = cairo_pdf, width = 14, height = 10)



