############################################################
## FROM SCRATCH R SCRIPT (Year-wise analysis + all plots)
## Kaniadakis divergence + mixture/exponential displacements
## + exponential-coordinate barycenter + deviation
############################################################

## =========================
## 0) Setup
## =========================
setwd("/Users/muhammadshoaib/Desktop/IG notes")   # <-- change if needed

library(showtext)
library(sysfonts)
library(readxl)     # read Excel
library(dplyr)      # data wrangling
library(tidyr)      # pivot_longer
library(ggplot2)    # plotting
library(reshape2)   # melt for heatmaps

## =========================
## 1) Load data
## =========================
data_raw <- read_excel("CoDa.xlsx")                                # read your file

countries <- c("Belgium","Denmark","France","Germany","Greece",
               "Italy","Netherland","Spain","Switzerland","UK")    # column order

library(dplyr)

data <- data_raw %>%
  mutate(Year = as.integer(Year)) %>%
  dplyr::select(Year, dplyr::all_of(countries))
## Convert comma decimals to numeric (safe for Excel files that store as text)
data[countries] <- lapply(data[countries], function(x) {
  as.numeric(gsub(",", ".", as.character(x)))
})

## =========================
## 2) Create compositional data (proportions)
## =========================
rowsum <- rowSums(data[countries])                                 # row totals
com_data <- data                                                   # copy
com_data[countries] <- data[countries] / rowsum                    # proportions per year

years <- com_data$Year                                             # year labels
P <- as.matrix(com_data[, countries])                              # matrix of compositions (n x d)
n <- nrow(P)                                                       # number of years
d <- ncol(P)                                                       # number of countries

## Long data for composition plots
com_long <- com_data %>%
  pivot_longer(cols = all_of(countries),
               names_to = "Country",
               values_to = "Value")

## =========================
## 3) Define ALL functions (paper-consistent)
## =========================

## ---- normalize vector to composition (for safety) ----
comp_normalize <- function(x, eps = 1e-12) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  x <- pmax(x, 0)
  s <- sum(x)
  if (s <= 0) return(rep(1/length(x), length(x)))
  x <- x / s
  x <- pmax(x, eps)
  x / sum(x)
}

## ---- Kaniadakis ingredients ----
A <- function(p) 2*(p^2)/(1+p^2)                                   # Nudit’s function

ptilde <- function(p) {                                            # escort distribution
  ap <- A(p)
  ap / sum(ap)
}

Ek <- function(p, u) sum(u * ptilde(p))                            # escort expectation eE_p[u]

Kl <- function(p, eps = 1e-12) {                                   # Kaniadakis log (your definition)
  p <- pmax(as.numeric(p), eps)
  0.5 * (p - 1/p)
}

Kexp <- function(x) {                                              # inverse of Kl: kexp
  x <- as.numeric(x)
  x + sqrt(x^2 + 1)
}

Kdiv <- function(p, q) {                                           # Kaniadakis divergence
  Ek(p, Kl(p)) - Ek(p, Kl(q))
}

## =========================
## Displacements (Eq. 20 and Eq. 21)
## =========================
mixture_disp <- function(p, q) {
  (q - p) / A(p)
}

exponential_disp <- function(p, q) {
  u <- Kl(q) - Kl(p)
  u - Ek(p, u)
}

## =========================
## 2) Mixture barycenter + deviation (affine mean in eta-coordinates)
##    eta_p(q) = (q-p)/(A∘p)  =>  q = p + (A∘p) ∘ eta
## =========================
Kbary_dev_mixture <- function(P, years = NULL, p_ref = NULL,
                             deviation = c("Kbar_to_f","Kf_to_bar"),
                             eps = 1e-12) {

  deviation <- match.arg(deviation)
  P <- as.matrix(P)
  n <- nrow(P)

  for (j in 1:n) P[j,] <- comp_normalize(P[j,], eps = eps)

  if (is.null(p_ref)) p_ref <- P[1,]
  p_ref <- comp_normalize(p_ref, eps = eps)

  ## compute eta_p(f_j)
  Eta <- t(apply(P, 1, function(q) mixture_disp(p_ref, q)))

  ## mean affine coordinate
  etabar <- colMeans(Eta)

  ## reconstruct barycenter in the mixture chart: fbar = p + A(p) ∘ etabar
  fbar_raw <- p_ref + A(p_ref) * etabar
  fbar <- comp_normalize(fbar_raw, eps = eps)

  ## deviations (1D summaries)
  dev <- numeric(n)
  if (deviation == "Kbar_to_f") {
    for (j in 1:n) dev[j] <- Kdiv(fbar, P[j,])   # DK(fbar || f_j)
  } else {
    for (j in 1:n) dev[j] <- Kdiv(P[j,], fbar)   # DK(f_j || fbar)
  }

  list(
    barycenter = fbar,
    mean_coord = etabar,
    deviation  = data.frame(Year = years, Deviation = dev),
    reference  = p_ref,
    deviation_type = deviation
  )
}

## =========================
## 3) Exponential barycenter + deviation (your LaTeX section)
##    fbar ∝ kexp( sbar + klog(p_ref) ) then normalize
## =========================
Kbary_dev_exponential <- function(P, years = NULL, p_ref = NULL,
                                 deviation = c("Kbar_to_f","Kf_to_bar"),
                                 eps = 1e-12) {

  deviation <- match.arg(deviation)
  P <- as.matrix(P)
  n <- nrow(P)

  for (j in 1:n) P[j,] <- comp_normalize(P[j,], eps = eps)

  if (is.null(p_ref)) p_ref <- P[1,]
  p_ref <- comp_normalize(p_ref, eps = eps)

  ## s_p(f_j)
  S <- t(apply(P, 1, function(q) exponential_disp(p_ref, q)))

  ## mean affine coordinate
  sbar <- colMeans(S)

  ## reconstruct barycenter: fbar ∝ kexp(sbar + klog(p_ref))
  fbar <- comp_normalize(Kexp(sbar + Kl(p_ref, eps = eps)), eps = eps)

  ## deviations
  dev <- numeric(n)
  if (deviation == "Kbar_to_f") {
    for (j in 1:n) dev[j] <- Kdiv(fbar, P[j,])
  } else {
    for (j in 1:n) dev[j] <- Kdiv(P[j,], fbar)
  }

  list(
    barycenter = fbar,
    mean_coord = sbar,
    deviation  = data.frame(Year = years, Deviation = dev),
    reference  = p_ref,
    deviation_type = deviation
  )
}
## =========================
## 4) VISUALIZATIONS (Year-wise)
## =========================

## ---- Plot 1: Stacked composition over years ----
p1 <- ggplot(com_long, aes(x = Year, y = Value, fill = Country)) +
  geom_area() +
  theme_minimal() +
  labs(title = "Composition over years (stacked)",
       x = "Year", y = "Proportion")
print(p1)

## ---- Plot 2: Each country separately (easy to read) -

font_add("Times New Roman", regular = "/Library/Fonts/Times New Roman.ttf")
showtext_auto(TRUE)

p2 <- ggplot(com_long |> dplyr::filter(!is.na(Year), !is.na(Value))) +
  geom_line(aes(x = Year, y = Value, color = Country), linewidth = 1.2) +
  geom_point(aes(x = Year, y = Value, color = Country), size = 2) +
  facet_wrap(~ Country, nrow = 2, scales = "free_y") +
  theme_bw() +
  labs(title = "Country Wise Proportions Over Time",
       x = "", y = "", color = "") +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    
    ## title: bold + Times New Roman
    plot.title = element_text(family = "Times New Roman", face = "bold", size = 16, hjust = 0.5),
    
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  guides(color = guide_legend(nrow = 1))

print(p2)

ggsave("CWPOT.pdf", plot = p2, device = cairo_pdf, width = 12, height = 6)
## ---- Plot 3: Heatmap of Kdiv between all pairs of years ----
Kmat <- matrix(0, n, n)
for (i in 1:n) for (j in 1:n) Kmat[i, j] <- Kdiv(P[i,], P[j,])

rownames(Kmat) <- years
colnames(Kmat) <- years

K_long <- melt(Kmat)
colnames(K_long) <- c("Year_i","Year_j","Kdiv")

p3d_allvals <- ggplot(K_long, aes(x = factor(Year_j), y = factor(Year_i))) +
  
  ## circles
  geom_point(aes(size = Kdiv, fill = Kdiv),
             shape = 21, color = "black", stroke = 0.8, alpha = 1) +
  
  ## labels
  geom_label(aes(label = sprintf("%.3f", Kdiv)),
             size = 3.6, fontface = "bold",
             label.size = 0.15, fill = "white", alpha = 0.90,
             vjust = -0.8) +
  
  ## nicer legend names
  scale_fill_gradient(low = "lightblue", high = "navy",
                      name = "Kaniadakis divergence") +
  scale_size(range = c(3, 14),
             name = "Kaniadakis divergence") +
  
  coord_fixed() +
  theme_bw(base_size = 12) +
  theme(
    ## make year tick labels MORE bold + bigger
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    
    ## (optional) make axis titles bold too
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
  ) +
  labs(
    title = "",
    x = "", y = ""
  )

print(p3d_allvals)
ggsave("kdiv.pdf", plot = p3d_allvals, device = cairo_pdf, width = 10, height = 10)

## ---- Plot 4: Drift from baseline year (first year, usually 2008) ----
## ---- Drift from baseline ----
baseline_year <- min(years)
p0 <- P[which(years == baseline_year), ]

K_from_base <- sapply(1:n, function(i) Kdiv(p0, P[i,]))

## ---- Year-to-year change ----
K_prev <- rep(NA_real_, n)
for (i in 2:n) K_prev[i] <- Kdiv(P[i-1,], P[i,])

## ---- Combine for one faceted plot ----
df_comp <- dplyr::bind_rows(
  data.frame(Year = years, Value = K_from_base, Type = paste0("Drift: Kdiv(", baseline_year, ", Year)")),
  data.frame(Year = years, Value = K_prev,      Type = "Year-to-year: Kdiv(Year[t-1], Year[t])")
) %>%
  dplyr::filter(!is.na(Value)) %>%
  dplyr::mutate(Type = factor(Type))   # keep stable ordering

p45 <- ggplot(df_comp, aes(x = Year, y = Value, color = Type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.2) +
  facet_wrap(~ Type, nrow = 1, scales = "free_y") +   # 2 graphs in one window
  scale_x_continuous(breaks = sort(unique(years))) +  # show each year on x-axis
  theme_bw(base_size = 12) +
  theme(
    strip.text   = element_text(face = "bold", size = 12),
    axis.text.x  = element_text(face = "bold", angle = 45, hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title   = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"                         # legend below
  ) +
  labs(
    title = "",
    x = "",
    y = "",
    color = ""
  ) +
  guides(color = guide_legend(nrow = 1))

print(p45)
ggsave(
  filename = "Kdiv_drift_and_yearly_change.pdf",
  plot     = p45,
  device   = cairo_pdf,
  width    = 12,
  height   = 6
)

## ---- Plot 5: Year-to-year change (consecutive years) ----
K_prev <- rep(NA_real_, n)
for (i in 2:n) K_prev[i] <- Kdiv(P[i-1,], P[i,])

df5 <- data.frame(Year = years, Kdiv_prev_year = K_prev)

p5 <- ggplot(df5, aes(x = Year, y = Kdiv_prev_year)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Year-to-year change: Kdiv(Year[t-1], Year[t])",
       x = "Year", y = "Kaniadakis divergence")
print(p5)

## ---- Plot 6: Mixture displacement heatmap η_baseline(Year) (Eq. 20) ----
## Mixture displacement matrix -> long
Mix_mat <- t(apply(P, 1, function(q) mixture_disp(p0, q)))
rownames(Mix_mat) <- years
colnames(Mix_mat) <- countries

Mix_long <- reshape2::melt(Mix_mat)
colnames(Mix_long) <- c("Year","Country","Value")
Mix_long$Type <- paste0("Mixture  Displacement ", baseline_year, "(Year)")

## Exponential displacement matrix -> long
Exp_mat <- t(apply(P, 1, function(q) exponential_disp(p0, q)))
rownames(Exp_mat) <- years
colnames(Exp_mat) <- countries

Exp_long <- reshape2::melt(Exp_mat)
colnames(Exp_long) <- c("Year","Country","Value")
Exp_long$Type <- paste0("Exponential  Displacement ", baseline_year, "(Year)")

## Combine
Disp_long <- dplyr::bind_rows(Mix_long, Exp_long)

## optional: fix ordering
Disp_long$Country <- factor(Disp_long$Country, levels = countries)
Disp_long$Year    <- factor(Disp_long$Year, levels = rev(sort(unique(as.integer(as.character(Disp_long$Year))))))

p_disp <- ggplot(Disp_long, aes(x = Country, y = Year, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", Value)), size = 3.2, fontface = "bold") +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0,
                       name = "") +
  facet_wrap(~ Type, nrow = 1) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    strip.text  = element_text(face = "bold"),
    plot.title  = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  ) +
  labs(title = " ",
       x = "", y = "")

print(p_disp)

ggsave("Displacements_mixture_exponential.pdf", p_disp, device = cairo_pdf, width = 14, height = 6)
## Register Times New Roman (macOS)
font_add("Times New Roman", regular = "/Library/Fonts/Times New Roman.ttf")
showtext_auto(TRUE)

p_norm <- ggplot(df_norm, aes(x = Year, y = Norm, color = Type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.2) +
  scale_x_continuous(breaks = sort(unique(years))) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x  = element_text(face = "bold", angle = 45, hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    
    ## Title in Times New Roman + bold
    plot.title   = element_text(family = "Times New Roman", face = "bold", size = 16, hjust = 0.5),
    
    legend.position = "bottom"
  ) +
  labs(
    title = paste0("Magnitude of displacements relative to ", baseline_year),
    x = "",
    y = "",
    color = ""
  ) +
  guides(color = guide_legend(nrow = 1))

print(p_norm)

## Save using cairo_pdf (important)
ggsave("Displacement_norms.pdf", p_norm, device = cairo_pdf, width = 12, height = 5)


library(sysfonts)
library(showtext)

## Register Times New Roman (macOS)
font_add("Times New Roman", regular = "/Library/Fonts/Times New Roman.ttf")
showtext_auto(TRUE)

p_norm <- ggplot(df_norm, aes(x = Year, y = Norm, color = Type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.2) +
  scale_x_continuous(breaks = sort(unique(years))) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x  = element_text(face = "bold", angle = 45, hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    
    ## Title in Times New Roman + bold
    plot.title   = element_text(family = "Times New Roman", face = "bold", size = 16, hjust = 0.5),
    
    legend.position = "bottom"
  ) +
  labs(
    title = paste0("Magnitude of displacements relative to ", baseline_year),
    x = "Year",
    y = "L2 norm",
    color = "Displacement"
  ) +
  guides(color = guide_legend(nrow = 1))

print(p_norm)

## Save using cairo_pdf (important)
#ggsave("Displacement_norms.pdf", p_norm, device = cairo_pdf, width = 12, height = 5)

P <- as.matrix(com_data[, countries])
years <- com_data$Year

mix_res <- Kbary_dev_mixture(P, years = years, deviation = "Kbar_to_f")
exp_res <- Kbary_dev_exponential(P, years = years, deviation = "Kbar_to_f")

mix_res$barycenter   # mixture barycenter
exp_res$barycenter   # exponential barycenter

mix_res$deviation    # DK(bary || Year) for mixture barycenter
exp_res$deviation    # DK(bary || Year) for exponential barycenter
df_dev <- dplyr::bind_rows(
  transform(mix_res$deviation, Method = "Mixture barycenter"),
  transform(exp_res$deviation, Method = "Exponential barycenter")
)

p_dev <- ggplot(df_dev, aes(x = Year, y = Deviation, color = Method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.2) +
  scale_x_continuous(breaks = sort(unique(years))) +
  theme_bw() +
  theme(
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    plot.title  = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  ) +
  labs(title = "Deviation from barycenter (DK(barycenter || Year))",
       x = "Year", y = "Kaniadakis divergence", color = "")

print(p_dev)



## ---- Optional: Plot barycenter composition itself ----
df_bary <- data.frame(Country = countries, Value = res_bary$barycenter)

p11 <- ggplot(df_bary, aes(x = Country, y = Value)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Barycenter composition (exponential-coordinate barycenter)",
       x = "", y = "Proportion")
print(p11)

############################################################
## END OF SCRIPT
############################################################