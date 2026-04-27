library("stagedtrees")
library("infotheo")
source("divergences.R")
conditional_information <- function (X, Y, S = NULL, method = "emp") 
{
  if (is.null(S)) {
    Ires <- mutinformation(X, Y, method)
  }
  else {
    U <- data.frame(S, X, Y)
    Hysx <- entropy(U, method)
    Hsx <- entropy(U[, c(1:NCOL(S), NCOL(S)+1)], method)
    Hys <- entropy(U[, c(1:NCOL(S), NCOL(U))], method)
    Hs <- entropy(S, method)
    Ires <- Hys - Hs - Hysx + Hsx
  }
  if(Ires < 0 & Ires > -0.001) {
    Ires <- 0
  }
  Ires
}


## maximising conditional mutual information with class
ordering_cmi <- function(train){
  n <- ncol(train) - 1
  free <- names(train)[-1]
  ms <- sapply(free, function(v){
    infotheo::mutinformation(X = train$answer, Y = train[[v]]) / infotheo::entropy(train[[v]])
  })
  selected <- names(which.max(ms))
  free <- free[-which.max(ms)]
  for (i in 2:n){
    ms <- sapply(free, function(v){
      conditional_information(X = train$answer, 
                              Y = train[[v]],
                              S = train[, selected]) / infotheo::entropy(train[[v]])
    })
    selected <- c(selected, names(which.max(ms)))
    free <- free[-which.max(ms)]
  }
  return(c("answer", selected))
}


## minimizing conditional entropy
ordering_ch <- function(train){
  n <- ncol(train) - 1
  free <- names(train)[-1]
  ms <- sapply(free, function(v){
    infotheo::condentropy(X = train[[v]], Y = train$answer)
  })
  selected <- c("answer", names(which.min(ms)))
  free <- free[-which.min(ms)]
  for (i in 2:n){
    ms <- sapply(free, function(v){
      infotheo::condentropy(X = train[[v]], 
                            Y = train[, selected])
    })
    selected <- c(selected, names(which.min(ms)))
    free <- free[-which.min(ms)]
    
  }
  return(selected)
}

## maximising normalized mi
ordering_mi <- function(train){
  ms <- sapply(names(train)[-1], function(v){
    infotheo::mutinformation(X = train$answer, Y = train[[v]]) / 
      infotheo::entropy(train[[v]])
  })
  names(train)[c(1, 1 + order(ms, decreasing = TRUE))]
}

st_full <- function(train, test, ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE)
  predict(model, newdata = test, class = "answer")
}


st_naive_totvar_cmi <- function(train, test, ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  model <- stages_hclust(model, k = 2, distance = "totvar", method = "ward.D2")
  predict(model, newdata = test, class = "answer")
}

st_naive_jefferys_cmi <- function(train, test,  ...) {
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  model <- stages_hclust(model, k = 2, distance = jeffreys, method = "ward.D2")
  predict(model, newdata = test, class = "answer")
}
st_naive_kaniadakis_cmi <- function(train, test, ...) {
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  model <- stages_hclust(model, k = 2, distance = kaniadakis, method = "ward.D2")
  predict(model, newdata = test, class = "answer")
}
st_naive_jensen_shannon_cmi <- function(train, test,  ...) {
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  model <- stages_hclust(model, k = 2, distance = jensenshannon, method = "ward.D2")
  predict(model, newdata = test, class = "answer")
}
st_naive_hellinger_cmi <- function(train, test,  ...) {
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  model <- stages_hclust(model, k = 2, distance = "hellinger", method = "ward.D2")
  predict(model, newdata = test, class = "answer")
}


st_full_mi <- function(train, test, ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_mi(train))
  predict(model, newdata = test, class = "answer")
}

st_full_ch <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_ch(train))
  predict(model, newdata = test, class = "answer")
}

st_full_cmi <- function(train, test, ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  predict(model, newdata = test, class = "answer")
}

st_indep <- function(train, test, ...){
  model <- indep(train, lambda = 1, join_unobserved = TRUE)
  predict(model, newdata = test, class = "answer")
}

st_indep_mi <- function(train, test,  ...){
  model <- indep(train, lambda = 1, join_unobserved = TRUE, order = ordering_mi(train))
  predict(model, newdata = test, class = "answer")
}

st_indep_ch <- function(train, test,  ...){
  model <- indep(train, lambda = 1, join_unobserved = TRUE, order = ordering_ch(train))
  predict(model, newdata = test, class = "answer")
}

st_indep_cmi <- function(train, test, ...){
  model <- indep(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  predict(model, newdata = test, class = "answer")
}

st_hc_indep_5 <- function(train, test, ...){
  order_var <- colnames(train)
  model <- indep(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 5) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 5) { order_var <- order_var[2:5] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_hc_indep_mi_5 <- function(train, test, ...){
  order_var <- ordering_mi(train)
  model <- indep(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 5) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 5) { order_var <- order_var[2:5] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_hc_indep_ch_5 <- function(train, test,  ...){
  order_var <- ordering_ch(train)
  model <- indep(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 5) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 5) { order_var <- order_var[2:5] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_hc_indep_cmi_5 <- function(train, test,  ...){
  order_var <- ordering_cmi(train)
  model <- indep(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 5) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 5) { order_var <- order_var[2:5] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_hc_indep_7 <- function(train, test,  ...){
  order_var <- colnames(train)
  model <- indep(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 7) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 7) { order_var <- order_var[2:7] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_hc_indep_mi_7 <- function(train, test,  ...){
  order_var <- ordering_mi(train)
  model <- indep(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 7) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 7) { order_var <- order_var[2:7] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_hc_indep_ch_7 <- function(train, test,  ...){
  order_var <- ordering_ch(train)
  model <- indep(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 7) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 7) { order_var <- order_var[2:7] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_hc_indep_cmi_7 <- function(train, test,  ...){
  order_var <- ordering_cmi(train)
  model <- indep(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 7) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 7) { order_var <- order_var[2:7] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_hc_full_5 <- function(train, test,  ...){
  order_var <- colnames(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 5) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 5) { order_var <- order_var[2:5] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_hc_full_mi_5 <- function(train, test,  ...){
  order_var <- ordering_mi(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 5) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 5) { order_var <- order_var[2:5] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_hc_full_cmi_5 <- function(train, test,  ...){
  order_var <- ordering_cmi(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 5) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 5) { order_var <- order_var[2:5] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_hc_full_ch_5 <- function(train, test,  ...){
  order_var <- ordering_ch(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 5) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 5) { order_var <- order_var[2:5] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_hc_full_7 <- function(train, test,  ...){
  order_var <- colnames(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 7) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 7) { order_var <- order_var[2:7] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_hc_full_cmi_7 <- function(train, test,  ...){
  order_var <- ordering_cmi(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 7) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 7) { order_var <- order_var[2:7] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_hc_full_mi_7 <- function(train, test,  ...){
  order_var <- ordering_mi(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 7) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 7) { order_var <- order_var[2:7] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_hc_full_ch_7 <- function(train, test,  ...){
  order_var <- ordering_ch(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 7) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 7) { order_var <- order_var[2:7] }
  model <- stages_hc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_fbhc <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE)
  model <- stages_fbhc(model)
  predict(model, newdata = test, class = "answer")
}

st_fbhc_ch <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_ch(train))
  model <- stages_fbhc(model)
  predict(model, newdata = test, class = "answer")
}

st_fbhc_cmi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  model <- stages_fbhc(model)
  predict(model, newdata = test, class = "answer")
}

st_fbhc_mi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_mi(train))
  model <- stages_fbhc(model)
  predict(model, newdata = test, class = "answer")
}

st_bhc_5 <- function(train, test,  ...){
  order_var <- colnames(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 5) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 5) { order_var <- order_var[2:5] }
  model <- stages_bhc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_bhc_mi_5 <- function(train, test,  ...){
  order_var <- ordering_mi(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 5) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 5) { order_var <- order_var[2:5] }
  model <- stages_bhc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}


st_bhc_cmi <- function(train, test,  ...){
  order_var <- ordering_cmi(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  model <- stages_bhc(model)
  predict(model, newdata = test, class = "answer")
}

st_bhc_cmi_5 <- function(train, test,  ...){
  order_var <- ordering_cmi(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 5) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 5) { order_var <- order_var[2:5] }
  model <- stages_bhc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_bhc_ch_5 <- function(train, test,  ...){
  order_var <- ordering_ch(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 5) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 5) { order_var <- order_var[2:5] }
  model <- stages_bhc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_bhc_7 <- function(train, test,  ...){
  order_var <- colnames(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 7) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 7) { order_var <- order_var[2:7] }
  model <- stages_bhc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_bhc_mi_7 <- function(train, test,  ...){
  order_var <- ordering_mi(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 7) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 7) { order_var <- order_var[2:7] }
  model <- stages_bhc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_bhc_cmi_7 <- function(train, test,  ...){
  order_var <- ordering_cmi(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 7) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 7) { order_var <- order_var[2:7] }
  model <- stages_bhc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_bhc_ch_7 <- function(train, test,  ...){
  order_var <- ordering_ch(train)
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = order_var)
  if(length(order_var) <= 7) order_var <- order_var[2:NCOL(train)]
  if(length(order_var) > 7) { order_var <- order_var[2:7] }
  model <- stages_bhc(model, scope = order_var)
  predict(model, newdata = test, class = "answer")
}

st_bj_kl <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE)
  model <- stages_bj(model, distance = "kullback", thr = 0.2)
  predict(model, newdata = test, class = "answer")
}

st_bj_kl_mi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_mi(train))
  model <- stages_bj(model, distance = "kullback", thr = 0.2)
  predict(model, newdata = test, class = "answer")
}

st_bj_kl_01_cmi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  model <- stages_bj(model, distance = "kullback", thr = 0.01)
  predict(model, newdata = test, class = "answer")
}

st_bj_kl_05_cmi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  model <- stages_bj(model, distance = "kullback", thr = 0.05)
  predict(model, newdata = test, class = "answer")
}

st_bj_kl_20_cmi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  model <- stages_bj(model, distance = "kullback", thr = 0.2)
  predict(model, newdata = test, class = "answer")
}

st_bj_kl_cmi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  model <- stages_bj(model, distance = "kullback", thr = 0.2)
  predict(model, newdata = test, class = "answer")
}

st_bj_kl_ch <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_ch(train))
  model <- stages_bj(model, distance = "kullback", thr = 0.2)
  predict(model, newdata = test, class = "answer")
}

st_bj_tv <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE)
  model <- stages_bj(model, distance = "totvar", thr = 0.2)
  predict(model, newdata = test, class = "answer")
}

st_bj_tv_mi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_mi(train))
  model <- stages_bj(model, distance = "totvar", thr = 0.2)
  predict(model, newdata = test, class = "answer")
}

st_bj_tv_cmi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  model <- stages_bj(model, distance = "totvar", thr = 0.2)
  predict(model, newdata = test, class = "answer")
}

st_bj_tv_ch <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_ch(train))
  model <- stages_bj(model, distance = "totvar", thr = 0.2)
  predict(model, newdata = test, class = "answer")
}

st_bj_cd <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE)
  model <- stages_bj(model, distance = "chandarw", thr = 0.2)
  predict(model, newdata = test, class = "answer")
}

st_bj_cd_mi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_mi(train))
  model <- stages_bj(model, distance = "chandarw", thr = 0.2)
  predict(model, newdata = test, class = "answer")
}

st_bj_cd_cmi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  model <- stages_bj(model, distance = "chandarw", thr = 0.2)
  predict(model, newdata = test, class = "answer")
}

st_bj_cd_ch <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_ch(train))
  model <- stages_bj(model, distance = "chandarw", thr = 0.2)
  predict(model, newdata = test, class = "answer")
}

st_hclust <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE)
  model <- stages_hclust(model, k = 2, distance = "totvar", method = "complete")
  predict(model, newdata = test, class = "answer")
}

st_hclust_mi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_mi(train))
  model <- stages_hclust(model, k = 2, distance = "totvar", method = "complete")
  predict(model, newdata = test, class = "answer")
}

st_hclust_cmi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  model <- stages_hclust(model, k = 2, distance = "totvar", method = "complete")
  predict(model, newdata = test, class = "answer")
}


st_hclust_wd2_cmi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  model <- stages_hclust(model, k = 2, distance = "totvar", method = "ward.D2")
  predict(model, newdata = test, class = "answer")
}

st_hclust_ch <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_ch(train))
  model <- stages_hclust(model, k = 2, distance = "totvar", method = "complete")
  predict(model, newdata = test, class = "answer")
}

st_kmeans <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE)
  model <- stages_kmeans(model, k = 2, algorithm = "Hartigan-Wong")
  predict(model, newdata = test, class = "answer")
}

st_kmeans_mi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_mi(train))
  model <- stages_kmeans(model, k = 2, algorithm = "Hartigan-Wong")
  predict(model, newdata = test, class = "answer")
}

st_kmeans_cmi <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_cmi(train))
  model <- stages_kmeans(model, k = 2, algorithm = "Hartigan-Wong")
  predict(model, newdata = test, class = "answer")
}

st_kmeans_ch <- function(train, test,  ...){
  model <- full(train, lambda = 1, join_unobserved = TRUE, order = ordering_ch(train))
  model <- stages_kmeans(model, k = 2, algorithm = "Hartigan-Wong")
  predict(model, newdata = test, class = "answer")
}

st_naive_order <- function(train, test,  order){
  model <- stages_hclust(full(train, join_unobserved = TRUE,
                              lambda = 1, order = order),
                         distance = "totvar", method = "complete")
  predict(model, newdata = test, class = "answer")
}

st_naive_hl_order <- function(train, test,  order){
  model <- stages_hclust(full(train, join_unobserved = TRUE,
                              lambda = 1, order = order),
                         distance = "hellinger", method = "complete")
  predict(model, newdata = test, class = "answer")
}

st_fbhc_order <- function(train, test,  order, ...){
  model <- stages_fbhc(full(train, join_unobserved = TRUE, lambda = 1, 
                            order = order))
  predict(model, newdata = test, class = "answer")
}

st_bj_kl_order <- function(train, test, order, ...){
  model <- stages_bj(full(train, join_unobserved = TRUE, lambda = 1, 
                          order = order), distance = "kullback", thr = 0.2)
  predict(model, newdata = test, class = "answer")
}

st_bhc_tan_cl <- function(train, test,  order){
  bnc <- bnclassify::lp(bnclassify::tan_cl("answer", train), train, smooth = 1)
  model <- as_sevt(bnlearn::as.bn.fit(bnclassify::as_grain(bnc)))
  model <- stages_bhc(sevt_fit(model, train, lambda = 1))
  predict(model, newdata = test, class = "answer")
}

st_bhc_tan_hc <- function(train, test){
  bnc <- bnclassify::lp(bnclassify::tan_hc("answer", train, k = 5), train, smooth = 1)
  model <- as_sevt(bnlearn::as.bn.fit(bnclassify::as_grain(bnc)))
  model <- stages_bhc(sevt_fit(model, train, lambda = 1))
  predict(model, newdata = test, class = "answer")
}

st_bhc_5db <- function(train, test){
  kdbk <- min(5, ncol(train)-1)
  bnc <- bnclassify::lp(bnclassify::kdb("answer", train, k = 5, kdbk = kdbk), train, smooth = 1)
  model <- as_sevt(bnlearn::as.bn.fit(bnclassify::as_grain(bnc)))
  model <- stages_bhc(sevt_fit(model, train, lambda = 1))
  predict(model, newdata = test, class = "answer")
}

st_bhc_3db <- function(train, test){
  kdbk <- min(3, ncol(train)-1)
  bnc <- bnclassify::lp(bnclassify::kdb("answer", train, k = 5, kdbk = kdbk), train, smooth = 1)
  model <- as_sevt(bnlearn::as.bn.fit(bnclassify::as_grain(bnc)))
  model <- stages_bhc(sevt_fit(model, train, lambda = 1))
  predict(model, newdata = test, class = "answer")
}

st_3db_cmi <- function(train, test){
  order <- ordering_cmi(train)
  allow <- cbind(order[1],order[-1])
  block <- bnlearn::ordering2blacklist(order)
  dag <- bnlearn::tabu(train, whitelist = allow, blacklist = block, 
                       maxp = 3 + 1)
  bn <- bnlearn::bn.fit(dag, train)
  tree <- as_sevt(bn, order = order)
  model <- stages_bhc(sevt_fit(tree, train, lambda = 1))
  predict(model, test, class = "answer")
}
