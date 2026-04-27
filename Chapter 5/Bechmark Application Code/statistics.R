library(pROC)

statistics <- c("time", "accuracy", "balanced_accuracy", "f1", "sens", "spec", 
                "fn", "fp", "error", "precision")


time <- function(res, true){
  res$time
}

accuracy <- function(res, true){
   sum(diag(table(res$predict, true))) / length(true)
}

balanced_accuracy <- function(res, true){
   P <- table(true)[1]
   N <- length(true) - P
   TT <- table(res$predict, true) 
   TPR <- ifelse(P > 0,TT[1,1] / P, 1)
   TNR <- ifelse(N > 0, sum(diag(TT)[-1]) / N, 1)  
   (TPR + TNR) / 2
} 

f1 <- function(res, true){
  TT <- table(res$predict, true) 
  TP <- TT[1,1] 
  FP <- sum(TT[1,]) - TP
  FN <- sum(TT) - sum(TT[1,]) - sum(diag(TT)[-1])
  (2*TP) / (2*TP + FP + FN)
}

sens <- function(res, true) {
   P <- table(true)[1]
   TT <- table(res$predict, true) 
   as.numeric(ifelse(P > 0,TT[1, 1] / P, 1))
}

spec <- function(res, true) {
   P <- table(true)[1]
   N <- length(true) - P
   TT <- table(res$predict, true) 
   as.numeric(ifelse(N > 0, sum(diag(TT)[-1]) / N, 1))
}

fn <- function(res, true) {
   P <- table(true)[1]
   TT <- table(res$predict, true) 
   as.numeric(ifelse(P > 0,TT[2, 1] / P, 1))
}

fp <- function(res, true) {
   P <- table(true)[1]
   N <- length(true) - P
   TT <- table(res$predict, true) 
   as.numeric(ifelse(N > 0, TT[1, 2] / N, 1))
}

error <- function(res, true) {
   1 - (sum(diag(table(res$predict, true))) / length(true))
}

precision <- function(res, true) {
   TT <- table(res$predict, true) 
   P <- rowSums(TT)[1]
   as.numeric(ifelse(P > 0, TT[1, 1] / P, 1))
}



