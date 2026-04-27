library("reshape2")
library("stagedtrees")

evaluate_models <- function(path) {
  models <- readRDS(path)
  data_path <- models[[1]][[1]]$data_path
  gt <- readRDS(data_path)
  true_model <- gt$true_model
  results <- lapply(models, function(models_1){
    lapply(models_1, function(res){
      if (!class(res$model) == "sevt"){
        data.frame(
          BIC = NA,
          AIC = NA,
          hamming = NA,
          cid = NA,
          time = NA,
          sample_size = res$sample_size,
          divergence = res$divergence,
          k = res$k,
          alg = res$alg,
          method = res$method,
          p = gt$P,
          q = gt$Q,
          ktrue = gt$Ktrue
        )
      }else{
        data.frame(
          BIC = BIC(res$model),
          AIC = AIC(res$model),
          hamming = hamming_stages(true_model, res$model, FUN = mean),
          cid = NA,
          #cid = cid(true_model, res$model, FUN = mean)$cid,
          time = res$time,
          sample_size = res$sample_size,
          divergence = res$divergence,
          k = res$k,
          alg = res$alg,
          method = res$method,
          p = gt$P,
          q = gt$Q,
          ktrue = gt$Ktrue
        )
      }

    })
  })

  reshape2::melt(results, id.vars = c("sample_size",
                                      "alg", "method",
                                      "divergence", "k",
                                      "p", "q", "ktrue"))
}
