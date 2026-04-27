source("evaluate_models.R")
library("pbapply")

### Evaluation ####
all_fitted <- list.files("fitted_model_1/simulated_data/", pattern = "*.rds",
                         recursive = TRUE,
                         full.names = TRUE)

DD <- reshape2::melt(pblapply(all_fitted, evaluate_models),
                       measure.vars = "value", cl = 5)

saveRDS(DD, file = "DD_all.rds")
