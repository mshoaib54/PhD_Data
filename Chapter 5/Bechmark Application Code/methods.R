### all classifers are a function of train and test (plus additional args)
### and return probabilities, predictions and cutoff.

source("methods/st_methods.R")
#source("methods/bn_methods.R")
#source("methods/bnc_methods.R")
#source("methods/nnet_methods.R")
#source("methods/glm_methods.R")
source("methods/rf_methods.R")
#source("methods/logistic_methods.R")
#source("methods/tree_methods.R")
#source("methods/discriminant_analysis_methods.R")
#source("methods/naive_methods.R")
#source("methods/naive_2_methods.R")
#source("methods/boosting_methods.R")
#source("methods/bagging_methods.R")
#source("methods/adaboost_methods.R")
#source("methods/svm_methods.R")
#source("methods/gam_methods.R")

classifiers <- c("st_naive_totvar_cmi", 
                 "st_naive_jefferys_cmi",
                 "st_naive_kaniadakis_cmi",
                 "st_naive_jensen_shannon_cmi",
                 "st_naive_hellinger_cmi")


