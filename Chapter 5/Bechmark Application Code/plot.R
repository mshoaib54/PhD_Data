library(ggplot2)
library(data.table)

datasets <- factor(read.table("binary_datasets_names_final.tsv")[, 1])

nreps <- 10

### select which methods to plot
classifiers <- c("st_naive_totvar_cmi", 
                 "st_naive_jefferys_cmi",
                 "st_naive_kaniadakis_cmi",
                 "st_naive_jensen_shannon_cmi",
                 "st_naive_hellinger_cmi")

good  <-  c("st_naive_totvar_cmi", 
            "st_naive_jefferys_cmi",
            "st_naive_kaniadakis_cmi",
            "st_naive_jensen_shannon_cmi",
            "st_naive_hellinger_cmi")


## read .rds files
AVG <- readRDS("AVG.rds")
AVG <- AVG[,,classifiers] 
dimnames(AVG)[3] <- list(classifier = good) 
data <- as.data.table(AVG)

## plot accuracies
PLT <- ggplot(data = data[stat %in% c("auc", "balanced_accuracy")], 
	      aes(y = data, x = value, shape = classifier, 
                        group = classifier, color = classifier)) + 
  geom_jitter(height = 0.2, width = 0, alpha = 0.5) + facet_grid(cols = vars(stat), scales = "free") + 
  scale_colour_discrete(name = "Classifier", 
		     labels = good) + 
  scale_shape_manual(name = "Classifier", 
              labels = good,
	      values = c(rep(15, 9), rep(19, 7))) + 
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 30),
    legend.box.spacing = unit(0.5, "lines"),
    legend.box.margin = ggplot2::margin(
      t = 0,
      r = 5,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    plot.margin = ggplot2::margin(
      t = 0,
      r = 5,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) + xlab("") 
ggsave("plot_accuracy.pdf", PLT,  width = 7, height = 6, units = "in")



## plot time
PLT <- ggplot(data = data[stat == "time"], aes(y = data, x = value, 
                                                                      group = classifier, color = classifier)) + 
  geom_jitter(height = 0.2, alpha = 0.5) + 
  theme_bw() +
  scale_x_log10()+
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 30),
    legend.box.spacing = unit(0.5, "lines"),
    legend.box.margin = ggplot2::margin(
      t = 0,
      r = 15,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    plot.margin = ggplot2::margin(
      t = 0,
      r = 15,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) + guides(color=guide_legend(nrow=3,byrow=TRUE)) +  
  xlab("seconds")

ggsave("plot_time.pdf", PLT, width = 5, height = 6, units = "in")


## plot spec - sens - fp - fn
PLT <- ggplot(data = data[stat %in% c("sens", "fn", "spec", "fp")], aes(y = data, x = value, 
                                                                                                          group = classifier, color = classifier)) + 
  geom_jitter(height = 0.2, width = 0, alpha = 0.5) + facet_grid(cols = vars(stat), scales = "free") + 
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 30),
    legend.box.spacing = unit(0.5, "lines"),
    legend.box.margin = ggplot2::margin(
      t = 0,
      r = 5,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    plot.margin = ggplot2::margin(
      t = 0,
      r = 5,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) + xlab("")

ggsave("plot_confusion_matrix.pdf", PLT, width = 7, height = 6, units = "in")




## plot auc - precision - f1
PLT <- ggplot(data = data[stat %in% c("accuracy", "precision", "f1")], aes(y = data, x = value, 
                                                                                               group = classifier, color = classifier)) + 
  geom_jitter(height = 0.2, width = 0, alpha = 0.5) + facet_grid(cols = vars(stat), scales = "free") + 
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 30),
    legend.box.spacing = unit(0.5, "lines"),
    legend.box.margin = ggplot2::margin(
      t = 0,
      r = 5,
      b = 0,
      l = 0,
      unit = "pt"
    ),
    plot.margin = ggplot2::margin(
      t = 0,
      r = 5,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) + xlab("")
 ggsave("plot_auc_precision_f1.pdf", PLT, width = 7, height = 6, units = "in")




library(xtable)

print(xtable(t(AVG['accuracy',, good]), label = "tab:acc", digits = 4,
	     caption = "Accuracy for all the considered classifiers over the datasets in the experiments."),
      floating = FALSE, 
      file = "accuracy_table.tex", booktabs = TRUE, scalebox = 0.5)
print(xtable(t(AVG['balanced_accuracy',, good]), label = "tab:bacc", digits = 4,
caption = "Balanced accuracy for all the considered classifiers over the datasets in the experiments."),
      floating = FALSE, 
      file = "balanced_accuracy_table.tex", booktabs = TRUE, scalebox = 0.5)
print(xtable(t(AVG['auc',, good]), label = "tab:auc", digits = 4,
caption = "Area under the ROC curve for all the considered classifiers over the datasets in the experiments."),
      floating = FALSE,
      file = "auc_table.tex", booktabs = TRUE, scalebox = 0.5)







plot_auc_precision_f1_like_accuracy <- function(
    data,
    legend_position = "bottom",
    free_y = TRUE
) {
  library(data.table)
  library(ggplot2)
  
  DT <- as.data.table(data)
  
  # --- checks
  req <- c("stat", "data", "classifier", "value")
  miss <- setdiff(req, names(DT))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))
  
  # --- stats
  DT <- DT[stat %in% c("accuracy", "f1")]
  if (nrow(DT) == 0) stop("No rows for stat in {auc, precision, f1}. Check unique(DT$stat).")
  
  # --- find Jensen-Shannon classifier key robustly
  js_name <- grep("jensen|shannon|js", unique(DT$classifier), ignore.case = TRUE, value = TRUE)
  js_name <- js_name[1]
  if (is.na(js_name) || length(js_name) == 0) {
    stop("Couldn't find Jensen-Shannon in DT$classifier. Run: unique(DT$classifier)")
  }
  
  # --- keep only requested classifiers
  keep <- c(
    "st_naive_jensen_shannon_cmi",
    "st_naive_jefferys_cmi",
    "st_naive_kaniadakis_cmi",
    "st_naive_hellinger_cmi",
    "st_naive_totvar_cmi",
    js_name
  )
  
  keep <- intersect(keep, unique(DT$classifier))
  if (length(keep) < 2) stop("Too few classifiers found. Check DT$classifier names.")
  
  labels <- c("Jefferys", "Kaniadakis", "Hellinger", "totvar", "JensenShannon")[seq_along(keep)]
  
  DT <- DT[classifier %in% keep]
  DT[, classifier := factor(classifier, levels = keep, labels = labels)]
  
  # stat order (optional)
  DT[, stat := fifelse(stat == "auc", "accuracy", stat)]
  DT[, stat := factor(stat, levels = c("accuracy", "f1"))]
  
  #DT[, stat := factor(stat, levels = c("auc", "precision", "f1"))]
  
  # --- plot (BOLD dots like plot_accuracy.pdf)
  PLT <- ggplot(
    DT,
    aes(
      x = classifier,
      y = value,
      fill = classifier
    )
  ) +
    geom_point(
      shape = 21,        # filled circle with outline
      size  = 2.8,       # make it bolder
      alpha = 1,         # fully opaque (darker)
      color = "black",   # outline color
      stroke = 0.35      # outline thickness
    ) +
    facet_grid(
      rows = vars(stat),
      cols = vars(data),
      scales = if (free_y) "free_y" else "fixed",
      switch = "y"
    ) +
    scale_fill_discrete(labels = labels) +
    theme_bw() +
    theme(
      legend.position = legend_position,
      legend.title = element_blank(),
      
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      
      # Hide x labels (legend carries names, like your accuracy PDF)
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      
      strip.placement = "outside",
      # ✅ bold dataset names + ✅ bold stat names
      strip.text.x = element_text(face = "bold"),
      strip.text.y = element_text(face = "bold"),
      
      legend.box.spacing = unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t=0, r=5, b=0, l=0, unit="pt"),
      plot.margin = ggplot2::margin(t=0, r=5, b=0, l=0, unit="pt")
    )
  
  return(PLT)
}
PLT <- plot_auc_precision_f1_like_accuracy(data)
PLT
ggsave("plot_auc_precision_f1.pdf", PLT, width = 14, height = 4)
