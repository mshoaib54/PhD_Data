args <- commandArgs(trailingOnly = TRUE)
datasets <- "monks1"
if (length(args) > 0){
    datasets <- args
}
library(ggplot2)
library(data.table)
stats <- c("accuracy", "balanced_accuracy", "auc")

alldata <- data.table()
for (dataset in datasets){
    TABLE <- readRDS(paste0("TABLE_ORDER_", dataset, ".rds"))
    ## compute averages
    AVG <- apply(TABLE, c(1,2,3), mean, na.rm = TRUE)
    ## transform to data.table
    data <- as.data.table(AVG)
    data$data <- dataset
    alldata <- rbind(alldata, data)
}

classifiers_cmi <- sub("order", "cmi" ,unique(alldata$classifier))
classifiers <- c(classifiers_cmi)
AVG_OR <- readRDS("AVG.rds")


for (st in stats){
data_cmi <- as.data.table(AVG_OR)[data %in% datasets & 
				    stat == st & 
				    classifier %in% classifiers]

ggplot(data = alldata[stat %in% st], 
       aes(y = sub("_order","",classifier), x = value, 
	   group = sub("_order","",classifier))) + 
  facet_grid(cols = vars(data), scales = "free_x") + 
  geom_violin() + 
#  geom_boxplot() + 
  geom_point(aes(y = sub("_cmi", "", classifier), 
                 x = value),
             size = 3, 
             shape = 4,
             data = data_cmi, show.legend =  TRUE) + 
  theme_bw() +  
  #theme(axis.ticks.y = element_blank(), 
  #	axis.text.y = element_blank()) + 
  ylab("") + 
  scale_y_discrete(labels = c("ST_FBHC", "ST_Naive"))+
  xlab(st) + 
  labs(color = "ordering", shape = "ordering") + 
  ggsave(paste0("plot_", st, "_order",".pdf"), 
	 width = 4.5, height = 2.5, units = "in")
}
