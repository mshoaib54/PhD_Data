#Libraries

library(robCompositions)
library("Hmisc")
library("readxl")
library(ggplot2)
library(tidyverse)
library(reshape2)
# Load the Dataset
data <- read_excel("/Users/muhammadshoaib/Desktop/IG notes/CoDa.xlsx")
# Convertng the "," to "." 
data[,-c(1,2,3,4)] <- sapply(data[,-c(1,2,3,4)], function(x) as.numeric(gsub(",", ".",x)))
head(data)
is.na(data)
CoDa <- data
# Computing the Proportions

rowsum <- rowSums(CoDa[,-c(1,2,3,4)])
CoDa <- cbind(CoDa, rowsum)
CoDa <- sweep(CoDa[,-c(1,2,3,4,15)], 1, rowSums(CoDa[,-c(1,2,3,4,15)]), FUN="/")
CoDa <- cbind(data[,c(1,2,3,4)],CoDa)
Total<- rowSums(CoDa[,-c(1,2,3,4)])
CoDa <- cbind(CoDa,Total)
CoDa_Rows_Com <- CoDa[-c(1,2,3,4,15)]


# Checking the E(u) = 0 
# Karniadakis expectation of u wrt p
Ek <- function(p,u){
  sum(u * ptilde(p))
}
## ## ## ## define useful funtions
# Karniadakis expectation of u wrt p
kan_div <- function(p,q){
  A <- 2*(p^2)/(1+p^2)
  Kl <- 1/2 * (p-1/p)
  ptilde <- A(p)/sum(A(p))
  k_div <- Ek( p , Kl ) - Ek( p , Kl )
  return(k_div)
}
## ## ## ## compute the Karniadakis divergence between any row of com_data (p) and any other row of Ecom_data (q)
nrow <-  nrow(CoDa_Rows_Com)
KanD_com_data <- matrix(ncol=nrow, nrow = nrow) # initialise the divergence matrix 
for (p in 1:nrow)  
  for (q in 1:nrow)  { 
    KanD_com_data[p,q] <- kan_div( CoDa_Rows_Com[p,] , CoDa_Rows_Com[q,] ) 
  }
# plot the heat map of Karniadakis divergence between rows of Ecom_data
rowNames <- data$Year  # Names of compositional data points
colnames(KanD_com_data) <- rowNames  
rownames(KanD_com_data) <- rowNames 
heatmap(KanD_com_data,  Colv = NA, Rowv = NULL )
heatmap(round( KanD_com_data, 1)  ,  symm = TRUE )

######
#order row/columns of the heat map
melted_kan <- melt(KanD_com_data)
p_1<-ggplot(data = melted_kan, aes(x=Var1, y=Var2,
                                 fill=value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient(low = "white", high = "blue",
                      name = "Kaniadakis\nDivergence")+
  theme_bw()+
  theme(axis.text.x = element_text(vjust = 1))  +
  theme(axis.title.y = element_text(vjust = 1))+
  theme(axis.text.x = element_text(angle = 0, vjust = 1))+
  scale_x_continuous(labels=as.character(rowNames),breaks=rowNames)+
  scale_y_continuous(labels=as.character(rowNames),breaks=rowNames)+
  #coord_fixed()+
  ggtitle("Kaniadakis Divergence Matrix heatmap") +
  theme(plot.title = element_text(hjust = 0.5))
p_1 + labs(x = "", y = "")
## ## ## ## compute the Karniadakis divergence between uniform distribution (p) and rows of Ecom_data (q)
p_uni <- rep(1,ncol(CoDa_Rows_Com))/ncol(CoDa_Rows_Com) 
Kdiv_com_Unif <- c()
for (q in 1:nrow)   {
  Kdiv_com_Unif[q] <- kan_div( p_uni , CoDa_Rows_Com[q,] )
}
Kdiv_com_Unif
#####
## ## ## ##

