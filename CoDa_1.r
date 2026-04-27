  
######  Supporting file for the paper 
######  Dually Flat Affine Geometry of the Probability Simplex Bases on Kaniadakis’ Logarithm for Applications to Compositional Data Statistics
######  by Giovanni Pistone and Muhammad Shoaib

######  CHANGE AS NECESSARY
setwd("/Users/muhammadshoaib/Desktop/IG notes")

## ## ## ## Load and prepare the Dataset
#install.packages("Hmisc")
#install.packages("readxl")
library("Hmisc")
library("readxl")
library(ggplot2)
library(tidyverse)
library(reshape2)

data <- read_excel("CoDa.xlsx")

# Convertng the "," to "." 
data[,-c(1,2,3,4)] <- sapply(data[,-c(1,2,3,4)], function(x) as.numeric(gsub(",", ".",x)))
CoDa <- data[,-c(1,2,3,4)]

# Computing the Proportions
rowsum <- rowSums(CoDa)
com_data <- CoDa/rowsum 

## ## ## ## define useful funtions
# the Nudit's function
A<- function(p) {
  2*(p^2)/(1+p^2)
}
# Escort function
ptilde <- function(p){
  A(p)/sum(A(p))
}
# Karniadakis expectation of u wrt p
Ek <- function(p,u){
   sum(u*ptilde(p))
}
# Karniadakis logarithm
Lk <- function(p){
  1/2 * (p-1/p)
}
# Karniadakis divergence between two compositional data points
Kdiv <- function(p,q) {
  Ek( p , Lk(p) ) - Ek( p , Lk(q) )
}

## ## ## ## compute the Karniadakis divergence between any row of Ecom_data (p) and any other row of Ecom_data (q)
nrow <-  nrow(com_data)
Kdiv_com_data <- matrix(ncol=nrow, nrow = nrow) # initialise the divergence matrix 
for (p in 1:nrow)  for (q in 1:nrow)  { 
  Kdiv_com_data[p,q] <- Kdiv( com_data[p,] , com_data[q,] ) 
  }

# plot the heat map of Karniadakis divergence between rows of Ecom_data
rowNames <- data$Year  # Names of compositional data points
colnames(Kdiv_com_data) <- rowNames  
rownames(Kdiv_com_data) <- rowNames 
heatmap(Kdiv_com_data,  Colv = NA, Rowv = NULL )
heatmap(round( Kdiv_com_data, 1)  ,  symm = TRUE )

######
#order row/columns of the heat map
melted_kan <- melt(Kdiv_com_data)
p<-ggplot(data = melted_kan, aes(x=Var1, y=Var2,
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
p + labs(x = "", y = "")

## ## ## ## compute the Karniadakis divergence between uniform distribution (p) and rows of Ecom_data (q)
p <- rep(1,ncol(com_data))/ncol(com_data) 
Kdiv_com_Unif <- c()
for (q in 1:nrow)   {
  Kdiv_com_Unif[q] <- Kdiv( p , com_data[q,] )
  }
Kdiv_com_Unif
#####
## ## ## ##



#check analitically whether logk(p/q) = or not to logk(p)-logk(q)



## ## ## ## notes for future 
library(robCompositions)

