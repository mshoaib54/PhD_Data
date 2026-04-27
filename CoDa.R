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

# Kaniadakis logarithm 
kan_log <- function(x){
  1/2*(x -(1/x))
}
## define the Nudit's function
A<- function(x) {
  (2*(x^2)/(1+x^2))
}

# Apply Kaniadakis Logarithm on the original data
for(i in 1:nrow(CoDa_Rows_Com)) {       # for-loop over rows
  CoDa_Rows_Com[i, ] <- lapply(CoDa_Rows_Com[i,],kan_log)
}
# Apply escort probability on transformed data(Kaniadakis Log)
a.p<-ap<-lapply(CoDa_Rows_Com[1,],A)
a.p<-unlist(a.p)
a.p<- a.p/sum(a.p)
# Kaniadakis Divergence on transformed data
div_1<-CoDa_Rows_Com - c(CoDa_Rows_Com[1, ])
#div_1<-div_1[-1,]
#div_1<-rbind(CoDa_Rows_Com[1,],div_1)
D_1 <- mapply('*',div_1,a.p)
d_1<-rowSums(D_1)
######
# Taking Uniform Distribution as reference for the transformed data
#uni_dist1<-runif(10, min = 1/10, max = 1/10)
#uni_div1 <- CoDa_Rows_Com - c(CoDa_Rows_Com)
#uni_D1 <- mapply('*',uni_div1,a.p )


#####
#Kaniadakis Divergence on original Compositinal data 

com_data <- CoDa[-c(1,2,3,4,15)]
# Finding escort probability
ap<-lapply(com_data[1,],A)
ap<-unlist(ap)
ap<- ap/sum(ap)
# kaniadakis divergence
div<-com_data - c(com_data[1, ])
#div<-div[-1,]
#div<-rbind(com_data[1,],div)
# Divergence matrix
D <- mapply('*',div,ap)
d_2 <- rowSums(D)

######
# Taking Uniform Distribution as reference
uni_dist<-runif(10, min = 1/10, max = 1/10)
uni_div <- com_data - c(uni_dist)
uni_D <- mapply('*',uni_div,ap )
uni_row_sum <- rowSums(uni_D)

#####
#s(q) = u
s_p_q <- -exp(div_1 - d_2)

# Checking the E(u) = 0 

## ## ## ## define useful funtions
# Karniadakis expectation of u wrt p
kan_div <- function(p,q){
  A <- 2*(p^2)/(1+p^2)
  L_p <- 1/2 * (p-1/p)
  L_q <- 1/2 * (q-1/q)
  ptilde <- A*L_p/sum(A*L_p)
  k_div <- sum(p*ptilde*L_p) - sum(p*L_q*ptilde)
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
Nrow <-  nrow(CoDa_Rows_Com)
Kdiv_Uni_data <- matrix(ncol=Nrow, nrow = Nrow) # initialise the divergence matrix 
for (p in 1:nrow)  for (q in 1:nrow)  { 
  Kdiv_Uni_data[p,q] <- kan_div( p_uni , CoDa_Rows_Com[q,] ) 
}
for( p in 1:Norw){
uni_div[p,q] <- Kdiv_Uni_data[q,] - c(p_uni)
}
