  
######  Supporting file for the paper 
######  Dually Flat Affine Geometry of the Probability Simplex Bases on Kaniadakis’ Logarithm for Applications to Compositional Data Statistics
######  by Giovanni Pistone and Muhammad Shoaib

######  CHANGE AS NECESSARY
setwd("/Users/muhammadshoaib/Desktop/IG notes")

## ## ## ## Load and prepare the Dataset

library("Hmisc")
library("readxl")
library(dplyr)
## for the graphics 
library(ggplot2)
library(tidyverse)
library(reshape2)
# For Trenary Diagrams
#devtools::install_git('https://bitbucket.org/nicholasehamilton/ggtern')
library(ggtern)
library(Ternary)
library(ggalt)
library(plotly)
library(ggpubr)
library(ggthemes)
data <- read_excel("CoDa.xlsx")

# Convertng the "," to "." 
data[,-c(1,2,3,4)] <- sapply(data[,-c(1,2,3,4)], function(x) as.numeric(gsub(",", ".",x)))
CoDa <- data[,-c(1,2,3,4)]

# Computing the Proportions
rowsum <- rowSums(CoDa)
com_data <- CoDa/rowsum 

## ## ## ## define useful functions
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
   sum(u * ptilde(p))
}
# Karniadakis logarithm
Kl <- function(p){
  1/2 * (p-1/p)
}
# Karniadakis divergence between two compositional data points
Kdiv <- function(p,q) {
  Ek( p , Kl(p) ) - Ek( p , Kl(q) )
}
# 

s_q <- function(p,q){
  (Kl(q) - Kl(p)) + Kdiv(p,q)
}
# Defining the transformed q
tran_q <- function(p,q){
  exp(Kl(p) - Kdiv(p,q)+ s_q(p,q))
}
# exponential displacement 
exponential_disp <- function(p,q) {
  -(Kl(p)-Kl(q) + Kdiv(p,q))
}
# mixture displacement 
mixture_disp <- function(p,q) {
 (q-p)/A(p)
}
## ## ## ## end of definition of useful functions

## ## ## ## compute the Karniadakis divergence between any row of Ecom_data (p) and any other row of Ecom_data (q)
nrow <-  nrow(com_data)
Kdiv_com_data <- matrix(ncol=nrow, nrow = nrow) # initialise the divergence matrix 
for (p in 1:nrow)  for (q in 1:nrow)  { 
  Kdiv_com_data[p,q] <- Kdiv(com_data[p,],com_data[q,] )
}


# plot the heat map of Karniadakis divergence between rows of Ecom_data
rowNames <- data$Year  # Names of compositional data points
colnames(Kdiv_com_data) <- rowNames  
rownames(Kdiv_com_data) <- rowNames 
#heatmap(Kdiv_com_data,  Colv = NA, Rowv = NULL )
#heatmap(round( Kdiv_com_data, 1)  ,  symm = TRUE )
#order row/columns of the heat map
melted_kan <- melt(Kdiv_com_data)
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
  theme(axis.text = element_text(face="bold"))+
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  scale_x_continuous(labels=as.character(rowNames),breaks=rowNames)+
  scale_y_continuous(labels=as.character(rowNames),breaks=rowNames)+
  #coord_fixed()+
  ggtitle("Heatmap of Kaniadakis Divergence Matrix") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(family = "Times New Roman"))
p_1 + labs(x = "", y = "")

## ## ## ## compute the Karniadakis divergence between uniform distribution (p) and rows of Ecom_data (q)
p <- rep(1,ncol(com_data))/ncol(com_data) 
Kdiv_com_Unif <- c()
for (q in 1:nrow)    Kdiv_com_Unif[q] <- Kdiv( p , com_data[q,] )
Kdiv_com_Unif

## ## ##
# Mixture Displacement on the data
mix_disp<-mapply(mixture_disp,com_data[1,], com_data[1:14,])
rowNames <- data$Year  # Names of compositional data points
colnames(mix_disp) <- c("Belgium", "Denmark" ,"France", "Germany","Greece",     
                        "Italy", "Netherland", "Spain","Switzerland","UK")  
rownames(mix_disp) <- rowNames 
#order row/columns of the heat map
melt_mix_dis <- melt(mix_disp)
p_2<-ggplot(data = melt_mix_dis, aes(x=Var2, y=Var1,
                                 fill=value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient(low = "white", high = "blue",
                      name = "Legend")+
  theme_bw()+
  theme(axis.text.x = element_text(vjust = 1))  +
  theme(axis.title.y = element_text(vjust = 1))+
  theme(axis.text.x = element_text(angle = 0, vjust = 1))+
  theme(axis.text = element_text(face="bold"))+
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  #scale_x_continuous(labels=as.character(rowNames),breaks=rowNames)+
  scale_y_continuous(labels=as.factor(rowNames),breaks=rowNames)+
  #coord_fixed()+
  ggtitle("Mixture Displacement") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(family = "Times New Roman"))
p_2  + labs(x = "", y = "")
## ## ## 
# Exponential Dispalcement
Exp_disp<-mapply(exponential_disp,com_data[1,], com_data[1:14,])

my_list <- list()
for( q in 1:nrow(com_data-1)){

    my_list[[q]] = exponential_disp(com_data[1,] ,com_data[q,])
   
}
Exp_data_m <- matrix(unlist(my_list), ncol=10, byrow=TRUE)
# Exponential Displacement on the data
rowNames <- data$Year  # Names of compositional data points
colnames(Exp_data_m) <- c("Belgium", "Denmark" ,"France", "Germany","Greece",     
                        "Italy", "Netherland", "Spain","Switzerland","UK")  
rownames(Exp_data_m) <- rowNames
#order row/columns of the heat map
melt_exp_dis <- melt(Exp_data_m)
p_3<-ggplot(data = melt_exp_dis, aes(x=Var2, y=Var1,
                                   fill=value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient(low = "white", high = "blue",
                      name = "Legend")+
  theme_bw()+
  theme(axis.text.x = element_text(vjust = 1))  +
  theme(axis.title.y = element_text(vjust = 1))+
  theme(axis.text.x = element_text(angle = 0, vjust = 1))+
  
  theme(axis.text = element_text(face="bold"))+
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  #scale_x_continuous(labels=as.character(rowNames),breaks=rowNames)+
  scale_y_continuous(labels=as.factor(rowNames),breaks=rowNames)+
  #coord_fixed()+
  ggtitle("Exponential Displacement") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(family = "Times New Roman"))
p_3  + labs(x = "", y = "")

#####
## Computing the Displacement from mean value
mean<- colMeans(com_data)
ex_di <- list()
for( q in 1:nrow(com_data-1)){
  
  ex_di[[q]] = exponential_disp(mean ,com_data[q,])
  
}
exp_displ <- matrix(unlist(ex_di), ncol=10, byrow=TRUE)
# Exponential Displacement from mean on the data
rowNames <- data$Year  # Names of compositional data points
colnames(exp_displ) <- c("Belgium", "Denmark" ,"France", "Germany","Greece",     
                          "Italy", "Netherland", "Spain","Switzerland","UK")  
rownames(exp_displ) <- rowNames
#order row/columns of the heat map
melt_exp_disp <- melt(exp_displ)
p_4<-ggplot(data = melt_exp_disp, aes(x=Var2, y=Var1,
                                   fill=value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient(low = "white", high = "blue",
                      name = "Legend")+
  theme_bw()+
  theme(axis.text.x = element_text(vjust = 1))  +
  theme(axis.title.y = element_text(vjust = 1))+
  theme(axis.text.x = element_text(angle = 0, vjust = 1))+
  
  theme(axis.text = element_text(face="bold"))+
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  #scale_x_continuous(labels=as.character(rowNames),breaks=rowNames)+
  scale_y_continuous(labels=as.factor(rowNames),breaks=rowNames)+
  #coord_fixed()+
  ggtitle("Exponential Displacement from Mean") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(family = "Times New Roman"))
p_4  + labs(x = "Countries", y = "Years")

#####
# Mixture Displacement from mean
m_dis <- mapply(mixture_disp,mean,com_data[1:14,])


rowNames <- data$Year  # Names of compositional data points
colnames(m_dis) <- c("Belgium", "Denmark" ,"France", "Germany","Greece",     
                         "Italy", "Netherland", "Spain","Switzerland","UK")  
rownames(m_dis) <- rowNames
#order row/columns of the heat map
melt_mix_disp <- melt(m_dis)
p_5<-ggplot(data = melt_mix_disp, aes(x=Var2, y=Var1,
                                    fill=value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient(low = "white", high = "blue",
                      name = "Legend")+
  theme_bw()+
  theme(axis.text.x = element_text(vjust = 1))  +
  theme(axis.title.y = element_text(vjust = 1))+
  theme(axis.text.x = element_text(angle = 0, vjust = 1))+
  
  theme(axis.text = element_text(face="bold"))+
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  #scale_x_continuous(labels=as.character(rowNames),breaks=rowNames)+
  scale_y_continuous(labels=as.factor(rowNames),breaks=rowNames)+
  #coord_fixed()+
  ggtitle("Mixture Displacement from Mean") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(family = "Times New Roman"))
p_5  + labs(x = "Countries", y = "Years")
  

ggarrange(p_2,p_5 +
            rremove("ylab"),
          ncol = 2, nrow = 1, labels  = "AUTO", hjust = c(-5, -2.5))
ggarrange(p_3,p_4 +
            rremove("ylab"),
          ncol = 2, nrow = 1, labels  = "AUTO", hjust = c(-5, -2.5))
          # labels = c("A", "B", "C","D"))























#######
# Plotting Ternary Graphs
new_com_data <- cbind(com_data,data[c(1,3)])
colnames(new_com_data)[12] <- "PC"
ggtern(data=new_com_data, aes(x=UK,y=Switzerland, z=Denmark)) +
  #geom_point() +
 # labs(title="Population structure, 2015") +
 # theme_rgbw()
geom_point(aes(fill = PC),
           size = 6,
           shape = 21,
           color = "black") +
 ggtitle("Trenary Diagram") +
 labs(fill = "Period Code") +
 theme_rgbw() +
 theme(legend.position = c(0,1),
 legend.justification = c(1, 1))

ggtern(data=new_com_data, aes(x=UK,y=Switzerland, z=Denmark
                                     ,fill=PC,shape=PC)) + 
  theme_bw() + 
  theme_legend_position('tr') + 
  geom_encircle(alpha=0.5,size=1) + 
  geom_point() +
  labs(title    = "Trenary Diagram",
       subtitle = "")
data(Fragments)
arrangement = list()
for(base in c('identity','ilr')){
  x = ggtern(Fragments,aes(Qm,Qp,M)) +
    stat_density_tern(geom='polygon',
                      aes(fill=..level..),
                      base=base,  ###NB Base Specification
                      colour='grey50') + 
    theme_dark() + geom_point() + 
    ggtitle(sprintf("Basis: %s",base)) +
    scale_fill_gradient(low='green',high='red') +
    limit_tern(.5,1,.5)
  arrangement[[length(arrangement) + 1]] = x
}
grid.arrange(grobs = arrangement,nrow=1)

# we create the ternary plot using plotly
p3 <- plot_ly(
  new_com_data, a = ~'% UK', b = ~'Switzerland', c = ~'Denmark',
  frame=~Year,
  color = ~`Planning Region`, type = "scatterternary", colors=~colors,
  size = ~Total,
  text = ~paste('Switzerland',sep='', round('Switzerland',1),'%',
                '<br>UK',
                round('UK',1),'%', '<br>Old:',
                round('Denmark',1),'%','<br>Subzone:', Subzone, hoverinfo="text",
                '<br>Planning Area:', `Planning Area`),
  marker = list(symbol = 'circle', opacity=0.4, sizemode="diameter", sizeref=2,
                line = list(width = 2, color = '#FFFFFF'))) %>%
  layout(
    title="Singapore Population 2002-2017",
    ternary=list(aaxis=list(title="UK",min=0.5),
                 baxis = list(title="Switzerland",min=0.2),
                 caxis = list(title="Denmark")),
    paper_bgcolor = 'rgb(243, 243, 243)',
    plot_bgcolor = 'rgb(243, 243, 243)'
  ) %>%
  animation_slider(
    currentvalue = list(prefix = "YEAR ", font = list(color="red"))
  ) %>%
  animation_opts(
    2000, redraw = FALSE
  )
p3

## ## ## ## notes for future 
library(robCompositions)
