##################
##### PACKAGES ###
##################
library(ggplot2)
library(scales)
library(grid)
library(gridExtra)
library(reshape2)

#########################
#### LOADING DATA #######
#########################

# Change to directory where bottleneck data is stored (comment out depending on if you're using OSX or Ubuntu)
setwd("~/Documents/OneDrive - University of Exeter/Data/Bottlenecks/Phage_bn_exp/") #Work OSX

# Need the 'count_t3_csv.csv' file
t3 = read.csv(file.choose(), header=T)
t3 = na.exclude(t3)
t3$bottleneck = as.factor(t3$bottleneck)
t3$ID = as.factor(t3$ID)

t3 = melt(t3, measure.vars = c('pfu', 'cfu'))

quartz()
ggplot(t3, aes(y=value, x=bottleneck))+
  geom_boxplot(aes(colour=variable), position=pd)+
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  labs(y=expression(bold("p.f.u/c.f.u. ml"*{}^{-1}*"")), x=expression(bold('Bottleneck size')))+
  theme_bw()+
  scale_colour_discrete(name = 'Measurement',
                         breaks = c('pfu', 'cfu'),
                         labels = c('Phage titer', 'Bacterial density'))+
  theme(legend.title = element_text(face='bold'))
