##################
##### PACKAGES ###
##################
library(ggplot2)
library(scales)
library(grid)
library(gridExtra)
library(reshape2)
library(ggpubr)

###################
### FUNCTIONS #####
###################

# Make printing out some stats for summary tables/figures quicker and easier
model_stats = function(model){
  sum = coef(model)
  conf = confint(model, level=c(0.95))
  print(c(sum[1]))
  print(c(conf[1,1]))
  print(c(conf[1,2]))
  
  stats = data.frame(sum[1], conf[1,1], conf[1,2])
  clip = pipe('pbcopy', 'w')
  write.table(stats, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  print('Coefficients copied to the clipboard')
  
}

#### Compare AIC values
compare_AICs = function(df){
  print(df)
  col_len = length(df[,2])
  AIC_min = abs(min(df[,2]))
  for (i in seq(1, col_len, 1)){
    print( (abs(df[i,2])) - AIC_min)
  }
}


#### GRAPHING PARAMETERS
OG_bottleneck = c('100', '1000', '10000', '100000', '1000000','10000000', '100000000', '1000000000')

bottleneck_names_facet = list(
  '100'      = expression('10'^2*''),
  '1000'     = expression('10'^3*''),
  '10000'    = expression('10'^4*''),
  '1e+05'   =  expression('10'^5*''),
  '1e+06'  = expression('10'^6*''),
  '1e+07' = expression('10'^7*''),
  '1e+08' = expression('10'^8*''),
  '1e+09' = expression('10'^9*'')
)

bottleneck_names_legend = c(
  expression('10'^2*''), 
  expression('10'^3*''), 
  expression('10'^4*''),
  expression('10'^5*''),
  expression('10'^6*''),
  expression('10'^7*''),
  expression('10'^8*''),
  expression('10'^9*'')
)


#########################
#### LOADING DATA #######
#########################

# Change to directory where bottleneck data is stored (comment out depending on if you're using OSX or Ubuntu)
#setwd("~/Documents/OneDrive - University of Exeter/Data/Bottlenecks/Bacteria_and_phage_bn_exp/") #Work OSX
#setwd("~/OneDrive/Data/Bottlenecks/") #Ubuntu
setwd("~/OneDrive - University of Exeter/Data/Bottlenecks/Phage_bn_exp/") #Home OSX

crispr = read.csv('./data/all_data_t3.csv')
crispr$bottleneck = as.factor(crispr$bottleneck)
crispr$mean.CR1.spacers = as.numeric(crispr$mean.CR1.spacers)
crispr$mean.CR2.spacers = as.numeric(crispr$mean.CR2.spacers)
crispr$mean.total.spacers = as.numeric(crispr$mean.total.spacers)

str(crispr)

##### Summary figures to eyeball the data

# Boxplot to get an impression of the differences between bottleneck treatments
CR2_all_plot_phage_bn = ggplot(crispr, aes(y=mean.CR2.spacers, x=bottleneck))+
  geom_boxplot()+
  labs(x='Bottleneck size', y="Spacer number")+
  ggtitle('CRISPR2 spacers in all clones')+
  theme_bw()+
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 20))+
  theme(axis.title = element_text(face='bold', size=16))+
  
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text = element_text(size=12))+
  
  scale_y_continuous(breaks=c(seq(0,2.5,0.5)))

CR2_all_plot_phage_bn

CR1_all_plot_phage_bn = ggplot(crispr, aes(y=mean.CR1.spacers, x=bottleneck))+
  geom_boxplot()+
  labs(x='Bottleneck size', y="Spacer number")+
  ggtitle('CRISPR1 spacers in all clones')+
  theme_bw()+
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 20))+
  theme(axis.title = element_text(face='bold', size=16))+
  
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text = element_text(size=12))+
  
  coord_cartesian(ylim=c(0.000,0.375))
  
  
CR1_all_plot_phage_bn

total_all_plot_phage_bn = ggplot(crispr, aes(y=mean.total.spacers, x=bottleneck))+
  geom_boxplot()+
  geom_boxplot()+
  labs(x='Bottleneck size', y="Spacer number")+
  ggtitle('Total spacers in all clones')+
  theme_bw()+
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 20))+
  theme(axis.title = element_text(face='bold', size=16))+
  
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text = element_text(size=12))+
  
  scale_y_continuous(breaks=c(seq(0,2.5,0.5)))
total_all_plot_phage_bn

ggsave('CR1_all_plot.png', plot=CR1_all_plot_phage_bn, device='png',
       path='./figs/', width=7, height = 5, units=c('in'), dpi=600)  
ggsave('CR2_all_plot.png', plot=CR2_all_plot_phage_bn, device='png',
       path='./figs/', width=7, height = 5, units=c('in'), dpi=600)  
ggsave('total_all_plot.png', plot=total_all_plot_phage_bn, device='png',
       path='./figs/', width=7, height = 5, units=c('in'), dpi=600)  
##############
### MODELS ###
##############

#### CR1 SPACERS
m1 = glm(mean.CR1.spacers~bottleneck, data=crispr)
m2 = glm(mean.CR1.spacers~pfu, data=crispr)
m3 = glm(mean.CR1.spacers~bottleneck*pfu, data=crispr)

anova(m3, test='LRT')
AICs = AIC(m1, m2, m3)
compare_AICs(AICs)

par(mfrow=c(2,2))
plot(m1)
plot(m2)
plot(m3)

crispr$bottleneck = relevel(crispr$bottleneck, ref='10000000')
crispr$bottleneck = relevel(crispr$bottleneck, ref='1000000')
crispr$bottleneck = relevel(crispr$bottleneck, ref='100000')
crispr$bottleneck = relevel(crispr$bottleneck, ref='10000')
crispr$bottleneck = relevel(crispr$bottleneck, ref='1000')
crispr$bottleneck = relevel(crispr$bottleneck, ref='100')

m3 = glm(mean.CR1.spacers~bottleneck*pfu, data=crispr)
summary(m3)
model_stats(m3)

#### CR2 SPACRS
m1 = glm(mean.CR2.spacers~bottleneck, data=crispr)
m2 = glm(mean.CR2.spacers~pfu, data=crispr)
m3 = glm(mean.CR2.spacers~bottleneck*pfu, data=crispr)

anova(m3, test='LRT')
AICs = AIC(m1, m2, m3)
compare_AICs(AICs)

par(mfrow=c(2,2))
plot(m1)
plot(m2)
plot(m3)

crispr$bottleneck = relevel(crispr$bottleneck, ref='10000000')
crispr$bottleneck = relevel(crispr$bottleneck, ref='1000000')
crispr$bottleneck = relevel(crispr$bottleneck, ref='100000')
crispr$bottleneck = relevel(crispr$bottleneck, ref='10000')
crispr$bottleneck = relevel(crispr$bottleneck, ref='1000')
crispr$bottleneck = relevel(crispr$bottleneck, ref='100')

m3 = glm(mean.CR2.spacers~bottleneck*pfu, data=crispr)
summary(m3)
model_stats(m3)

#### TOTAL SPACERS
m1 = glm(mean.total.spacers~bottleneck, data=crispr)
m2 = glm(mean.total.spacers~pfu, data=crispr)
m3 = glm(mean.total.spacers~bottleneck*pfu, data=crispr)

anova(m3, test='LRT')
AICs = AIC(m1, m2, m3)
compare_AICs(AICs)

par(mfrow=c(2,2))
plot(m1)
plot(m2)
plot(m3)

crispr$bottleneck = relevel(crispr$bottleneck, ref='10000000')
crispr$bottleneck = relevel(crispr$bottleneck, ref='1000000')
crispr$bottleneck = relevel(crispr$bottleneck, ref='100000')
crispr$bottleneck = relevel(crispr$bottleneck, ref='10000')
crispr$bottleneck = relevel(crispr$bottleneck, ref='1000')
crispr$bottleneck = relevel(crispr$bottleneck, ref='100')

m3 = glm(mean.total.spacers~bottleneck*pfu, data=crispr)
summary(m3)
model_stats(m3)


######### 
## Load summary data

crispr_sum = read.csv('./data/PCR/spacer_summary.csv', header=T)
crispr_sum$bottleneck = as.factor(crispr_sum$bottleneck)
crispr_sum = subset(crispr_sum, measurement != 'phage')
cr1_sum = subset(crispr_sum, measurement    == 'CR1')
cr2_sum = subset(crispr_sum, measurement    == 'CR2')
total_sum = subset(crispr_sum, measurement  == 'total')

pd = position_dodge(0.1)

cr1_sum_plot_all = ggplot(data=cr1_sum, aes(x=bottleneck, y=Mean, group=measurement))+
  #geom_bar(stat="identity", size=0.4, position=pd, 
  #         fill='white', colour='grey1')+
  geom_point(stat='identity', size=3)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.4, size=0.7, position=pd)+
  labs(x='Bottleneck size', y=expression(bold("Mean no. of CRISPR1 spacers")))+
  theme_bw()+
  theme(legend.key.size = unit(0.5, 'cm'))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face="bold", size=14))+
  theme(legend.title = element_text(face='bold', size=12))+
  coord_cartesian(ylim = c(0.0, 0.416537112))+
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text=element_text(size=12))
cr1_sum_plot_all

cr2_sum_plot_all <- ggplot(data=cr2_sum, aes(x=bottleneck, y=Mean, group=measurement))+
  #geom_bar(stat="identity", size=0.4, position=pd, 
  #         fill='white', colour='grey1')+
  geom_point(stat='identity', size=3)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.4, size=0.7, position=pd)+
  
  labs(x='Bottleneck size', y=expression(bold("Mean no. of CRISPR2 spacers")))+
  theme_bw()+
  
  theme(legend.key.size = unit(0.5, 'cm'))+
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face="bold", size=14))+
  theme(legend.title = element_text(face='bold', size=12))+
  
  coord_cartesian(ylim = c( 0, 2.3237750))+
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text=element_text(size=12))
  
  
cr2_sum_plot_all

total_sum_plot_all <- ggplot(data=total_sum, aes(x=bottleneck, y=Mean, group=measurement))+
  #geom_bar(stat="identity", size=0.4, position=pd, 
  #         fill='white', colour='grey1')+
  geom_point(stat='identity', size=3)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.2, size=0.7, position=pd)+
  labs(x='Bottleneck size', y=expression(bold("Total CRISPR spacers")))+
  theme_bw()+
  theme(legend.key.size = unit(0.5, 'cm'))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face="bold", size=14))+
  theme(legend.title = element_text(face='bold', size=12))+
  coord_cartesian(ylim = c( 0, 3.4473236))+
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text=element_text(size=12))+
  scale_y_continuous(breaks=c(seq(0,3.5,0.5)))

total_sum_plot_all

comparison_plot_all = ggplot(data=crispr_sum, aes(x=bottleneck, y=Mean, group=measurement))+
  geom_point(stat='identity', size=3, aes(shape=measurement, colour=measurement), position=pd)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper, colour=measurement), width=0.4, size=0.7, position=pd)+
  labs(x='Bottleneck size', y=expression(bold('Number of spacers')))+
  theme_bw()+

  theme(plot.title = element_text(face='bold', hjust=0.5, size = 16))+
  theme(axis.title = element_text(face='bold', size=16))+
  theme(legend.title=element_text(face='bold', size=14))+
  
  scale_shape_discrete(name=c('Spacer type'),
                       breaks=c('CR1', 'CR2', 'total'),
                       labels=c('CRISPR 1', 'CRISPR 2', 'Total number\nof spacers'))+
  scale_colour_discrete(name=c('Spacer type'),
                       breaks=c('CR1', 'CR2', 'total'),
                       labels=c('CRISPR 1', 'CRISPR 2', 'Total number\nof spacers'))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.key.width = unit(1.5, 'cm'))+
  theme(legend.text = element_text(size=12))+
  theme(legend.title.align = 0.5)+
  
  theme(axis.text = element_text(size=12))+
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  
  coord_cartesian(ylim=c(0,3.4473236))

comparison_plot_all



#####################
#### ANALYSIS OF JUST THE CRISPR CLONES #####
#####################

just = read.csv('./data/PCR/just_CRISPR_clones_phage_bn.csv')
just$bottleneck = as.factor(just$bottleneck)
cr1_just = subset(just, mean.CR1.spacers >= 1)
cr2_just = subset(just, mean.CR2.spacers >=1)
total_just = subset(just, mean.total.spacers >=1)
str(just)

##### Summary figures to eyeball the data

# Boxplot to get an impression of the differences between bottleneck treatments
CR2_just_plot_phage_bn = ggplot(cr2_just, aes(y=mean.CR2.spacers, x=bottleneck))+
  geom_boxplot()+
  labs(x='Bottleneck size', y="Spacer number")+
  ggtitle('CRISPR2 spacers in CRISPR clones')+
  theme_bw()+
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 20))+
  theme(axis.title = element_text(face='bold', size=16))+
  
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text = element_text(size=12))+
  
  scale_y_continuous(breaks = c(seq(0,3.5,0.5)))
  
CR2_just_plot_phage_bn

CR1_just_plot_phage_bn = ggplot(cr1_just, aes(y=mean.CR1.spacers, x=bottleneck))+
  geom_boxplot()+
  labs(x='Bottleneck size', y="Spacer number")+
  ggtitle('CRISPR1 spacers in CRISPR clones')+
  theme_bw()+
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 20))+
  theme(axis.title = element_text(face='bold', size=16))+
  
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text = element_text(size=12))
CR1_just_plot_phage_bn

total_just_plot_phage_bn = ggplot(total_just, aes(y=mean.total.spacers, x=bottleneck))+
  geom_boxplot(na.rm = F)+
  labs(x='Bottleneck size', y="Spacer number")+
  #ggtitle('Spacers in both CRISPR loci\nin CRISPR clones')+
  theme_bw()+
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 20))+
  theme(axis.title = element_text(face='bold', size=16))+
  
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text = element_text(size=12))+
  
  scale_y_continuous(breaks = c(seq(0,4.5,0.5)))
total_just_plot_phage_bn

ggsave('CR1_just_plot.png', plot=CR1_just_boxplot, device='png',
       path='./figs/', width=7, height = 5, units=c('in'), dpi=600)  
ggsave('CR2_just_plot.png', plot=CR2_just_boxplot, device='png',
       path='./figs/', width=7, height = 5, units=c('in'), dpi=600)  
ggsave('total_just_plot.png', plot=total_just_boxplot, device='png',
       path='./figs/', width=7, height = 5, units=c('in'), dpi=600) 


##############
### MODELS ###
##############

#### CR1 SPACERS
m1 = glm(mean.CR1.spacers~bottleneck, data=cr1_just)
m2 = glm(mean.CR1.spacers~pfu, data=cr1_just)
m3 = glm(mean.CR1.spacers~bottleneck*pfu, data=cr1_just)

anova(m3, test='LRT')
AICs = AIC(m1, m2, m3)
compare_AICs(AICs)

par(mfrow=c(2,2))
plot(m1)
plot(m2)
plot(m3)

cr1_just$bottleneck = relevel(cr1_just$bottleneck, ref='10000000')
cr1_just$bottleneck = relevel(cr1_just$bottleneck, ref='1000000')
cr1_just$bottleneck = relevel(cr1_just$bottleneck, ref='100000')
cr1_just$bottleneck = relevel(cr1_just$bottleneck, ref='10000')
cr1_just$bottleneck = relevel(cr1_just$bottleneck, ref='1000')
cr1_just$bottleneck = relevel(cr1_just$bottleneck, ref='100')

m3 = glm(mean.CR1.spacers~bottleneck*pfu, data=cr1_just)
summary(m3)
model_stats(m3)

#### CR2 SPACRS
m1 = glm(mean.CR2.spacers~bottleneck, data=just)
m2 = glm(mean.CR2.spacers~pfu, data=just)
m3 = glm(mean.CR2.spacers~bottleneck*pfu, data=just)

anova(m3, test='LRT')
AICs = AIC(m1, m2, m3)
compare_AICs(AICs)

par(mfrow=c(2,2))
plot(m1)
plot(m2)
plot(m3)

just$bottleneck = relevel(just$bottleneck, ref='10000000')
just$bottleneck = relevel(just$bottleneck, ref='1000000')
just$bottleneck = relevel(just$bottleneck, ref='100000')
just$bottleneck = relevel(just$bottleneck, ref='10000')
just$bottleneck = relevel(just$bottleneck, ref='1000')
just$bottleneck = relevel(just$bottleneck, ref='100')

m1 = glm(mean.CR2.spacers~bottleneck, data=just)
summary(m1)
model_stats(m1)

#### TOTAL SPACERS
m1 = glm(mean.total.spacers~bottleneck, data=just)
m2 = glm(mean.total.spacers~pfu, data=just)
m3 = glm(mean.total.spacers~bottleneck*pfu, data=just)

anova(m3, test='LRT')
AICs = AIC(m1, m2, m3)
compare_AICs(AICs)

par(mfrow=c(2,2))
plot(m1)
plot(m2)
plot(m3)

just$bottleneck = relevel(just$bottleneck, ref='10000000')
just$bottleneck = relevel(just$bottleneck, ref='1000000')
just$bottleneck = relevel(just$bottleneck, ref='100000')
just$bottleneck = relevel(just$bottleneck, ref='10000')
just$bottleneck = relevel(just$bottleneck, ref='1000')
just$bottleneck = relevel(just$bottleneck, ref='100')

m1 = glm(mean.total.spacers~bottleneck, data=just)
summary(m1)
model_stats(m1)

#### SUMMARY FIGURES

just_sum = read.csv('./data/PCR/spacer_summary_just_CRISPR.csv', header=T)
just_sum$bottleneck = as.factor(just_sum$bottleneck)
just_sum = subset(just_sum, measurement != 'phage')
cr1_sum_just = subset(just_sum, measurement    == 'CR1')
cr2_sum_just = subset(just_sum, measurement    == 'CR2')
total_sum_just = subset(just_sum, measurement  == 'total')

pd = position_dodge(0.1)

cr1_sum_plot_just = ggplot(data=cr1_sum_just, aes(x=bottleneck, y=Mean, group=measurement))+
  #geom_bar(stat="identity", size=0.4, position=pd, 
  #         fill='white', colour='grey1')+
  geom_point(stat='identity', size=3)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.4, size=0.7, position=pd)+
  labs(x='Bottleneck size', y=expression(bold("Mean no. of CRISPR1 spacers")))+
  theme_bw()+
  theme(legend.key.size = unit(0.5, 'cm'))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face="bold", size=14))+
  theme(legend.title = element_text(face='bold', size=12))+
  
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text=element_text(size=12))

cr1_sum_plot_just

cr2_sum_plot_just <- ggplot(data=cr2_sum_just, aes(x=bottleneck, y=Mean, group=measurement))+
  #geom_bar(stat="identity", size=0.4, position=pd, 
  #         fill='white', colour='grey1')+
  geom_point(stat='identity', size=3)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.4, size=0.7, position=pd)+
  labs(x='Bottleneck size', y=expression(bold("Mean no. of CRISPR2 spacers")))+
  theme_bw()+
  theme(legend.key.size = unit(0.5, 'cm'))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face="bold", size=14))+
  theme(legend.title = element_text(face='bold', size=12))+
  coord_cartesian(ylim=c(0, 2.3237750))+
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text=element_text(size=12))
cr2_sum_plot_just

total_sum_plot_just <- ggplot(data=total_sum_just, aes(x=bottleneck, y=Mean, group=measurement))+
  #geom_bar(stat="identity", size=0.4, position=pd, 
  #         fill='white', colour='grey1')+
  geom_point(stat='identity', size=3)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.2, size=0.7, position=pd)+
  labs(x='Bottleneck size', y=expression(bold("Total CRISPR spacers")))+
  theme_bw()+
  theme(legend.key.size = unit(0.5, 'cm'))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face="bold", size=14))+
  theme(legend.title = element_text(face='bold', size=12))+
  coord_cartesian(ylim=c(0, 3.4473236))+
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text=element_text(size=12))+
  scale_y_continuous(breaks=c(seq(0,3.5,0.5)))

total_sum_plot_just

layer_scales(total_sum_plot_just)$y$range$range

comparison_plot_just = ggplot(data=just_sum, aes(x=bottleneck, y=Mean, group=measurement))+
  geom_point(stat='identity', size=3, aes(shape=measurement, colour=measurement), position=pd)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper, colour=measurement), width=0.4, size=0.7, position=pd)+
  labs(x='Bottleneck size', y=expression(bold('Number of spacers')))+
  theme_bw()+
  
  theme(plot.title = element_text(face='bold', hjust=0.5, size = 16))+
  theme(axis.title = element_text(face='bold', size=16))+
  theme(legend.title=element_text(face='bold', size=14))+
  
  scale_shape_discrete(name=c('Spacer type'),
                       breaks=c('CR1', 'CR2', 'total'),
                       labels=c('CRISPR 1', 'CRISPR 2', 'Total number\nof spacers'))+
  scale_colour_discrete(name=c('Spacer type'),
                        breaks=c('CR1', 'CR2', 'total'),
                        labels=c('CRISPR 1', 'CRISPR 2', 'Total number\nof spacers'))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.key.width = unit(1.5, 'cm'))+
  theme(legend.text = element_text(size=12))+
  theme(legend.title.align = 0.5)+
  
  theme(axis.text = element_text(size=12))+
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)

comparison_plot_just

layer_scales(comparison_plot_just)$y$range$range

all_sum_plots = ggarrange(cr1_sum_plot_all+labs(x='', y='CRISPR1')+ggtitle('All clones'), 
                          cr1_sum_plot_just+labs(x='', y='')+ggtitle("CRISPR clones"),
                          
             cr2_sum_plot_all+labs(x='', y='CRISPR2'), 
             cr2_sum_plot_just+labs(x='', y=''),
             
             total_sum_plot_all+labs(x='', y='Both loci'), 
             total_sum_plot_just+labs(y='', x=''),
             
             ncol=2, nrow=3,
             align='v')

all_sum_plots = annotate_figure(all_sum_plots,
                                bottom=text_grob("Bottleneck size", face='bold', size=18, hjust=0.5, vjust=-0.8),
                                top=text_grob('Mean spacer number at 3 d.p.i. in\nfull bottleneck experiment', 
                                              face='bold', size=18, hjust=0.5, color = 'red'))
quartz()
all_sum_plots

ggsave('test.png', plot=all_sum_plots, device='png',
       path='.', width=8.27, height = 11.69, units=c('in'), dpi=1000)

all_comp_plots = ggarrange(comparison_plot_all+labs(x='')+scale_y_continuous(breaks=c(seq(0,3.5,0.5))),
                           comparison_plot_just+labs(x='', y='')+scale_y_continuous(breaks=c(seq(0,3.5,0.5)))+
                             coord_cartesian(ylim=c(0,3.4473236)),
                           common.legend = T, legend='right')
all_comp_plots = annotate_figure(all_comp_plots,
                                 bottom=text_grob("Bottleneck size", face='bold', size=18, hjust=1, vjust=-0.8))
quartz()
all_comp_plots

ggsave('test.png', plot=all_comp_plots, device='png',
       path='.', width=11.69, height = 6, units=c('in'), dpi=600)
