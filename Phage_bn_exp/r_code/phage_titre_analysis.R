##################
##### PACKAGES ###
##################
install.packages('gridExtra')

library(ggplot2)
library(scales)
library(grid)
library(gridExtra)
library(reshape2)
library(cowplot)

#########################
#### LOADING DATA #######
#########################

# Change to directory where bottleneck data is stored (commend out depending on if you're using OSX or Ubuntu)
setwd("~/Documents/OneDrive - University of Exeter/Data/Bottlenecks/Phage_bn_exp/") #Work OSX
#setwd("~/OneDrive/Data/Bottlenecks/") #Ubuntu
#setwd("~/OneDrive - University of Exeter/Data/Bottlenecks/Phage_bn_exp/") #Home OSX


# Load in the bacteria-phage count data and remove NAs
phage = read.csv("./data/counts/counts_master_csv.csv", header=T)
phage = na.exclude(phage)

#set the bottleneck treatment as a factor
phage$bottleneck = as.factor(phage$bottleneck)
phage$ID = as.factor(phage$ID)

# make new columns in the data for log-transformed pfu and cfu counts
# adds one to account for zero values
phage$pfu               = phage$pfu + 1
phage$log.pfu           = log10(phage$pfu)

########################
### DEFINE FUNCTIONS ###
########################

# A list of new bottleneck treatment names, plus a function, for better graphs
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


OG_timepoints = c('t0', 't1', 't2', 't3', 't4', 't5')

timepoint_names_facet = list(
  't0' = '0 d.p.i.',
  't1' = '1 d.p.i.',
  't2' = '2 d.p.i.',
  't3' = '3 d.p.i.',
  't4' = '4 d.p.i.',
  't5' = '5 d.p.i'
)

timepoint_names_legend = c('0', '1', '2', '3', '4', '5')


#Make the facet labels for the graphs better
bottleneck_labeller = function(variable, value) {
  return(bottleneck_names_facet[value])
}

timepoint_labeller = function(variable, value) {
  return(timepoint_names_facet[value])
}

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


##### RAW DATA FIGS #####

pd = position_dodge(0.1)

raw_sum = ggplot(aes(y=pfu, x=timepoint), data=phage)+
  geom_boxplot(aes(colour=bottleneck))+
  labs(x=expression(bold('Timepoint')), y=expression(bold("P.f.u ml"*{}^{-1}*"")))+
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
raw_sum

raw_sum = ggplot(aes(y=pfu, x=timepoint, group=ID), data=phage)+
  geom_point(aes(colour=bottleneck), position=pd)+
  geom_path(stat='identity', aes(colour=bottleneck), position=pd)+
  labs(x=expression(bold('Timepoint')), y=expression(bold("Log-10 p.f.u ml"*{}^{-1}*"")))+
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
raw_sum

raw_sum = ggplot(aes(y=pfu, x=timepoint, group=ID), data=phage)+
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  labs(x=expression(bold('Timepoint')), y=expression(bold("Log-10 p.f.u ml"*{}^{-1}*"")))+
  facet_wrap(~bottleneck, labeller = bottleneck_labeller)+
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
quartz()
raw_sum

raw_sum = ggplot(aes(y=pfu, x=bottleneck), data=subset(phage, timepoint=='t3'))+
  geom_boxplot(position=pd)+
  labs(x=expression(bold('Bottleneck')), y=expression(bold("Log-10 p.f.u ml"*{}^{-1}*"")))+
  #facet_wrap(~timepoint)+
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
quartz()
raw_sum

## Just have a look at the data for the 1:100 bottleneck treatment
phage_108 = subset(phage, bottleneck == '100000000')
str(phage_108)
phage_108 = na.exclude(phage_108)

raw_sum_108 = ggplot(aes(y=pfu, x=timepoint, group=ID), data=phage_108)+
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  labs(x=expression(bold('Timepoint')), y=expression(bold("Log-10 p.f.u ml"*{}^{-1}*"")))+
  ggtitle('Phage titers in 1:100 bottleneck treatment')+
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))
raw_sum_108


##### MODELS

mod_1 = glm(log.pfu~bottleneck, data=phage)
mod_2 = glm(log.pfu~timepoint, data=phage)
mod_3 = glm(log.pfu~bottleneck*timepoint, data=phage)

par(mfrow=c(2,2))
plot(mod_3)
summary(mod_3)
model_stats(mod_3)

anova(mod_1, mod_2, mod_3, test='LRT')
mod_AICs = AIC(mod_1, mod_2, mod_3)
mod_AICs
compare_AICs(mod_AICs)

phage$bottleneck = relevel(phage$bottleneck, ref="10000000")
phage$bottleneck = relevel(phage$bottleneck, ref='1000000')
phage$bottleneck = relevel(phage$bottleneck, ref='100000')
phage$bottleneck = relevel(phage$bottleneck, ref='10000')
phage$bottleneck = relevel(phage$bottleneck, ref='1000')
phage$bottleneck = relevel(phage$bottleneck, ref='100')

phage$timepoint = relevel(phage$timepoint, ref='t5')
phage$timepoint = relevel(phage$timepoint, ref='t4')
phage$timepoint = relevel(phage$timepoint, ref='t3')
phage$timepoint = relevel(phage$timepoint, ref='t2')
phage$timepoint = relevel(phage$timepoint, ref='t1')
phage$timepoint = relevel(phage$timepoint, ref='t0')

mod_3 = glm(log.pfu~bottleneck*timepoint, data=phage)
summary(mod_3)
model_stats(mod_3)


##### SUMMARY PLOTS

phage_sum = read.csv(file.choose(), header=T)
phage_sum = subset(phage_sum, log.pfu.CI =! 0)
phage_sum$bottleneck = as.factor(phage_sum$bottleneck)

all_line_plot = ggplot(data=phage_sum, aes(x=timepoint, y=log.pfu.mean, group=bottleneck))+
  geom_errorbar(aes(ymin=log.pfu.lower, ymax=log.pfu.upper, colour=bottleneck), width=0.1, size=0.7, position=pd)+
  geom_point(stat="identity", size = 5, aes(colour=bottleneck), position=pd)+
  geom_path(stat='identity', aes(colour=bottleneck), position=pd)+
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log-10 p.f.u ml"*{}^{-1}*"")))+
  ggtitle("")+
  theme_bw()+
  scale_colour_discrete(name='Bottleneck\nsize',
                        breaks = OG_bottleneck,
                        labels = bottleneck_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 20))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5'),
                   labels=c('0', '1', '2', '3', '4', '5'))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

quartz()
t3_facet_plot = ggplot(data=subset(phage_sum, timepoint=='t3'), aes(x=bottleneck, y=log.pfu.mean))+
  geom_errorbar(aes(ymin=log.pfu.lower, ymax=log.pfu.upper), width=0.1, size=0.7, position=pd)+
  geom_point(stat="identity", size = 4, position=pd)+
  labs(x='Bottleneck size', y=expression(bold("Log-10 p.f.u ml"*{}^{-1}*"")))+
  ggtitle("")+
  #facet_wrap(~timepoint, labeller = timepoint_labeller)+
  theme_bw()+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 20))+
  theme(axis.title = element_text(face="bold", size=18))+
  theme(strip.text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  scale_x_discrete(breaks = OG_bottleneck, labels=bottleneck_names_legend)+
  scale_y_continuous(breaks=c(seq(0,5,0.5)))

all_facet_plot = ggplot(data=subset(phage_sum, timepoint!='t0'), aes(x=bottleneck, y=log.pfu.mean))+
  geom_errorbar(aes(ymin=log.pfu.lower, ymax=log.pfu.upper), width=0.1, size=0.7, position=pd)+
  geom_point(stat="identity", size = 3, position=pd)+
  labs(x='Bottleneck size', y=expression(bold("Log-10 p.f.u ml"*{}^{-1}*"")))+
  ggtitle("")+
  facet_wrap(~timepoint, labeller = timepoint_labeller)+
  theme_bw()+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 20))+
  theme(axis.title = element_text(face="bold", size=18))+
  theme(strip.text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  scale_x_discrete(breaks = OG_bottleneck, labels=bottleneck_names_legend)

quartz()
all_facet_plot

quartz()
ggplot(data=phage_sum, aes(x=timepoint, y=log.pfu.mean, group=bottleneck))+
  geom_errorbar(aes(ymin=log.pfu.lower, ymax=log.pfu.upper), width=0.1, size=0.7, position=pd)+
  geom_point(stat="identity", aes(size = 0.3), position=pd)+
  geom_path(stat='identity', position=pd)+
  labs(x='Timepoint', y=expression(bold("Log-10 p.f.u ml"*{}^{-1}*"")))+
  ggtitle("")+
  facet_wrap(~bottleneck, labeller = bottleneck_labeller)+
  theme_bw()+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face="bold", size=14))+
  theme(legend.position = 'none')+
  theme(strip.text = element_text(size=14))+
  theme(axis.text = element_text(size=12))

all_line_plot
t3_facet_plot

ggsave('phage_bn_line_titre_summary.png', plot=all_line_plot, device='png',
       path='./figs/', width=12, height = 8, units=c('in'), dpi=600)
ggsave('phage_bn_t3_titre_summary.png', plot=t3_facet_plot, device='png',
       path='./figs/', width=10, height = 8, units=c('in'), dpi=600)
ggsave('phage_bn_all_titre_summary.png', plot=all_facet_plot, device='png',
       path = './figs/', width = 10, height = 6, units=c('in'), dpi=600)



