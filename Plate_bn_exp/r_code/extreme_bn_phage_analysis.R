library(ggplot2)
library(scales)
library(grid)
library(gridExtra)
library(reshape2)
library(ggpubr)

#########################
#### LOADING DATA #######
#########################

# Change to directory where bottleneck data is stored (commend out depending on if you're using OSX or Ubuntu)
setwd("~/Documents/OneDrive - University of Exeter/Data/Bottlenecks/Plate_bn_exp/") #Work OSX
#setwd("~/OneDrive/Data/Bottlenecks/") #Ubuntu
#setwd("~/OneDrive - University of Exeter/Data/Bottlenecks/Phage_bn_exp/") #Home OSX


# Load in the bacteria-phage count data and remove NAs
extreme = read.csv("./data/counts/counts_master_extreme_bn_csv.csv", header=T)
extreme = subset(extreme, timepoint != '')
#extreme = na.exclude(extreme)

#set the bottleneck treatment as a factor
extreme$bottleneck = as.factor(extreme$bottleneck)
extreme$ID = as.factor(extreme$ID)

# make new columns in the data for log-transformed pfu and cfu counts
# adds one to account for zero values
extreme$pfu               = extreme$pfu + 1
extreme$log.pfu           = log10(extreme$pfu)

extreme$bottleneck = relevel(extreme$bottleneck, ref='5-clone')
extreme$bottleneck = relevel(extreme$bottleneck, ref='monoculture')

########################
### DEFINE FUNCTIONS ###
########################

# A list of new bottleneck treatment names, plus a function, for better graphs
OG_bottleneck = c('monoculture', '5-clone', '50-clone')

bottleneck_names_facet = list(
  'monoculture'      = expression('Monoculture'),
  '5-clone'     = expression('5-clone'),
  '50-clone'    = expression('50-clone')
)

bottleneck_names_legend = c(
  expression('Monoculture'), 
  expression('5-clone'), 
  expression('50-clone')
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

raw_sum = ggplot(aes(y=pfu, x=timepoint), data=extreme)+
  geom_boxplot(aes(colour=bottleneck))+
  labs(x=expression(bold('Timepoint')), y=expression(bold("P.f.u ml"*{}^{-1}*"")))+
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
raw_sum

raw_sum = ggplot(aes(y=pfu, x=timepoint, group=ID), data=extreme)+
  geom_point(aes(colour=bottleneck), position=pd)+
  geom_path(stat='identity', aes(colour=bottleneck), position=pd)+
  labs(x=expression(bold('Timepoint')), y=expression(bold("Log-10 p.f.u ml"*{}^{-1}*"")))+
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
raw_sum

raw_sum = ggplot(aes(y=pfu, x=timepoint, group=ID), data=extreme)+
  geom_point(stat='identity', position=pd, aes(colour=ID))+
  geom_path(stat='identity', position=pd, aes(colour=ID))+
  geom_text(aes(label=ID), hjust=-0.3, vjust=0, position = position_jitter())+
  labs(x=expression(bold('Timepoint')), y=expression(bold("Log-10 p.f.u ml"*{}^{-1}*"")))+
  facet_wrap(~bottleneck, labeller = bottleneck_labeller)+
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()+
  theme(legend.position = 'none')
quartz()
raw_sum

raw_sum = ggplot(aes(y=pfu, x=bottleneck), data=extreme)+
  geom_boxplot(position=pd)+
  labs(x=expression(bold('Timepoint')), y=expression(bold("Log-10 p.f.u ml"*{}^{-1}*"")))+
  facet_wrap(~timepoint)+
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
raw_sum

### Look at cfu and pfu as a function of each other

t3= read.csv('data/all_data_t3.csv', header = T)
t3$pfu = t3$pfu + 1
t3$prop.Resistant = t3$prop.CRISPR + t3$prop.SM

t3$bottleneck = relevel(t3$bottleneck, ref='5-clone')
t3$bottleneck = relevel(t3$bottleneck, ref='monoculture')

m1 = glm(pfu~cfu*bottleneck, data=t3, family=gaussian(link='log'))
m2 = glm(pfu~cfu, data=t3)
m3 = glm(pfu~1, data=t3)
anova(m1, test='LRT')
AICs_t3 = AIC(m3, m2, m1)
AICs_t3
compare_AICs(AICs_t3)

summary(m1)
par(mfrow=c(2,2))
plot(m1)

t3_mono = subset(t3, bottleneck=='monoculture')


slopes_1 = ggplot(aes(x=cfu, y=pfu), data=t3)+
  geom_point()+
  geom_text(aes(label=ID), hjust=0, vjust=-0.3, position = position_jitter())+
  geom_smooth(method='glm', formula=y~x, aes(linetype=bottleneck, colour=bottleneck), se=F)+
  
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  labs(y=expression(bold('pfu ml'*{}^{-1}*'')), x=expression(bold("c.f.u ml"*{}^{-1}*"")))+
  theme_bw()

slopes_2= ggplot(aes(x=cfu, y=pfu), data=t3)+
  geom_point()+
  geom_text(aes(label=ID), hjust=0, vjust=-0.3, position = position_jitter())+
  geom_smooth(method='glm', formula=y~x, se=F)+
  
  facet_wrap(~bottleneck,scales = 'free')+
  
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  labs(y=expression(bold('pfu ml'*{}^{-1}*'')), x=expression(bold("c.f.u ml"*{}^{-1}*"")))+
  theme_bw()

quartz()  
slopes_1
slopes_2

#### PFU as a function of the proportion of different genotypes at t=3
str(t3)

m1 = glm(pfu~bottleneck*prop.CRISPR, data=t3, family=gaussian(link='log'))
m2 = glm(pfu~prop.CRISPR, data=t3, family=gaussian(link='log'))
m3 = glm(pfu~1, data=t3)
anova(m1, test='LRT')
AICs_t3 = AIC(m3, m2, m1)
AICs_t3
compare_AICs(AICs_t3)

t3$bottleneck = relevel(t3$bottleneck, ref='50-clone')
t3$bottleneck = relevel(t3$bottleneck, ref='5-clone')
t3$bottleneck = relevel(t3$bottleneck, ref='monoculture')

summary(m1)
par(mfrow=c(2,2))
plot(m1)

slopes_CRISPR_1 = ggplot(aes(x=prop.CRISPR, y=pfu), data=t3)+
  geom_point()+
  geom_text(aes(label=ID), hjust=0, vjust=-0.3, position = position_jitter())+
  geom_smooth(method='glm', formula=y~x, aes(linetype=bottleneck, colour=bottleneck), se=F)+
  
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  
  labs(y=expression(bold('pfu ml'*{}^{-1}*'')), x='Proportion of CRISPR clones')+
  theme_bw()

slopes_CRISPR_2 = ggplot(aes(x=prop.CRISPR, y=pfu), data=t3)+
  geom_point()+
  geom_text(aes(label=ID), hjust=0, vjust=-0.3, position = position_jitter())+
  geom_smooth(method='glm', formula=y~x, se=F)+
  
  facet_wrap(~bottleneck,scales = 'free')+
  
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+

  labs(y=expression(bold('pfu ml'*{}^{-1}*'')), x="Proportion")+
  ggtitle('CRISPR')+
  theme_bw()+
  theme(plot.title = element_text(face='bold', size='16'))

quartz()
slopes_CRISPR_1
slopes_CRISPR_2

### SM

str(t3)

m1 = glm(pfu~bottleneck*prop.SM, data=t3, family=gaussian(link='log'))
m2 = glm(pfu~prop.SM, data=t3, family=gaussian(link='log'))
m3 = glm(pfu~1, data=t3)
anova(m1, test='LRT')
AICs_t3 = AIC(m3, m2, m1)
AICs_t3
compare_AICs(AICs_t3)

t3$bottleneck = relevel(t3$bottleneck, ref='50-clone')
t3$bottleneck = relevel(t3$bottleneck, ref='5-clone')
t3$bottleneck = relevel(t3$bottleneck, ref='monoculture')

summary(m1)
  par(mfrow=c(2,2))
plot(m1)

slopes_SM_1 = ggplot(aes(x=prop.SM, y=pfu), data=t3)+
  geom_point()+
  geom_text(aes(label=ID), hjust=0, vjust=-0.3, position = position_jitter())+
  geom_smooth(method='glm', formula=y~x, aes(linetype=bottleneck, colour=bottleneck), se=F)+
  
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  
  labs(y=expression(bold('pfu ml'*{}^{-1}*'')), x='Proportion of SM clones')+
  theme_bw()

slopes_SM_2 = ggplot(aes(x=prop.SM, y=pfu), data=t3)+
  geom_point()+
  geom_text(aes(label=ID), hjust=0, vjust=-0.3, position = position_jitter())+
  geom_smooth(method='glm', formula=y~x, se=F)+
  
  facet_wrap(~bottleneck,scales = 'free')+
  
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  
  labs(y=expression(bold('pfu ml'*{}^{-1}*'')), x="Proportion")+
  ggtitle('Surface mutants')+
  theme_bw()+
  theme(plot.title = element_text(face='bold', size='16'))

quartz()
slopes_SM_1
slopes_SM_2

#### Resistant
str(t3)

m1 = glm(pfu~bottleneck*prop.Resistant, data=t3, family=gaussian(link='log'))
m2 = glm(pfu~prop.Resistant, data=t3, family=gaussian(link='log'))
m3 = glm(pfu~1, data=t3)
anova(m1, test='LRT')
AICs_t3 = AIC(m3, m2, m1)
AICs_t3
compare_AICs(AICs_t3)

t3$bottleneck = relevel(t3$bottleneck, ref='50-clone')
t3$bottleneck = relevel(t3$bottleneck, ref='5-clone')
t3$bottleneck = relevel(t3$bottleneck, ref='monoculture')

summary(m1)
par(mfrow=c(2,2))
plot(m1)

slopes_Resistant_1 = ggplot(aes(x=prop.Resistant, y=pfu), data=t3)+
  geom_point()+
  geom_text(aes(label=ID), hjust=0, vjust=-0.3, position = position_jitter())+
  geom_smooth(method='glm', formula=y~x, aes(linetype=bottleneck, colour=bottleneck), se=F)+
  
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  
  labs(y=expression(bold('pfu ml'*{}^{-1}*'')), x='Proportion of Resistant clones')+
  theme_bw()

slopes_Resistant_2 = ggplot(aes(x=prop.Resistant, y=pfu), data=t3)+
  geom_point()+
  geom_text(aes(label=ID), hjust=0, vjust=-0.3, position = position_jitter())+
  geom_smooth(method='glm', formula=y~x, se=F)+
  
  facet_wrap(~bottleneck,scales = 'free')+
  
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  
  labs(y=expression(bold('pfu ml'*{}^{-1}*'')), x="Proportion of Resistant clones")+
  theme_bw()

quartz()
slopes_Resistant_1
slopes_Resistant_2

#### Sensitive
str(t3)

m1 = glm(pfu~bottleneck*prop.Sensitive, data=t3, family=gaussian(link='log'))
m2 = glm(pfu~prop.Sensitive, data=t3, family=gaussian(link='log'))
m3 = glm(pfu~1, data=t3)
anova(m1, test='LRT')
AICs_t3 = AIC(m3, m2, m1)
AICs_t3
compare_AICs(AICs_t3)

t3$bottleneck = relevel(t3$bottleneck, ref='50-clone')
t3$bottleneck = relevel(t3$bottleneck, ref='5-clone')
t3$bottleneck = relevel(t3$bottleneck, ref='monoculture')

summary(m1)
par(mfrow=c(2,2))
plot(m1)

slopes_Sensitive_1 = ggplot(aes(x=prop.Sensitive, y=pfu), data=t3)+
  geom_point()+
  geom_text(aes(label=ID), hjust=0, vjust=-0.3, position = position_jitter())+
  geom_smooth(method='glm', formula=y~x, aes(linetype=bottleneck, colour=bottleneck), se=F)+
  
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  
  labs(y=expression(bold('pfu ml'*{}^{-1}*'')), x='Proportion of Sensitive clones')+
  theme_bw()

slopes_Sensitive_2 = ggplot(aes(x=prop.Sensitive, y=pfu), data=t3)+
  geom_point()+
  geom_text(aes(label=ID), hjust=0, vjust=-0.3, position = position_jitter())+
  geom_smooth(method='glm', formula=y~x, se=F)+
  
  facet_wrap(~bottleneck,scales = 'free')+
  
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  
  labs(y=expression(bold('pfu ml'*{}^{-1}*'')), x="Proportion")+
  ggtitle('Sensitives')+
  theme_bw()+
  theme(plot.title = element_text(face='bold', size='16'))

quartz()
slopes_Sensitive_1
slopes_Sensitive_2

library(cowplot)
phenotype_phage_model_plots = plot_grid(slopes_CRISPR_2, slopes_SM_2,
                            ncol = 1, nrow=2, align = 'v', axis='l', rel_widths = c(.3,.3,.3,3), rel_heights = c(.3,.3,.3,3))
quartz()
phenotype_phage_model_plots

ggsave('plate_exp_pheno_mod_plots.png', phenotype_phage_model_plots, device='png',
       path = '.', width=8.27, height = 9, unit=c('in'), dpi=400)


#### relationship between PFU and total number of spacers in CRISPR clones
slopes_spacers_2 = ggplot(aes(x=mean.just_total.spacers, y=pfu), data=t3)+
  geom_point()+
  geom_text(aes(label=ID), hjust=0, vjust=-0.3, position = position_jitter())+
  geom_smooth(method='glm', formula=y~x, se=F)+
  
  facet_wrap(~bottleneck,scales = 'free')+
  
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  
  labs(y=expression(bold('pfu ml'*{}^{-1}*'')), x="Mean CRISPR1 & 2 spacers")+
  ggtitle('Spacers')+
  theme_bw()+
  theme(plot.title = element_text(face='bold', size='16'))

quartz()
slopes_spacers_2

t3$bottleneck = relevel(t3$bottleneck, ref='50-clone')
t3$bottleneck = relevel(t3$bottleneck, ref='5-clone')
t3$bottleneck = relevel(t3$bottleneck, ref='monoculture')

m1 = glm(pfu~bottleneck*mean.just_total.spacers, data=t3, family=gaussian(link='log'))
m2 = glm(pfu~mean.just_total.spacers, data=t3, family=gaussian(link='log'))
m3 = glm(pfu~1, data=t3)
anova(m1, test='LRT')
AICs_t3 = AIC(m3, m2, m1)
AICs_t3
compare_AICs(AICs_t3)

summary(m1)
confint(m1, level=c(0.95))

phenotype_phage_model_plots = plot_grid(slopes_CRISPR_2, slopes_SM_2, slopes_spacers_2,
                                        ncol = 1, nrow=3, align = 'v', axis='l')
quartz()
phenotype_phage_model_plots

ggsave('plate_exp_pheno_mod_plots_2.png', phenotype_phage_model_plots, device='png',
       path = '.', width=8.27, height = 9 , unit=c('in'), dpi=400)
