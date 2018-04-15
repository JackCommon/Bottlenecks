setwd('~/Documents/OneDrive - University of Exeter/Data/Bottlenecks/Phage_bn_exp/')
#setwd('~/OneDrive - University of Exeter/Data/Bottlenecks/Phage_bn_exp/')
list.files()

# Make the names of genotypes better when graphing
OG_bottleneck = c('100', '1000', '10000', '100000', '1000000','10000000', '100000000', '1000000000')
OG_timepoints = c('t0', 't1', 't2', 't3')
OG_genotypes = c('prop.CRISPR', 'prop.SM', 'prop.Sensitive')

genotype_names_facet = list(
  'prop.CRISPR' = 'CRISPR',
  'prop.Sensitive' = 'SM',
  'prop.SM' = 'Sensitive'
)

genotype_names_legend = c('CRISPR', 'SM', 'Sensitive')

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

genotype_labeller = function(variable, value) {
  return(genotype_names_facet[value])
}


pheno = read.csv('data/phenotype/phenotype_all.csv', header=T)

pheno$Bottleneck = as.factor(pheno$Bottleneck)
pheno$Replicate = as.factor(pheno$Replicate)
pheno$Clone = as.factor(pheno$Clone)

temp1 = subset(pheno, Bottleneck == '1000000000')
temp1 = subset(pheno, Bottleneck == '100000000')
temp1 = subset(pheno, Bottleneck == '10000000')
temp1 = subset(pheno, Bottleneck == '1000000')
temp1 = subset(pheno, Bottleneck == '100000')
temp1 = subset(pheno, Bottleneck == '10000')
temp1 = subset(pheno, Bottleneck == '1000')
temp1 = subset(pheno, Bottleneck == '100')

temp2 = subset(temp1, Replicate == '6')
temp2 = subset(temp1, Replicate == '5')
temp2 = subset(temp1, Replicate == '4')
temp2 = subset(temp1, Replicate == '3')
temp2 = subset(temp1, Replicate == '2')
temp2 = subset(temp1, Replicate == '1')

temp3 = subset(temp2, Phenotype == 'Sensitive')
length(temp3$Phenotype)
temp3 = subset(temp2, Phenotype == 'CRISPR')
length(temp3$Phenotype)
temp3 = subset(temp2, Phenotype == 'SM')
length(temp3$Phenotype)

##### PLOTS
# Load in the data of the phenotypes in each treatment
pheno = read.csv('./data/phenotype/phenotype_results.csv', header=T)
pheno = na.exclude(pheno)
str(pheno)

pheno_m = melt(pheno, id.vars = c('bottleneck', 'repeat.'), measure.vars = c('prop.CRISPR', 'prop.SM', 'prop.Sensitive'))
pheno_m$bottleneck = as.factor(pheno_m$bottleneck)
pheno_m$repeat. = as.factor(pheno_m$repeat.)
pheno_m$value = as.numeric(pheno_m$value)
str(pheno_m)

pd = position_dodge(0.1)

CRISPR_plot <- ggplot(data=subset(pheno_m, variable=='prop.CRISPR'), aes(x=bottleneck, y=value))+
  geom_boxplot(position=pd)+
  #geom_point(stat='identity', size=3)+
  #geom_path(stat="identity")+
  facet_wrap(~variable, labeller = genotype_labeller)+
  labs(x='Bottleneck size', y='Proportion')+
  #ggtitle("Phenotype summary")+
  theme_bw()+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text = element_text(size=12))+
  theme(strip.text = element_text(size=14))

SM_plot <- ggplot(data=subset(pheno_m, variable=='prop.SM'), aes(x=bottleneck, y=value))+
  geom_boxplot(position=pd)+
  #geom_point(stat='identity', size=3)+
  #geom_path(stat="identity")+
  facet_wrap(~variable, labeller = genotype_labeller)+
  labs(x='Bottleneck size', y='Proportion')+
  #ggtitle("Phenotype summary")+
  theme_bw()+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text = element_text(size=12))+
  theme(strip.text = element_text(size=14))

Sensitive_plot <- ggplot(data=subset(pheno_m, variable=='prop.Sensitive'), aes(x=bottleneck, y=value))+
  geom_boxplot(position=pd)+
  #geom_point(stat='identity', size=3)+
  #geom_path(stat="identity")+
  facet_wrap(~variable, labeller = genotype_labeller)+
  labs(x='Bottleneck size', y='Proportion')+
  #ggtitle("Phenotype summary")+
  theme_bw()+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text = element_text(size=12))+
  theme(strip.text = element_text(size=14))

pheno_plot_1 <- ggplot(data=pheno_m, aes(x=bottleneck, y=value))+
  geom_boxplot(position=pd)+
  #geom_point(stat='identity', size=3)+
  #geom_path(stat="identity")+
  facet_wrap(~variable, labeller = genotype_labeller)+
  labs(x='Bottleneck size', y='Proportion')+
  #ggtitle("Phenotype summary")+
  theme_bw()+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text = element_text(size=12))+
  theme(strip.text = element_text(size=14))

pheno_plot_2 <- ggplot(data=pheno_m, aes(x=bottleneck, y=value))+
  geom_boxplot(aes(fill=variable), na.rm=F)+
  labs(x='Bottleneck size', y='Proportion')+
  #ggtitle("Phenotype summary")+
  theme_bw()+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face='bold', size=16))+
  
  scale_x_discrete(breaks=OG_bottleneck, labels=bottleneck_names_legend)+
  theme(axis.text = element_text(size=12))+
  
  scale_fill_discrete(name='Genotype',
                      breaks=OG_genotypes,
                      labels=genotype_names_legend)+
  theme(legend.title = element_text(face='bold', size=16))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(1, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=12))


quartz()
CRISPR_plot
SM_plot
Sensitive_plot
quartz()
pheno_plot_1
pheno_plot_2

ggsave('phenotype_sum_facet.png', plot=pheno_plot_1, device='png',
       path='./figs/', width=12, height = 4.6, units=c('in'), dpi=600)  
ggsave('phenotype_sum_alltogether.png', plot=pheno_plot_2, device='png',
       path='./figs/', width=12, height = 8, units=c('in'), dpi=600)  
  