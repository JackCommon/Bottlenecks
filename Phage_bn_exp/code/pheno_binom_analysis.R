library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(magrittr)

setwd("./OneDrive/Data/Bottlenecks/Phage_bn_exp/")

pheno <- read.csv("./data/phenotype/phenotype_prop.csv")
pheno$Bottleneck %<>% as.factor
pheno %<>% select(-Phenotype)
pheno <- melt(pheno, measure.vars = c('CRISPR', 'SM', 'Sensitive'))

# Interested in the proportion of each phenotype in each bottleneck treatment, so the models use a binomial distribution with a logit link
mod.null <- glm(value~1, data=pheno, family=binomial)
mod.1 <- glm(value~Bottleneck, data=pheno, family=binomial)
mod.2 <- glm(value~variable, data=pheno, family=binomial)
mod.g <- glm(value~Bottleneck*variable, data=pheno, family=binomial)

# Compare using delta-AIC to get K-L values
AIC(mod.null, mod.1, mod.2, mod.g) %>% compare_AICs()

# And check R2 
pseudo_R2(mod.1, mod.null)
pseudo_R2(mod.2, mod.null)
pseudo_R2(mod.g, mod.null)

# Likelihood test
anova(mod.g, test="LRT")

# Summarise the global model
summary(mod.g)

# Check dispersion against the Chis-squared distribution
pchisq(summary(mod.g)$dispersion * mod.g$df.residual, mod.g$df.residual, lower=F)

# Get coefficients from the model
model.tables(aov(mod.g), "mean")

pheno$Bottleneck %<>% relevel(., ref="9")
pheno$Bottleneck %<>% relevel(., ref="8")
pheno$Bottleneck %<>% relevel(., ref="7")
pheno$Bottleneck %<>% relevel(., ref='6')
pheno$Bottleneck %<>% relevel(., ref='5')
pheno$Bottleneck %<>% relevel(., ref='4')
pheno$Bottleneck %<>% relevel(., ref='3')
pheno$Bottleneck %<>% relevel(., ref='2')

pheno$variable %<>% relevel(ref="CRISPR")
pheno$variable %<>% relevel(ref="SM")
pheno$variable %<>% relevel(ref="Sensitive")

mod.g <- glm(value~Bottleneck*variable, data=pheno, family=binomial)
logit2prob(confint(mod.g))
confint(mod.g)

exp2_pheno_sum <- read.csv('data/phenotype/phenotype_binom_summary.csv')
exp2_pheno_sum$bottleneck %<>% as.factor()
exp2_pheno_sum %<>% filter(., variable!='SM')
pd = position_dodge(0.1)

exp2_pheno_fig_binom<- ggplot(aes(y=mean, x=bottleneck, group=variable), data=exp2_pheno_sum)+
  geom_point(aes(shape=variable, colour=variable), stat='identity', size=5, position=pd)+
  geom_errorbar(aes(ymin=lower, ymax=upper, colour=variable), width=0.3, size=0.7, position=pd)+
  
  labs(x='Bottleneck strength', y='Relative frequency')+
  
  theme_bw()+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face='bold', size=16))+
  
  scale_x_discrete(breaks=sol_OG_bottleneck, labels=sol_bottleneck_names_legend)+
  theme(axis.text = element_text(size=12))+
  
  scale_color_discrete(name='Genotype',
                       breaks=c('CRISPR', 'Sensitive'),
                       labels=c('CRISPR', 'Sensitive'))+
  scale_shape_discrete(name='Genotype',
                       breaks=c('CRISPR', 'Sensitive'),
                       labels=c('CRISPR', 'Sensitive'))+
  theme(legend.title = element_text(face='bold', size=16))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(1, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=12))


exp2_pheno_fig + ggtitle("Normal Dist data")
exp2_pheno_fig_binom + ggtitle("Binomial Dist data")
