#### Dilution experiment: phenotype analysis
# Created: 22/04/18 Jack Common

rm(list=ls())

library(ggplot2)
library(scales)
library(reshape2)
library(dplyr)
library(tidyr)
library(magrittr)
library(cowplot)

#### AIC comparison function ####
compare_AICs = function(df){          # df is a dataframe of AIC values 
  print(df)                           # prints the origina AIC values 
  col_len = length(df[,2])            # extracts the number of number of models
  AIC_min = abs(min(df[,2]))          # finds the minimum AIC value
  for (i in seq(1, col_len, 1)){      # loop through the AIC values and prints the absolute differences from AIC_min
    print( (abs(df[i,2])) - AIC_min)
  }
}

### Convert log-odds to probabilities from binomial GLM output

# Data
pheno <- read.csv('./Plate_bn_exp_2/original_data/phenotype_1.csv', header=T)
pheno <- melt(pheno, measure.vars = c('CRISPR', 'SM', 'Sensitive'))
pheno$Replicate %<>% as.factor()

# Binomial model
mod.null <- glm(value~1, data=pheno, family=binomial)
mod.1 <- glm(value~bottleneck, data=pheno, family=binomial)
mod.2 <- glm(value~variable, data=pheno, family=binomial)
mod.g <- glm(value~bottleneck*variable, data=pheno, family=binomial)

# Compare using delta-AIC to get K-L values
AIC(mod.null, mod.1, mod.2, mod.g) %>% compare_AICs()

# And finally a likelihood test
anova(mod.g, test="Chisq")

summary(mod.g)
confint(mod.g)

# Get the actual probabilities from the model
model.tables(aov(mod.g), "mean")

# Confidence intervals
pheno$variable %<>% relevel(ref="SM")
pheno$variable %<>% relevel(ref="Sensitive")
pheno$variable %<>% relevel(ref="CRISPR")

pheno$bottleneck %<>% relevel(ref="50-clone")
pheno$bottleneck %<>% relevel(ref="5-clone")
pheno$bottleneck %<>% relevel(ref="1-clone")

mod.g <- glm(value~bottleneck*variable, data=pheno, family=binomial)

logit2prob(confint(mod.g))

### Summary figures
OG_genotypes = c('prop.CRISPR', 'prop.SM', 'prop.Sensitive')
genotype_names_facet = list(
  'prop.CRISPR' = 'CRISPR',
  'prop.Sensitive' = 'SM',
  'prop.SM' = 'Sensitive'
)

genotype_names_legend = c('CRISPR', 'SM', 'Sensitive')

genotype_labeller = function(variable, value) {
  return(genotype_names_facet[value])
}

pheno_sum <- read.csv('./Plate_bn_exp_2/summary_data/phenotype_summary.csv')
pheno_sum %<>% na.exclude()
pheno_sum$variable %<>% relevel(., ref='Sensitive')
pheno_sum$variable %<>% relevel(., ref='SM')
pheno_sum$variable %<>% relevel(., ref='CRISPR')

pheno_sum_fig <- ggplot(aes(y=mean, x=bottleneck, group=variable), data=pheno_sum)+
  geom_bar(aes(fill=variable), stat='identity', size=3.5, position=position_dodge(1))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.3, size=0.7, position=position_dodge(1))+
  
  labs(x='Bottleneck', y='Proportion')+
  
  theme_bw()+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face='bold', size=14))+
  theme(legend.title = element_text(face='bold', size=12))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = "right")+
  theme(legend.key.width = unit(1, 'cm'))+
  theme(legend.key.height = unit(0.8, 'cm'))+
  theme(legend.text = element_text(size=11))+
  
  scale_x_discrete(breaks=plate_OG_bottleneck, labels=plate_bottleneck_names_legend)+
  theme(axis.text = element_text(size=12))+
  
  scale_fill_manual(name='Genotype',
                    breaks = c('CRISPR', "Sensitive", "SM"),
                    labels = c("CRISPR", "Sensitive", "SM"),
                    values = c("#F8766D", "#619CFF", "#00BA38"))
pheno_sum_fig
