#### Plate bottleneck experiment 2: phenotype analysis
# Created: 22/04/18 Jack Common

rm(list=ls())

library(ggplot2)
library(scales)
library(reshape2)
library(dplyr)
library(tidyr)
library(magrittr)
library(cowplot)
library(lme4)

#### Functions ####
# AIC comparison function
compare_AICs = function(df){          # df is a dataframe of AIC values 
  print(df)                           # prints the origina AIC values 
  col_len = length(df[,2])            # extracts the number of number of models
  AIC_min = abs(min(df[,2]))          # finds the minimum AIC value
  for (i in seq(1, col_len, 1)){      # loop through the AIC values and prints the absolute differences from AIC_min
    print( (abs(df[i,2])) - AIC_min)
  }
}

## Convert log-odds to probabilities from binomial GLM output
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#### Data - CR-ve ####
pheno <- read.csv('./Revisions/original_data/CR_neg_phenotype.csv', header=T)
pheno <- melt(pheno, measure.vars = c('CRISPR', 'SM', 'Sensitive'))
pheno$Replicate %<>% as.factor()
pheno$Clone %<>% as.factor
pheno$bottleneck %<>% as.factor()
pheno <- plyr::rename(pheno, c("variable"="phenotype"))

#### Models - CR-ve ####
m.null <- glm(value~1, data=pheno, family=binomial)
m1 <- glm(value~bottleneck, data=pheno, family=binomial)
m2 <- glm(value~phenotype, data=pheno, family=binomial)
m.g <- glm(value~bottleneck*phenotype, data=pheno, family=binomial)

par(mfrow=c(2,2))
plot(m1)
plot(m2)
plot(m.g)

# Compare using delta-AIC to get K-L values
AIC(m.null, m1, m2, m.g) %>% compare_AICs()

# And finally a likelihood test
anova(m.g, test="Chisq")
drop1(m.g, test="Chisq")

summary(m.g)
confint(m.g)

# Get the actual probabilities from the model
model.tables(aov(m.g), "mean")

# Confidence intervals
pheno$phenotype %<>% relevel(ref="SM")
pheno$phenotype %<>% relevel(ref="Sensitive")
pheno$phenotype %<>% relevel(ref="CRISPR")

pheno$bottleneck %<>% relevel(ref="9")
pheno$bottleneck %<>% relevel(ref="8")
pheno$bottleneck %<>% relevel(ref="7")
pheno$bottleneck %<>% relevel(ref="6")
pheno$bottleneck %<>% relevel(ref="5")
pheno$bottleneck %<>% relevel(ref="4")
pheno$bottleneck %<>% relevel(ref="3")
pheno$bottleneck %<>% relevel(ref="2")

m.g <- glm(value~bottleneck*phenotype, data=pheno, family=binomial)

logit2prob(confint(m.g))

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

pheno_sum <- read.csv('./Revisions/summary_data/CRneg_pheno_summary.csv')
pheno_sum$bottleneck %<>% as.factor()
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
  
  scale_x_discrete(breaks=sol_OG_bottleneck, labels=sol_bottleneck_names_legend)+
  theme(axis.text = element_text(size=12))+
  
  scale_fill_manual(name='Genotype',
                    breaks = c('CRISPR', "Sensitive", "SM"),
                    labels = c("CRISPR", "Sensitive", "SM"),
                    values = c("#F8766D", "#619CFF", "#00BA38"))+
  NULL
pheno_sum_fig

layer_scales(pheno_sum_fig)$x$range$range
