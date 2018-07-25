#### Dilution experiment: phenotype analysis
# Created: 24/06/18 Jack Common

rm(list=ls())

library(ggplot2)
library(scales)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)
library(cowplot)
library(lme4)

#### Functions ####
compare_AICs = function(df){          # df is a dataframe of AIC values 
  print(df)                           # prints the origina AIC values 
  col_len = length(df[,2])            # extracts the number of number of models
  AIC_min = abs(min(df[,2]))          # finds the minimum AIC value
  for (i in seq(1, col_len, 1)){      # loop through the AIC values and prints the absolute differences from AIC_min
    print( (abs(df[i,2])) - AIC_min)
  }
}

### Convert log-odds to probabilities from binomial GLM output
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

CopyCIs <- function(model){
  conf = logit2prob(confint(model, level=c(0.95)))
  print(c(conf[1,1]))
  print(c(conf[1,2]))
  CIs <- data.frame(conf[1,1], conf[1,2])
  
  #clipboard(CIs)
  clip = pipe('pbcopy', 'w')
  write.table(CIs, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  print('Probability 95% CIs copied to the clipboard')
  
}

# Data ####
pheno <- read.csv('./Revisions/original_data/dilution_pheno_comp.csv', header=T)
pheno <- melt(pheno, measure.vars = c('CRISPR', 'SM', 'Sensitive'))
pheno$Replicate %<>% as.factor()
pheno$Clone %<>% as.factor
pheno$treatment %<>% relevel(ref="6")
pheno$treatment %<>% relevel(ref="L")
pheno$treatment %<>% relevel(ref="4")
pheno$treatment %<>% relevel(ref="S")
pheno %<>% plyr::rename(c("variable"="phenotype"))

# Binomial GLMs ####
m.null <- glm(value~1, data=pheno, family=binomial)
m1 <- glm(value~treatment, data=pheno, family=binomial)
m2 <- glm(value~phenotype, data=pheno, family=binomial)
m.global <- glm(value~treatment*phenotype, data=pheno, family=binomial)
me.global <- glmer(value~treatment*phenotype+(1|Replicate), data=pheno, family=binomial)

par(mfrow=c(2,2))
plot(m.null)
plot(m1)
plot(m2)
plot(m.global)
plot(me.global)

par(mfrow=c(1,1))
sresid <- resid(me.global, type="pearson")
hist(sresid)

# Compare using AIC and ANOVA
AIC(m.null, m1, m2, m.global) %>% compare_AICs()
anova(m.null, m1, m2, m.global, test="Chisq")

summary(m.global)
confint(mod.g)

# Get the actual probabilities from the model
model.tables(aov(m.global), "mean")

# Confidence intervals
pheno$phenotype %<>% relevel(ref="SM")
pheno$phenotype %<>% relevel(ref="Sensitive")
pheno$phenotype %<>% relevel(ref="CRISPR")

pheno$treatment %<>% relevel(ref="6")
pheno$treatment %<>% relevel(ref="4")
pheno$treatment %<>% relevel(ref="L")
pheno$treatment %<>% relevel(ref="S")

mod.g <- glm(value~treatment*phenotype, data=pheno, family=binomial)

logit2prob(confint(mod.g))

### Summary figures ####
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

pheno_sum <- read.csv('./Revisions/summary_data/dilution_pheno_summary.csv')
#pheno_sum %<>% na.exclude()
pheno_sum$phenotype %<>% relevel(., ref='SM')
pheno_sum$phenotype %<>% relevel(., ref='Sensitive')
pheno_sum$phenotype %<>% relevel(., ref='CRISPR')
pheno_sum$treatment %<>% relevel(ref="L")
pheno_sum$treatment %<>% relevel(ref="6")
pheno_sum$treatment %<>% relevel(ref="S")
pheno_sum$treatment %<>% relevel(ref="4")

pheno_sum_fig <- ggplot(aes(y=mean, x=treatment, group=phenotype), data=pheno_sum)+
  geom_bar(aes(fill=phenotype), stat='identity', size=3.5, position=position_dodge(1))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.3, size=0.7, position=position_dodge(1))+
  
  labs(x='Bottleneck', y='Proportion')+
  ggtitle("Proportion of different immune phenotypes compared between\nculture & phage bottleneck experiment vs. dilution experiment\n(10^-6 phage)")+
  
  theme_bw()+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face='bold', size=14))+
  theme(legend.title = element_text(face='bold', size=12))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = "right")+
  theme(legend.key.width = unit(1, 'cm'))+
  theme(legend.key.height = unit(0.8, 'cm'))+
  theme(legend.text = element_text(size=11))+
  theme(axis.text = element_text(size=12))+
  
  scale_fill_manual(name='Genotype',
                    breaks = c('CRISPR', "Sensitive", "SM"),
                    labels = c("CRISPR", "Sensitive", "SM"),
                    values = c("#F8766D","#00BA38", "#619CFF"))+
  scale_x_discrete(breaks=c("4", "S", "6", "L"),
                   labels=c(expression('10'^-4*''), "Small", expression("10"^-6*""), "Large"))+
  coord_cartesian(ylim=c(0,1))+
  NULL
pheno_sum_fig

ggsave("dil_exp_pheno_comp.png", pheno_sum_fig, path="./Revisions/figs",
       device="png", dpi=300, width=22, height=15, unit=c("cm"))
