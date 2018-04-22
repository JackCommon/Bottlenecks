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

# AIC comparison function
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
logit2prob(confint(mod.4))
