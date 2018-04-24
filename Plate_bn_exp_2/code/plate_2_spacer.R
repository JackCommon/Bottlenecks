### Plate bottleneck experimet 2: Spacer number analysis
# Created 24/4/18: Jack Common

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

#Data
spacer <- read.csv('./Plate_bn_exp_2/original_data/spacer_1.csv')
spacer %<>% na.exclude
spacer$Replicate %<>% as.factor()
spacer %<>% filter(Total.spacers > 0)
spacer %<>% select(-Zero)

spacer2 <- select(spacer, ID, Replicate, Bottleneck, Single, Multiple)
spacer2 <- melt(spacer2, measure.vars = c("Single", "Multiple"))

# Models ignoring replicate
spacer <- melt(spacer, measure.vars = c("One", "Two", "Three", "Four", "FiveOrMore", "Single", "Multiple"))

mod.null <- glm(value~1, data=spacer, family=binomial)
mod.1 <- glm(value~variable, data=spacer, family=binomial)
mod.2 <- glm(value~Bottleneck, data=spacer, family=binomial)
mod.g <- glm(value~variable*Bottleneck, data=spacer, family=binomial)

AIC(mod.null, mod.1, mod.2, mod.g) %>% compare_AICs()
anova(mod.g, test="Chisq")
summary(mod.g)

pchisq(mod.g$deviance, mod.g$df.residual)
model.tables(aov(mod.g), "mean")

# Models including replicate for Single or Multiple spacers

mod.null <- glm(value~1, data=spacer2, family = binomial)
mod.1 <- glm(value~variable, data=spacer2, family = binomial)
mod.2 <- glm(value~Bottleneck, data=spacer2, family = binomial)
mod.3 <- glm(value~Replicate, data=spacer2, family = binomial)
mod.4 <- glm(value~variable*Bottleneck, data=spacer2, family = binomial)
mod.5 <- glm(value~variable*Replicate, data=spacer2, family = binomial)
mod.6 <- glm(value~Bottleneck*Replicate, data=spacer2, family = binomial)
mod.g <- glm(value~variable*Bottleneck*Replicate, data=spacer2, family = binomial)

AIC(mod.null, mod.1, mod.2, mod.3,
    mod.4, mod.5, mod.6, mod.g) %>% compare_AICs()
anova(mod.g, test="Chisq")
pchisq(mod.2$deviance, mod.2$df.residual)

summary(mod.4)

model.tables(aov(mod.2), "mean")
