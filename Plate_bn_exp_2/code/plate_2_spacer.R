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

OG_spacers <- c("1", "2", "3", "4", ">=5")

spacer_names_facet <- list(
  "1" = "1",
  "2" = "2",
  "3" = "3",
  "4" = "4",
  "5" = "\u2265 5"
)

spacer_names_legend <- c("1", "2", "3", "4", "\u2265 5")

spacer_labeller = function(variable, value) {
  return(spacer_names_facet[value])
}

#### Data ####
spacer <- read.csv('./Plate_bn_exp_2/original_data/spacer_1.csv')
spacer %<>% na.exclude
spacer$Replicate %<>% as.factor()
spacer %<>% filter(Total.spacers > 0)
spacer %<>% select(-Zero)

spacer2 <- select(spacer, ID, Replicate, Bottleneck, Single, Multiple)
spacer2 <- melt(spacer2, measure.vars = c("Single", "Multiple"))

#### Models ignoring replicate ####
spacer <- melt(spacer, measure.vars = c("One", "Two", "Three", "Four", "FiveOrMore", "Single", "Multiple"))
spacer <- plyr::rename(spacer, c("variable"="Spacers"))

mod.null <- glm(value~1, data=spacer, family=binomial)
mod.1 <- glm(value~Spacers, data=spacer, family=binomial)
mod.2 <- glm(value~Bottleneck, data=spacer, family=binomial)
mod.g <- glm(value~Spacers*Bottleneck, data=spacer, family=binomial)

AIC(mod.null, mod.1, mod.2, mod.g) %>% compare_AICs()
anova(mod.g, test="Chisq")
summary(mod.g)

pchisq(mod.g$deviance, mod.g$df.residual)
model.tables(aov(mod.g), "mean")

#### Models including replicate for Single or Multiple spacers ####
spacer2 <- plyr::rename(spacer2, c("variable"="Spacers"))

mod.null <- glm(value~1, data=spacer2, family = binomial)
mod.1 <- glm(value~Spacers, data=spacer2, family = binomial)
mod.2 <- glm(value~Bottleneck, data=spacer2, family = binomial)
mod.3 <- glm(value~Replicate, data=spacer2, family = binomial)
mod.4 <- glm(value~Spacers*Bottleneck, data=spacer2, family = binomial)
mod.5 <- glm(value~Spacers*Replicate, data=spacer2, family = binomial)
mod.6 <- glm(value~Bottleneck*Replicate, data=spacer2, family = binomial)
mod.g <- glm(value~Spacers*Bottleneck*Replicate, data=spacer2, family = binomial)

AIC(mod.null, mod.1, mod.2, mod.3,
    mod.4, mod.5, mod.6, mod.g) %>% compare_AICs()
anova(mod.g, test="Chisq")
pchisq(mod.2$deviance, mod.2$df.residual)
r.squaredLR(mod.g)

summary(mod.g)

model.tables(aov(mod.g), "mean")

spacer2$Replicate %<>% relevel(ref="12")
spacer2$Replicate %<>% relevel(ref="11")
spacer2$Replicate %<>% relevel(ref="10")
spacer2$Replicate %<>% relevel(ref="9")
spacer2$Replicate %<>% relevel(ref="8")
spacer2$Replicate %<>% relevel(ref="7")
spacer2$Replicate %<>% relevel(ref="6")
spacer2$Replicate %<>% relevel(ref="5")
spacer2$Replicate %<>% relevel(ref="4")
spacer2$Replicate %<>% relevel(ref="3")
spacer2$Replicate %<>% relevel(ref="2")
spacer2$Replicate %<>% relevel(ref="1")

spacer2$Bottleneck %<>% relevel(ref="50-clone")
spacer2$Bottleneck %<>% relevel(ref="5-clone")
spacer2$Bottleneck %<>% relevel(ref="1-clone")

spacer2$Spacers %<>% relevel(ref="Multiple")
spacer2$Spacers %<>% relevel(ref="Single")

mod.g <- glm(value~Spacers*Bottleneck*Replicate, data=spacer2, family = binomial)
CopyCIs(mod.g)

#### Look at phage titre as function of spacer acquisition ####
spacer_sum <- read.csv("./Plate_bn_exp_2/summary_data/spacer_summary.csv", header=T) 
spacer_sum$Replicate %<>% as.factor
spacer_sum$Spacer.Number %<>% relevel(ref="Single")
spacer_sum$pfu %<>% +1
spacer_sum$log.pfu <- log10(spacer_sum$pfu)
spacer_sum %<>% select(-X)

mod.null <- glm(pfu~1, data=spacer_sum,
                family=gaussian(link="log"))
mod.1 <- glm(pfu~Spacer.Number, data=spacer_sum,
                family=gaussian(link="log"))

mod.1 <- glmer(log.pfu~Spacer.Number*Bot, data=spacer_sum,
             family=gaussian(link="identity"))

mod.2 <- glm(log.pfu~Spacer.Number+Bottleneck+Replicate, data=spacer_sum,
             family=gaussian(link="identity"))
par(mfrow=c(2,2))
plot(mod.2)

summary(mod.2)
drop1(mod.2)

#### Figures ####
spacer_plot <- ggplot(aes(y=mean, x=Replicate, group=Spacer.Number), data=spacer_sum)+
  geom_bar(stat='identity', aes(fill=Spacer.Number), size=0.5, colour="black", position = position_dodge(1))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2, size=0.7, position = position_dodge(1))+
  # geom_point(stat="identity", position = position_dodge(.5))+
  labs(x='Replicate', y="Proportion")+
  # ggtitle('Figure S4: Proportion of CRISPR clones\nwith different numbers of acquired spacers')+
  theme_bw()+
  
  facet_wrap(~Bottleneck)+
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 12, colour = "red"))+
  theme(axis.title = element_text(face='bold', size=16))+
  theme(strip.text = element_text(size=12))+
  theme(axis.text = element_text(size=12))+
  
  #scale_fill_manual(name=c(""),
  #                  values=c("white", "grey35"))+
  theme(legend.title = element_text(face="bold", size=12))+
  theme(legend.text = element_text(size = 12))+
  theme(legend.position = "top")+
  
  coord_cartesian(ylim=c(0,1))+
  scale_y_continuous(breaks = c(seq(0,1,0.25)))

spacer_plot

# Drop the "single" level
spacer_sum %<>% filter(Spacer.Number == "Single")

phage_plot <- ggplot(aes(y=log.pfu, x=Replicate), data=spacer_sum)+
  geom_bar(stat="identity", width=0.5)+
  labs(x='Replicate', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  #ggtitle('Figure S5: Phage titres at 3 d.p.i.')+
  theme_bw()+
  
  facet_wrap(~Bottleneck)+
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 12, colour = "red"))+
  theme(axis.title = element_text(face='bold', size=16))+
  theme(strip.text = element_text(size=12))+
  
 # scale_y_continuous(breaks=c(seq(0,10,1)))+
  
  theme(axis.text = element_text(size=12))
phage_plot

# Combine the plots so they can be compared
combined <- plot_grid(spacer_plot+labs(x=""), phage_plot,
                   ncol = 1, nrow=2, labels = c('A', 'B'),
                   align='hv', axis = 'b', rel_heights = c(1.2, 1))
combined
