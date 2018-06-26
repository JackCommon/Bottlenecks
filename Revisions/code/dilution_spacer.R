### Dilution experiment: Spacer analysis
# Created 26/6/18: Jack Common

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
spacer <- read.csv('./Revisions/original_data/dilution_spacers.csv')
spacer %<>% na.exclude
spacer$Replicate %<>% as.factor()
spacer$Clone %<>% as.factor
spacer %<>% filter(Total.spacers > 0)
spacer %<>% select(-Zero)

spacer2 <- select(spacer, Clone, Replicate, treatment, Single, Multiple)
spacer2 <- melt(spacer2, measure.vars = c("Single", "Multiple"))
spacer2 <- plyr::rename(spacer2, c("variable"="spacers"))

#### Models ignoring replicate ####
spacer <- melt(spacer, measure.vars = c("One", "Two", "Three", "Four", "FiveOrMore", "Single", "Multiple"))
spacer <- plyr::rename(spacer, c("variable"="spacers"))

m.null <- glm(value~1, data=spacer, family=binomial)
m1 <- glm(value~spacers, data=spacer, family=binomial)
m2 <- glm(value~treatment, data=spacer, family=binomial)
mg <- glm(value~spacers*treatment, data=spacer, family=binomial)

par(mfrow=c(2,2))
plot(m1)
plot(m2)
plot(mg)

AIC(m.null, m1, m2, mg) %>% compare_AICs()
anova(mg, test="Chisq")
summary(mg)

pchisq(mg$deviance, mg$df.residual)
model.tables(aov(mg), "mean")

#### Models including replicate for Single or Multiple spacers ####

m.null <- glmer(value~1+(1|Replicate), data=spacer2, family = binomial)
m1 <- glmer(value~spacers+(1|Replicate), data=spacer2, family = binomial)
m2 <- glmer(value~treatment+(1|Replicate), data=spacer2, family = binomial)
#m3 <- glmer(value~Replicate+(1|Replicate), data=spacer2, family = binomial)
m4 <- glmer(value~spacers*treatment+(1|Replicate), data=spacer2, family = binomial)
#m5 <- glmer(value~spacers*Replicate+(1|Replicate), data=spacer2, family = binomial)
#m6 <- glmer(value~treatment*Replicate+(1|Replicate), data=spacer2, family = binomial)
m.global <- glmer(value~spacers*treatment+(1|Replicate), data=spacer2, family = binomial)

AIC(m.null, m1, m2, #m3,
    m4, m.global) %>% compare_AICs()
anova(m.null, m1, m2, m4, m.global, test="Chisq")

summary(m.global)


#### Figures -  NOT UPDATED FOR DILUTION DATA ####
spacer_plot <- ggplot(aes(y=mean, x=Replicate, group=Spacer.Number), data=spacer_sum)+
  geom_bar(stat='identity', aes(fill=Spacer.Number), size=0.5, position = position_dodge(1))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2, size=0.7, position = position_dodge(1))+
  # geom_point(stat="identity", position = position_dodge(.5))+
  labs(x='Replicate', y="Frequency")+
  # ggtitle('Figure S4: Proportion of CRISPR clones\nwith different numbers of acquired spacers')+
  theme_bw()+
  
  facet_wrap(~Bottleneck)+
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 12, colour = "red"))+
  theme(axis.title = element_text(face='bold', size=16))+
  theme(strip.text = element_text(size=12))+
  theme(axis.text = element_text(size=12))+
  
  scale_fill_discrete(name=c("Spacer number"))+
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

ggsave("Spacer_Titre_comp.png", combined, path="../figs/",
       device="png", dpi=300,
       width=25, height=18, units=c("cm"))
 