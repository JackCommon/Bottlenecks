### Dilution Experiment - Phage & Host Population Analysis
# Created: 15/6/18 by Jack Common

rm(list=ls())

#### Dependencies ####
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
# Extracts the intercept coefficient (mean) and 95% CIs from the GLM objects, with added functionality to copy those to the clipboard for easier input to summary dataframes

# Enable this if using Ubuntu
#clipboard <- function(x, sep="\t", row.names=FALSE, col.names=FALSE){
#    con <- pipe("xclip -selection clipboard -i", open="w")
#    write.table(x, con, sep=sep, row.names=row.names, col.names=col.names)
#    close(con)
#}

model_stats = function(model){
  sum = coef(model)
  conf = confint(model, level=c(0.95))
  print(c(sum[1]))
  print(c(conf[1,1]))
  print(c(conf[1,2]))
  
  stats = data.frame(conf[1,1], conf[1,2])
  #clipboard(stats)
  clip = pipe('pbcopy', 'w')
  write.table(stats, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  print('Coefficients copied to the clipboard')
  
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

## Compare AIC values for model fit
# This function extracts the AIC for each GLM, and then compares the absolute relative differences for each AIC to the model with the lowest AIC. This acts as a measure of model fit. More can be found at: http://faculty.washington.edu/skalski/classes/QERM597/papers_xtra/Burnham%20and%20Anderson.pdf

compare_AICs = function(df){          # df is a dataframe of AIC values 
  print(df)                           # prints the origina AIC values 
  col_len = length(df[,2])            # extracts the number of number of models
  AIC_min = abs(min(df[,2]))          # finds the minimum AIC value
  for (i in seq(1, col_len, 1)){      # loop through the AIC values and prints the absolute differences from AIC_min
    print( (abs(df[i,2])) - AIC_min)
  }
}

## Pseudo-R2 from residual and null deviance
# Although there are many ways of extracting an R2 from a GLM, two are used here. For GLMs with a normal error distribution, I use a simple calculation of residual.deviance/null.deviance. For logistic regressions (i.e. a binomial error structure), I use McFadden's pseudo-R2 (1-log likelihood(model)/log likelihood(null model)) [see Hosmer et al. (2013) "Applied Logistic Regression" 3rd ed.]

pseudo_R2 = function(test.model, null.model){
  if (test.model$family$family == 'gaussian'){
    print(test.model$formula, quote = F)
    print('McFadden\'s pseudo-R2 for normally-distributed data:', quote=F)
    print( (test.model$deviance/test.model$null.deviance) , quote = F)
  }
  if (test.model$family == 'binomial'){
    print(test.model$formula, quote=F)
    print('McFadden\'s pseudo-R2 for logistic regression:', quote=F)
    print( (1-logLik(test.model)/logLik(null.model) ), quote=F)
  }
  if (test.model$family == 'poisson'){
    print(test.model$formula, quote=F)
    print('McFaddens\'s pseudo-R2:', quote=F)
    print( (test.model$deviance/test.model$df.residual), quote =F)
    print('Ratio of deviance to residual DF (to check fit to (quasi)poisson distribution):', quote=F)
    print( test.model$null.deviance/test.model$deviance, quote=F)
  }
}

### Convert log-odds to probabilities from binomial GLM output
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

## Functions and graphical stuff
# Next I've written some of simple functions that make the graph labels more human-readable. There are also some objects that enable these functions to work properly. *Importantly*, there are different ones for the different experiments, which have been indicated in the comments.

# Lists of the original names (as found in the source data) of the levels of each factor. These are used for the breaks() function in ggplots

## Timepoints are shared across experiments
OG_timepoints = c('t0', 't1', 't2', 't3', 't4', 't5')

timepoint_names_facet = list(
  't0' = '0 d.p.i.',
  't1' = '1 d.p.i.',
  't2' = '2 d.p.i.',
  't3' = '3 d.p.i.',
  't4' = '4 d.p.i.',
  't5' = '5 d.p.i.'
)

timepoint_names_legend = c('0', '1', '2', '3', '4', '5')

timepoint_labeller = function(variable, value) {
  return(timepoint_names_facet[value])
}

pd <- position_dodge(0.1)

#### Load and format data ####
data <- read.csv("./Revisions/original_data/dilution_counts_106phage.csv", header = T)
data <- select(data, -raw, -dilution, -MOI, -X)
data$ID %<>% as.factor()
data %<>% na.exclude
#$log.pfu <- log10(phage$pfu+1)

## Melt data for raw value plots
dataM <- melt(data, measure.vars = c("pfu", "cfu"))
dataM <- plyr::rename(dataM, c("variable"="measurement"))
dataM$treatment %<>% relevel(ref="S")

phageIDs <- c("p1", "p2", "p3",
              "p4", "p5", "p6") %>% 
  rep( (length(dataM$ID)/2)/6 )

hostIDs <- c("h1", "h2", "h3",
             "h4", "h5", "h6") %>% 
  rep( (length(dataM$ID)/2)/6 )

ID2 <- c(phageIDs, hostIDs)

dataM$ID2 <- as.factor(ID2)

## Raw figure ####
treat_names <- list("S"="Small",
                    "L"="Large")
treat_labeller <- function(variable,value){
  return(treat_names[value])
}

raw_plot <- ggplot(aes(y=value+1, x=timepoint, group=ID2), 
                    data=dataM)+
  
  geom_path(stat='identity', 
            aes(colour=measurement, linetype=measurement),
            position=pd)+
  scale_colour_manual(name='Measurement',
                      values =c("black", "grey"),
                      breaks = c("pfu", "cfu"),
                      labels = c("Phage", "Host"))+
  scale_linetype_manual(name='Measurement',
                        values=c(1, 2),
                        breaks = c("pfu", "cfu"),
                        labels = c("Phage", "Host"))+
  geom_point(stat='identity', 
             aes(shape=measurement, fill=measurement), colour="transparent",
             position=pd)+
  scale_fill_manual(name='Measurement',
                    values =c("black", "grey"),
                    breaks = c("pfu", "cfu"),
                    labels = c("Phage", "Host"))+
  scale_shape_manual(name='Measurement',
                     values = c(21,24),
                     breaks = c("pfu", "cfu"),
                     labels = c("Phage", "Host"))+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("P.f.u. ml"*{}^{-1}*"/ C.f.u. ml"*{}^{-1}*"")))+
  ggtitle('')+
  facet_grid(~treatment, labeller = treat_labeller)+
  
  theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5'),
                   labels=c('0', '1', '2', '3', '4', '5'))+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))+
  
  geom_hline(yintercept = 1e+2, linetype=2, colour="red")+
  annotate("text", 1.5, 1e+2, vjust=-1, label="Phage detection limit", colour="red")+
  NULL
quartz()
raw_plot

ggsave("dil_dynamics.png", raw_plot, path="./Revisions/figs/",
       device="png", dpi=300,
       width=30, height=15, units=c("cm"))


