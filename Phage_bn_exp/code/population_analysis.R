### Culture bottleneck experiment - population analysis ####
# Created: 18/6/18 by Jack Common

rm(list=ls())

#### Dependencies ####
library(ggplot2)
library(scales)
library(reshape2)
library(plyr)
library(dplyr)
library(broom)
library(tidyr)
library(magrittr)
library(cowplot)
library(lme4)
library(MuMIn)

#### Stats functions ####
stats <- data.frame(coef=rep(0,16), lower=rep(0,16), upper=rep(0,16))

make_coef_table <- function(model, e=NULL, treat){
  temp <- model_stats(model, e)
  stats$coef[treat-1] <- temp[1,1]
  stats$lower[treat-1] <- temp[1,2]
  stats$upper[treat-1] <- temp[1,3]
  return(stats)
}

model_stats <- function(model, e=NULL){
  sum <- fixef(model)
  conf <- confint(model, level=c(0.95), parm="beta_")
  
  stats <- data.frame(sum[1], conf[1,1], conf[1,2])
  if(e==TRUE){
    cat(c(exp(sum[1])), "\n")
    cat(c(exp(conf[1,1])), "\n")
    cat(c(exp(conf[1,2])), "\n")
    stats <- exp(stats)
  } else {
    cat(c(sum[1]), "\n")
    cat(c(conf[1,1]), "\n")
    cat(c(conf[1,2]), "\n")
  }
  clip <- pipe('pbcopy', 'w')
  write.table(stats, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  cat('Coefficients copied to the clipboard')
  return(stats)
  
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


#### Graphing functions ####
## Functions and graphical stuff
# Next I've written some of simple functions that make the graph labels more human-readable. There are also some objects that enable these functions to work properly. *Importantly*, there are different ones for the different experiments, which have been indicated in the comments.

# Lists of the original names (as found in the source data) of the levels of each factor. These are used for the breaks() function in ggplots

# The original genotype names are share across the experiments
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

### EXPERIMENTS 1 & @
## Experiments 1 & 2 had the same bottleneck treatments ('sol' is short for 'solution')
sol_OG_bottleneck = c('2', '3', '4', '5', '6','7', '8', '9')

# The new label objects need a list() for when they apply to facet headers (strip.text), and a vector when they apply to legend labels or axis text
sol_bottleneck_names_facet = list(
  '2'      = expression('10'^-2*''),
  '3'     = expression('10'^-3*''),
  '4'    = expression('10'^-4*''),
  '5'   =  expression('10'^-5*''),
  '6'  = expression('10'^-6*''),
  '7' = expression('10'^-7*''),
  '8' = expression('10'^-8*''),
  '9' = expression('10'^-9*'')
)

sol_bottleneck_names_legend = c(
  expression('10'^-2*''), 
  expression('10'^-3*''), 
  expression('10'^-4*''),
  expression('10'^-5*''),
  expression('10'^-6*''),
  expression('10'^-7*''),
  expression('10'^-8*''),
  expression('10'^-9*'')
)


# The functions which set these up
sol_bottleneck_labeller = function(variable, value) {
  return(sol_bottleneck_names_facet[value])
}


## Finally, a little ggplot object that makes the position of geoms more sensible

pd = position_dodge(0.1)


#### Data ####
data <- read.csv("Phage_bn_exp/original_data/population_counts.csv", header=T)
data <- select(data, -raw, -dilution)
data$ID %<>% as.factor()
data <- plyr::rename(data, c("ID"="replicate"))
#data$bottleneck %<>% as.factor
data$log.pfu <- log(data$pfu+1)
data$log.cfu <- log(data$cfu+1)
#data %<>% na.exclude

#### Model - correlation between cfu and pfu ####

m.null <- lmer(log.pfu~1+(1|replicate), data=data)
m1 <- lmer(log.pfu~log.cfu+(1|replicate), data=data)
m2 <- lmer(log.pfu~log.cfu*bottleneck+(1|replicate), data=data)
m3 <- lmer(log.pfu~log.cfu*timepoint+(1|replicate), data=data)
m.global <- lmer(log.pfu~log.cfu*bottleneck*timepoint+(1|replicate), data=data)

plot(m.null)
plot(m1)
plot(m2)
plot(m3)
plot(m.global)

AIC(m.null, m1, m2, m3, m.global) %>% compare_AICs()

anova(m.null, m1, m2, m3, m.global, test="F")
drop1(m.global, test="Chisq")
sresid <- resid(m.global, type="pearson")
hist(sresid)

summary(lm.global)


model_stats(m.global)
# F-test of the hierarchical model suggests there is a correlation between PFU and CFU

#### Model - cfu covarying with botleneck and timepoint ####
model_stats <- function(model){
  fixed <- fixef(model) %>% 
    as.data.frame() %>% 
    slice(1:8) %>% 
    as.data.frame()
  
  fixed$coef <- rep(0, length(fixed$.))
  fixed$coef[1] = fixed$.[1]
  
  for(i in seq(2, length(fixed$coef))){
    fixed$coef[i] = fixed$coef[1]+fixed$.[i]
  }
  
  CI <- confint(model, parm="beta_")
  
  CIdf <- as.data.frame(CI) %>% 
    slice(1:8) %>% 
    as.data.frame()
  
  CIdf$lowerCI <- rep(0, length(CIdf$`2.5 %`))
  CIdf$upperCI <- rep(0, length(CIdf$`2.5 %`))
  CIdf$lowerCI[1] <- CIdf$`2.5 %`[1]; CIdf$upperCI[1] <- CIdf$`97.5 %`[1]
  
  for(i in seq(2, length(CIdf$`2.5 %`))){
    CIdf$lowerCI[i] = CIdf$`2.5 %`[1]+CIdf$`2.5 %`[i]
    CIdf$upperCI[i] = CIdf$`97.5 %`[1]+CIdf$`97.5 %`[i]
  }
  
  final <- bind_cols(fixed, CIdf) %>% 
    select(coef, lowerCI, upperCI) %>% 
    exp
  
  clip <- pipe('pbcopy', 'w')
  write.table(final, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  cat('Exponentiated fixed effect coefficients & 95% confidence intervals of bottlenecks 2-9 at specified timepoint copied to the clipboard\n')
}

m.null <- lmer(log.cfu~1+(1|replicate), data=data)
m1 <- lmer(log.cfu~timepoint+(1|replicate), data=data)
m2 <- lmer(log.cfu~bottleneck+(1|replicate), data=data)
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)

plot(m.null)
plot(m1)
plot(m2)
plot(m.global)

AIC(m.null, m1, m2, m.global) %>% compare_AICs()

anova(m.null, m1, m2, m.global, test="Chisq")
drop1(m.global, test="Chisq")
summary(m.global)
tidy(m.global) %>% head(16)

data$timepoint %<>% relevel(ref="t5")
data$timepoint %<>% relevel(ref="t4")
data$timepoint %<>% relevel(ref="t3")
data$timepoint %<>% relevel(ref="t2")
data$timepoint %<>% relevel(ref="t1")
data$timepoint %<>% relevel(ref="t0")

# Copy coefficients ####

data$bottleneck %<>% relevel(ref="2")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=2)

data$bottleneck %<>% relevel(ref="3")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=3)

data$bottleneck %<>% relevel(ref="4")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=4)

data$bottleneck %<>% relevel(ref="5")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=5)

data$bottleneck %<>% relevel(ref="6")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=6)

data$bottleneck %<>% relevel(ref="7")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=7)

data$bottleneck %<>% relevel(ref="8")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=8)

data$bottleneck %<>% relevel(ref="9")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=9)

data$bottleneck %<>% relevel(ref="p2")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=10)

data$bottleneck %<>% relevel(ref="p3")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=11)

data$bottleneck %<>% relevel(ref="p4")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=12)

data$bottleneck %<>% relevel(ref="p5")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=13)

data$bottleneck %<>% relevel(ref="p6")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=14)

data$bottleneck %<>% relevel(ref="p7")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=15)

data$bottleneck %<>% relevel(ref="p8")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=16)

data$bottleneck %<>% relevel(ref="p9")
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)
stats <- make_coef_table(m.global, e=T, treat=17)

clip <- pipe('pbcopy', 'w')
write.table(stats, file=clip, sep='\t', row.names = F, col.names = F)
close(clip)
cat('Coefficients copied to the clipboard')


#### Raw fig - just cfu ####
ID <- c(seq(2.1,2.6,.1), seq(3.1,3.6,.1), seq(4.1,4.6,.1),
        seq(5.1,5.6,.1), seq(6.1,6.6,.1), seq(7.1,7.6,.1),
        seq(8.1,8.6,.1), seq(9.1,9.6,.1)) %>% 
  rep(6)
data$ID <- as.factor(ID)

raw_plot <- ggplot(aes(y=cfu+1, x=timepoint, group=ID), 
                   data=data)+
  
  geom_path(stat='identity', 
            aes(colour=bottleneck),
            position=pd)+
  geom_point(stat='identity', 
             aes(fill=bottleneck, colour=bottleneck),
             position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("C.f.u. ml"*{}^{-1}*"")))+
  ggtitle('')+
  #facet_wrap(~bottleneck)+
  
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
  theme(legend.text = element_text(size=12))

quartz()
raw_plot

#### Summary figure - data formatting ####
data <- read.csv("./Phage_bn_exp/data/counts/counts_summary_nbinom.csv")

dataM <- melt(data, measure.vars = c("pfu", "cfu"))
dataM <- plyr::rename(dataM, c("variable"="measurement"))
lower1 <- slice(dataM, 1:48) %>% 
  select(pfu.lower) %>% 
  plyr::rename(c("pfu.lower"="lower"))
lower2 <- slice(dataM, 49:96) %>% 
  select(cfu.lower) %>% 
  plyr::rename(c("cfu.lower"="lower"))
upper1 <- slice(dataM, 1:48) %>% 
  select(pfu.upper) %>% 
  plyr::rename(c("pfu.upper"="upper"))
upper2 <- slice(dataM, 49:96) %>% 
  select(cfu.upper) %>% 
  plyr::rename(c("cfu.upper"="upper"))
dataM$lower <- bind_rows(lower1, lower2)
dataM$upper <- bind_rows(upper1, upper2)

rm(lower1, lower2, upper1, upper2)

dataM <- select(dataM, bottleneck, timepoint, measurement, value, lower, upper)
dataM$bottleneck %<>% as.factor()

phageIDs <- c("p1", "p2", "p3",
              "p4", "p5", "p6") %>% 
  rep( (length(dataM$bottleneck)/2)/6 )
hostIDs <- c("h1", "h2", "h3",
             "h4", "h5", "h6") %>% 
  rep( (length(dataM$bottleenck)/2)/6 )

dataM$ID <- as.factor(c(phageIDs, hostIDs))

#### Summary figure - make the bastard ####
raw_phage_plot <- ggplot(aes(y=pfu+1, x=timepoint, group=bottleneck), 
                   data=data)+
  
  #geom_path(stat='identity', 
    #       aes(colour=bottleneck, linetype=measurement),
   #         position=pd)+
   # scale_linetype_manual(name='Measurement',
   #                     values=c(1, 2),
    #                    breaks = c("pfu", "cfu"),
    #                    labels = c("Phage", "Host"))+
  geom_point(stat='identity', 
             aes(shape=measurement, fill=bottleneck),
             position=pd)+
  #scale_shape_manual(name='Measurement',
  #                   values = c(21,24),
   #                  breaks = c("pfu", "cfu"),
    #                 labels = c("Phage", "Host"))+
  geom_errorbar(stat="identity", aes(ymin=lower+1, ymax=upper+1, linetype=measurement, colour=bottleneck),
                position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("P.f.u. ml"*{}^{-1}*"/ C.f.u. ml"*{}^{-1}*"")))+
  ggtitle('')+
  #facet_grid(~bottleneck)+
  
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
  annotate("text", 1.5, 1e+2, vjust=-1, label="Phage detection limit", colour="red")
quartz()
raw_phage_plot


#### Summary fig - just pfu ####
data$bottleneck %<>% as.factor()
data <- plyr::rename(data, c("bottleneck"="Bottleneck"))

phage_plot <- ggplot(aes(y=pfu+1, x=timepoint, group=Bottleneck), 
                   data=data)+
  
  geom_path(stat='identity', 
          aes(colour=Bottleneck),
          position=pd)+
  geom_point(stat='identity', 
             aes(fill=Bottleneck, colour=Bottleneck),
             position=pd)+
  geom_errorbar(stat="identity", aes(ymin=pfu.lower+1, ymax=pfu.upper+1, 
                                     colour=Bottleneck),
                position=pd, width=0)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("P.f.u. ml"*{}^{-1}*"")))+
  ggtitle('')+
  #facet_wrap(~Bottleneck)+
  
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
  annotate("text", 1.5, 1e+2, vjust=-1, label="Phage detection limit", colour="red")
quartz()
phage_plot

host_plot <- ggplot(aes(y=cfu+1, x=timepoint, group=Bottleneck), 
                     data=data)+
  
  geom_path(stat='identity', 
            aes(colour=Bottleneck),
            position=pd,linetype=2)+
  geom_point(stat='identity', 
             aes(fill=Bottleneck, colour=Bottleneck),
             position=pd, shape=24)+
  geom_errorbar(stat="identity", aes(ymin=cfu.lower+1, ymax=cfu.upper+1, 
                                     colour=Bottleneck),
                position=pd, width=0)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("C.f.u. ml"*{}^{-1}*"")))+
  ggtitle('')+
  #facet_wrap(~Bottleneck)+
  
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
  theme(legend.text = element_text(size=12))
quartz()
host_plot
