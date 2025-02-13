### Full bottleneck experiment - population analysis ####
# Created: 18/6/18 by Jack Common

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

#### Stats functions ####
# Extracts the intercept coefficient (mean) and 95% CIs from the GLM objects, with added functionality to copy those to the clipboard for easier input to summary dataframes

# Enable this if using Ubuntu
#clipboard <- function(x, sep="\t", row.names=FALSE, col.names=FALSE){
#    con <- pipe("xclip -selection clipboard -i", open="w")
#    write.table(x, con, sep=sep, row.names=row.names, col.names=col.names)
#    close(con)
#}

stats <- data.frame(coef=rep(0,8), lower=rep(0,8), upper=rep(0,8))

make_coef_table <- function(model, e=NULL, treat){
  temp <- model_stats(model, e)
  stats$coef[treat-1] <- temp[1,1]
  stats$lower[treat-1] <- temp[1,2]
  stats$upper[treat-1] <- temp[1,3]
  return(stats)
}

model_stats <- function(model, e=NULL){
  sum <- coef(m.global)
  conf <- confint(m.global, level=c(0.95))

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
data <- read.csv("Full_bn_2/original_data/population_counts.csv", header=T)
data <- select(data, -raw, -dilution)
data$ID %<>% as.factor()
data <- plyr::rename(data, c("ID"="replicate"))
data$bottleneck %<>% as.factor
data$log.pfu <- log(data$pfu+1)
data$log.cfu <- log(data$cfu+1)
#data %<>% na.exclude


#### Model - correlation between cfu and pfu ####

m.null <- lmer(log.pfu~1+(1|replicate), data=data)
m1 <- lmer(log.pfu~bottleneck+(1|replicate), data=data)
m2 <- lmer(log.pfu~timepoint+(1|replicate), data=data)
m3 <- lmer(log.pfu~bottleneck*timepoint+(1|replicate), data=data)
m4 <- lmer(log.pfu~log.cfu+(1|replicate), data=data)
m5 <- lmer(log.pfu~log.cfu*bottleneck+(1|replicate), data=data)
m6 <- lmer(log.pfu~log.cfu*timepoint+(1|replicate), data=data)
m.global <- lmer(log.pfu~log.cfu*bottleneck*timepoint+(1|replicate), data=data)

plot(m.null)
plot(m1)
plot(m2)
plot(m3)
plot(m4)
plot(m5)
plot(m6)
plot(m.global)

AIC(m.null, m1, m2, m3, 
    m4, m5, m6, m.global) %>% compare_AICs()

anova(m.null, m1, m2, m3, 
      m4, m5, m6, m.global, test="LRT")
drop1(m.global, test="Chisq")

summary(m.global)
model_stats(m4)
# F-test of the hierarchical model suggests there isn't a correlation between PFU and CFU

#### Model - cfu covarying with bottleneck and timepoint ####
m.null <- lmer(log.cfu~1+(1|replicate), data=data)
m1 <- lmer(log.cfu~timepoint+(1|replicate), data=data)
m2 <- lmer(log.cfu~bottleneck+(1|replicate), data=data)
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data)

#par(mfrow=c(2,2))
plot(m.null)
plot(m1)
plot(m2)
plot(m.global)

AIC(m.null, m1, m2, m.global) %>% compare_AICs()

anova(m.null, m1, m2, m.global, test="Chisq")
drop1(m.global, test="Chisq")
sresid <- resid(m.global, type = "pearson")  # Extract the standardised residuals
hist(sresid)
r.squaredGLMM(m.global) # all of the conditional R2 is explained by the marginal R2, so the random term is irrelevant. Hence,
                        # a linear model should be adequate
m.global <- lm(log.cfu~bottleneck*timepoint, data=data)

summary(m.global)
plot(m.global)
model_stats_linear(m.global)

data$treatment %<>% relevel(ref="phageN")
data$treatment %<>% relevel(ref="crisprN")
data$treatment %<>% relevel(ref="exp")

data$timepoint %<>% relevel(ref="t5")
data$timepoint %<>% relevel(ref="t4")
data$timepoint %<>% relevel(ref="t3")
data$timepoint %<>% relevel(ref="t2")
data$timepoint %<>% relevel(ref="t1")
data$timepoint %<>% relevel(ref="t0")

# Copy coefficients ####

data$bottleneck %<>% relevel(ref="2")
m.global <- lm(log.cfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=2)

data$bottleneck %<>% relevel(ref="3")
m.global <- lm(log.cfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=3)

data$bottleneck %<>% relevel(ref="4")
m.global <- lm(log.cfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=4)

data$bottleneck %<>% relevel(ref="5")
m.global <- lm(log.cfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=5)

data$bottleneck %<>% relevel(ref="6")
m.global <- lm(log.cfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=6)

data$bottleneck %<>% relevel(ref="7")
m.global <- lm(log.cfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=7)

data$bottleneck %<>% relevel(ref="8")
m.global <- lm(log.cfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=8)

data$bottleneck %<>% relevel(ref="9")
m.global <- lm(log.cfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=9)

clip <- pipe('pbcopy', 'w')
write.table(stats, file=clip, sep='\t', row.names = F, col.names = F)
close(clip)
cat('Coefficients copied to the clipboard')


#### Model - CFU*PFU for control data ####
data_cont <- filter(data, bottleneck%in%c("c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9"))

m.null <- lmer(log.pfu~1+(1|replicate), data=data_cont)
m1 <- lmer(log.pfu~bottleneck+(1|replicate), data=data_cont)
m2 <- lmer(log.pfu~timepoint+(1|replicate), data=data_cont)
m3 <- lmer(log.pfu~bottleneck*timepoint+(1|replicate), data=data_cont)
m4 <- lmer(log.pfu~log.cfu+(1|replicate), data=data_cont)
m5 <- lmer(log.pfu~log.cfu*bottleneck+(1|replicate), data=data_cont)
m6 <- lmer(log.pfu~log.cfu*timepoint+(1|replicate), data=data_cont)
m.global <- lmer(log.pfu~log.cfu*bottleneck*timepoint+(1|replicate), data=data_cont)

plot(m.null)
plot(m1)
plot(m2)
plot(m3)
plot(m4)
plot(m5)
plot(m6)
plot(m.global)

AIC(m.null, m1, m2, m3, m4, m5, m6, m.global) %>% compare_AICs()

anova(m.null, m1, m2, m3, m4, m5, m6, m.global, test="F")
drop1(m.global, test="Chisq")

# PFU covaries with CFU in the control data, but not in the experimental data

#### Model - CFU*PFU for experiment data ####
data_exp <-  filter(data, bottleneck%in%c(seq(2,9,1)))

m.null <- lmer(log.pfu~1+(1|replicate), data=data_exp)
m1 <- lmer(log.pfu~bottleneck+(1|replicate), data=data_exp)
m2 <- lmer(log.pfu~timepoint+(1|replicate), data=data_exp)
m3 <- lmer(log.pfu~bottleneck*timepoint+(1|replicate), data=data_exp)
m4 <- lmer(log.pfu~log.cfu+(1|replicate), data=data_exp)
m5 <- lmer(log.pfu~log.cfu*bottleneck+(1|replicate), data=data_exp)
m6 <- lmer(log.pfu~log.cfu*timepoint+(1|replicate), data=data_exp)
m.global <- lmer(log.pfu~log.cfu*bottleneck*timepoint+(1|replicate), data=data_exp)

plot(m.null)
plot(m1)
plot(m2)
plot(m3)
plot(m4)
plot(m5)
plot(m6)
plot(m.global)

AIC(m.null, m1, m2, m3, m4, m5, m6, m.global) %>% compare_AICs()

sresid <- resid(m.global, type = "pearson")  # Extract the standardised residuals
hist(sresid)

par(mfrow=c(1,2))
plot(sresid ~ data_exp$bottleneck*data_exp$timepoint) 

anova(m.null, m1, m2, m3, m4, m5, m6, m.global, test="F")
drop1(m.global, test="Chisq")
summary(m.global)
r.squaredGLMM(m.global)

#### Model - cfu covarying with bottleneck and timepoint EXPERIMENT DATA ####
m.null <- lmer(log.cfu~1+(1|replicate), data=data_exp)
m1 <- lmer(log.cfu~timepoint+(1|replicate), data=data_exp)
m2 <- lmer(log.cfu~bottleneck+(1|replicate), data=data_exp)
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data_exp)

#par(mfrow=c(2,2))
plot(m.null)
plot(m1)
plot(m2)
plot(m.global)

AIC(m.null, m1, m2, m.global) %>% compare_AICs()

anova(m.null, m1, m2, m.global, test="Chisq")
drop1(m.global, test="Chisq")
sresid <- resid(m.global, type = "pearson")  # Extract the standardised residuals
hist(sresid)
r.squaredGLMM(m.global) # all of the conditional R2 (random) is explained by the marginal R2 (fixed), so the random term is irrelevant. Hence,
                        # a linear model should be adequate

m.global <- lm(log.cfu~bottleneck*timepoint, data=data_exp)
summary(m.global)

#### Model - cfu covarying with bottleneck and timepoint CRISPR -VE DATA ####
m.null <- lmer(log.cfu~1+(1|replicate), data=data_cont)
m1 <- lmer(log.cfu~timepoint+(1|replicate), data=data_cont)
m2 <- lmer(log.cfu~bottleneck+(1|replicate), data=data_cont)
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|replicate), data=data_cont)

#par(mfrow=c(2,2))
plot(m.null)
plot(m1)
plot(m2)
plot(m.global)

AIC(m.null, m1, m2, m.global) %>% compare_AICs()

anova(m.null, m1, m2, m.global, test="Chisq")
drop1(m.global, test="Chisq")
sresid <- resid(m.global, type = "pearson")  # Extract the standardised residuals
hist(sresid)
r.squaredGLMM(m.global) # all of the conditional R2 (random) is explained by the marginal R2 (fixed), so the random term is irrelevant. Hence,
# a linear model should be adequate

m.global <- lm(log.cfu~bottleneck*timepoint, data=data_exp)
summary(m.global)

#### Model - pfu covarying with bottleneck & timepoint ####
m1 <- lm(log.pfu~1, data=data)
m2 <- lm(log.pfu~timepoint, data=data)
m3 <- lm(log.pfu~bottleneck, data=data)
m4 <- lm(log.pfu~treatment, data=data)
m5 <- lm(log.pfu~bottleneck*timepoint, data=data)
m6 <- lm(log.pfu~bottleneck*treatment, data=data)
m7 <- lm(log.pfu~timepoint*treatment, data=data)
m.global <- lm(log.pfu~bottleneck*timepoint*treatment, data=data)

par(mfrow=c(2,2))
plot(m1)
plot(m2)
plot(m3)
plot(m4)
plot(m5)
plot(m6)
plot(m7)
plot(m.global)

AIC(m1, m2, m3, m4,
    m5, m6, m7, m.global) %>% compare_AICs()

anova(m1, m2, m3, m4,
      m5, m6, m7, m.global, test="Chisq")
summary(m.global)

data$treatment %<>% relevel(ref="phageN")
data$treatment %<>% relevel(ref="crisprN")
data$treatment %<>% relevel(ref="exp")

data$timepoint %<>% relevel(ref="t5")
data$timepoint %<>% relevel(ref="t4")
data$timepoint %<>% relevel(ref="t3")
data$timepoint %<>% relevel(ref="t2")
data$timepoint %<>% relevel(ref="t1")
data$timepoint %<>% relevel(ref="t0")

# Copy coefficients ####

data$bottleneck %<>% relevel(ref="2")
m.global <- lm(log.pfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=2)

data$bottleneck %<>% relevel(ref="3")
m.global <- lm(log.pfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=3)

data$bottleneck %<>% relevel(ref="4")
m.global <- lm(log.pfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=4)

data$bottleneck %<>% relevel(ref="5")
m.global <- lm(log.pfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=5)

data$bottleneck %<>% relevel(ref="6")
m.global <- lm(log.pfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=6)

data$bottleneck %<>% relevel(ref="7")
m.global <- lm(log.pfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=7)

data$bottleneck %<>% relevel(ref="8")
m.global <- lm(log.pfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=8)

data$bottleneck %<>% relevel(ref="9")
m.global <- lm(log.pfu~bottleneck*timepoint*treatment, data=data)
stats <- make_coef_table(m.global, e=T, treat=9)

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
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("/ C.f.u. ml"*{}^{-1}*"")))+
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
data <- read.csv("./Full_bn_2/summary_data/population_summary.csv")

dataM <- melt(data, measure.vars = c("pfu", "cfu"))
dataM <- plyr::rename(dataM, c("variable"="measurement"))
lower1 <- slice(dataM, 1:144) %>% 
  select(pfu.lower) %>% 
  plyr::rename(c("pfu.lower"="lower"))
lower2 <- slice(dataM, 145:288) %>% 
  select(cfu.lower) %>% 
  plyr::rename(c("cfu.lower"="lower"))
upper1 <- slice(dataM, 1:144) %>% 
  select(pfu.upper) %>% 
  plyr::rename(c("pfu.upper"="upper"))
upper2 <- slice(dataM, 145:288) %>% 
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
raw_plot <- ggplot(aes(y=value+1, x=timepoint, group=bottleneck), 
                   data=dataM)+
  
 # geom_path(stat='identity', 
  #         aes(colour=bottleneck, linetype=measurement),
 #           position=pd)+
    #scale_linetype_manual(name='Measurement',
     #                   values=c(1, 2),
      #                  breaks = c("pfu", "cfu"),
       #                 labels = c("Phage", "Host"))+
  geom_point(stat='identity', 
             aes(shape=measurement, fill=bottleneck),
             position=pd)+
  scale_shape_manual(name='Measurement',
                     values = c(21,24),
                     breaks = c("pfu", "cfu"),
                     labels = c("Phage", "Host"))+
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
raw_plot


#### Summary fig - just pfu ####
data$bottleneck %<>% as.factor()
data <- plyr::rename(data, c("bottleneck"="Bottleneck"))

phage_plot <- ggplot(aes(y=pfu, x=timepoint, group=treatment), 
                   data=data)+
  
  geom_path(stat='identity', 
          aes(linetype=treatment),
          position=pd)+
  geom_point(stat='identity', 
             aes(shape=treatment),
             position=pd)+
  geom_errorbar(stat="identity", aes(ymin=pfu.lower+1, ymax=pfu.upper+1),
                position=pd, width=0)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("P.f.u. ml"*{}^{-1}*"")))+
  ggtitle('')+
  facet_wrap(~Bottleneck)+
  
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

  #geom_hline(yintercept = 1e+2, linetype=2, colour="red")+
  #annotate("text", 1.5, 1e+2, vjust=-1, label="Phage detection limit", colour="red")+
  NULL
quartz()
phage_plot

ggsave("phage_summary.png", phage_plot, path="./Revisions/",
       device="png", dpi=300, width=27, height=20, unit=c("cm"))

host_plot <- ggplot(aes(y=cfu+1, x=timepoint, group=treatment), 
                     data=data)+
  
  geom_path(stat='identity', 
            aes(colour=treatment),
            position=pd)+
  geom_point(stat='identity', 
             aes(colour=treatment),
             position=pd)+
  geom_errorbar(stat="identity", aes(ymin=cfu.lower+1, ymax=cfu.upper+1, colour=treatment),
                position=pd, width=0)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("C.f.u. ml"*{}^{-1}*"")))+
  ggtitle('')+
  facet_wrap(~Bottleneck)+
  
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
  NULL
quartz()
host_plot

ggsave("host_summary.png", host_plot, path="./Revisions/",
       device="png", dpi=300, width=27, height=20, unit=c("cm"))
