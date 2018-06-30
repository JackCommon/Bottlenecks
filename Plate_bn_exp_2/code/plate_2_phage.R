### Plate Bottleneck Experiment 2 Analysis
# Created: 15/4/18 by Jack Common

rm(list=ls())

#### Dependencies ####
library(ggplot2)
library(scales)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(broom)
library(magrittr)
library(cowplot)
library(lme4)
library(MuMIn)

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

# These vectors just give tidy labels for the figures for experiment 3
oneclone_IDs <- c('M.1', 'M.2', 'M.3', 'M.4', 'M.5', 'M.6')
fiveclone_IDs <- c('5.1', '5.2', '5.3', '5.4', '5.5', '5.6')
fiftyclone_IDs <- c('50.1', '50.2', '50.3', '50.4', '50.5', '50.6')

plate_replicate_names_legend = c('1', '2', '3', '4', '5', '6')

plate_OG_bottleneck = c('1-clone', '5-clone', '50-clone')

plate_bottleneck_names_facet = list(
  '1-clone'      = expression('1-clone'),
  '5-clone'          = expression('5-clone'),
  '50-clone'         = expression('50-clone')
)

plate_bottleneck_names_legend = c(
  expression('1-clone'), 
  expression('5-clone'), 
  expression('50-clone')
)

# Function for experiment 3 bottleneck labelling
plate_bottleneck_labeller = function(variable, value) {
  return(plate_bottleneck_names_facet[value])
}

pd = position_dodge(0.1)

#### Load and format data ####
phage <- read.csv("./Plate_bn_exp_2/original_data/plate2_counts.csv", header = T)
phage <- select(phage, -raw.count, -dilution)
phage$ID %<>% as.factor()
#phage %<>% na.exclude
#$log.pfu <- log10(phage$pfu+1)

## Melt data for raw value plots
phageM <- melt(phage, measure.vars = c("pfu", "cfu"))
phageM <- plyr::rename(phageM, c("variable"="measurement"))

phageIDs <- c("p1", "p2", "p3", "p4",
              "p5", "p6", "p7", "p8",
              "p9", "p10", "p11", "p12") %>% 
  rep( (length(phageM$ID)/2)/12 )
hostIDs <- c("h1", "h2", "h3", "h4",
             "h5", "h6", "h7", "h8",
             "h9", "h10", "h11", "h12") %>% 
  rep( (length(phageM$ID)/2)/12 )

ID2 <- c(phageIDs, hostIDs)

phageM$ID2 <- as.factor(ID2)

## 1-clone plot ####
mono_plot <- ggplot(aes(y=value+1, x=timepoint, group=ID2), 
                   data=subset(phageM, bottleneck == '1-clone'))+
  
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
  ggtitle('1-clone')+
  
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
  annotate("text", 1.5, 1e+2, vjust=-1, label="Detection limit", colour="red")

quartz()
mono_plot

## 5-clone plot ####
fiveclone_plot <- ggplot(aes(y=value+1, x=timepoint, group=ID2), 
                                         data=subset(phageM, bottleneck == '5-clone'))+
  
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
  ggtitle('5-clone')+
  
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
  annotate("text", 1.5, 1e+2, vjust=-1, label="Detection limit", colour="red")


fiveclone_plot
## 50-clone plot ####
fiftyclone_plot <- ggplot(aes(y=value+1, x=timepoint, group=ID2), 
                                          data=subset(phageM, bottleneck == '50-clone'))+
  
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
  ggtitle('50-clone')+
  
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
  annotate("text", 1.5, 1e+2, vjust=-1, label="Detection limit", colour="red")


fiftyclone_plot

## Arrange and save raw plots ####
mono_plot <- mono_plot + theme(plot.margin = unit(c(2,2,0,1), 'pt'))
fiveclone_plot <- fiveclone_plot + theme(plot.margin = unit(c(2,2,0,1), 'pt'))
fiftyclone_plot <- fiftyclone_plot + theme(plot.margin = unit(c(0,2,2,1), 'pt'))

raw_fig <- plot_grid(mono_plot+labs(x='')+theme(legend.position = 'none'), 
                   fiveclone_plot+labs(y='', x='')+theme(legend.position = 'none'), 
                   fiftyclone_plot+theme(legend.position = 'none'))
quartz()
raw_fig

detach("package:cowplot")

ggsave('all_raw.png', raw_fig, device = 'png',
       path = './Plate_bn_exp_2/figs/', width=27, height=17, unit=c('cm'), dpi=300)


#### Phage survival analysis ####
library(survival)
library(rms)
library(car)
library(multcomp)
library(relaimpo)

phage<-read.csv("./Plate_bn_exp_2/summary_data/survival_data.csv", header=T)
attach(phage)
names(phage)

summary(KM<-survfit(Surv(time_to_death,status)~1))
plot(KM, ylab="Survivorship", xlab="Transfer")

# KM ~ group
summary(KM<-survfit(Surv(time_to_death,status)~bottleneck))

jpeg("./figs/survplot.jpg", width=20, height=15, units="in", res=300)
par(mfrow=c(1,1), xpd=TRUE, oma=c(1.5,2.5,1,1), mai=c(1,1,1,1.2), bty="l", pty="s")

plot(survfit(Surv(phage$time_to_death,phage$status)~bottleneck), lty=c(1,3,5), lwd=c(5,5,5), ylab="", xlab="", axes=FALSE, ylim=c(0,1), xlim=c(0,5))

axis(1, tcl=-0.1, pos=0, cex.axis=1, lwd=c(3), cex.axis=2)
axis(1, at=2.5, lab="Days post-infection (d.p.i.)", tcl=0, line=2, cex.axis=3)

axis(2, tcl=-0.1, pos=-0, cex.axis=1, las=2, lwd=c(3), cex.axis = 2)
axis(2, at=0.5, lab="Proportion of phage\npopulations surviving", line=4, cex.axis=3, tcl=0)

legend(0.8,0.5, title=c("Bottleneck"),
       legend=c("1-clone", "5-clone", "50-clone"), 
       bty="o", lty=c(1,3,5), lwd=c(5,5,5), cex=3, adj=0)
dev.off()

# Cox proportional hazards model
model3<-coxph(Surv(time_to_death,status)~bottleneck)
summary(model3)

model3$loglik

anova(model3)
tapply(predict(model3),bottleneck,mean)

exp1.tukey <- summary(glht(model3, linfct = mcp(bottleneck = "Tukey")))
class(exp(exp1.tukey$test$coefficients))
exp1.tukey
exp1.tukey$test$coefficients
exp1.tukey$test$tstat
exp1.tukey$test$pvalues

HRs <- exp(exp1.tukey$test$coefficients)
SEs <- exp(exp1.tukey$test$sigma)
Z <- exp1.tukey$test$tstat
P <- exp1.tukey$test$pvalues

exp1.HRs <- data.frame(HRs, SEs, Z, P)
clip = pipe('pbcopy', 'w')
write.table(exp1.HRs, file=clip, sep='\t', row.names = F, col.names = F)
close(clip)

plot(survfit(model3), lty=c(1,2,3))


#### Model - PFU only ####
phage$log.pfu <- log(phage$pfu+1)

m.null <- lmer(log.pfu~1+(1|ID), data=phage)
m1 <- lmer(log.pfu~timepoint+(1|ID), data=phage)
m2 <- lmer(log.pfu~bottleneck+(1|ID), data=phage)
m.global <- lmer(log.pfu~bottleneck*timepoint+(1|ID), data=phage)
#m4 <- lmer(log.pfu~bottleneck*timepoint+(1|timepoint), data=phage)
#m5 <- lmer(log.pfu~bottleneck*timepoint+(1|bottleneck), data=phage)
#m6 <- lmer(log.pfu~bottleneck*timepoint+(ID|timepoint), data=phage)
#m7 <- lmer(log.pfu~bottleneck*timepoint+(ID|bottleneck), data=phage)
#m8 <- lmer(log.pfu~bottleneck*timepoint+(timepoint|bottleneck), data=phage)
#m.global <- lmer(log.pfu~bottleneck*timepoint+(timepoint|bottleneck)+(1|ID), data=phage)

plot(m.null)
plot(m1)
plot(m2)
#plot(m3)
#plot(m4)
#plot(m5)
#plot(m6)
#plot(m7)
#plot(m8)
plot(m.global)
AIC(m.null, m1, m2, m.global) %>% compare_AICs()
# Model 3 (bottleneck*timepoint as fixed with replicate as random) reduces AIC the most
#m.global <- lmer(log.pfu~bottleneck*timepoint+(1|ID), data=phage)
summary(m.global)
#anova(m.null, m1, m2, m.global, test="Chisq")
drop1(m.global, test="Chisq")
sresid <- resid(m.global, type = "pearson")  # Extract the standardised residuals
hist(sresid)
R2 <- r.squaredGLMM(m.global)
R2[1]/R2[2]*100

m.global <- lm(log.pfu~bottleneck*timepoint, data=phage)
tidy(m.global)

#### Model - PFU*CFU ####
phage$log.cfu <- log(phage$cfu+1)

m.null <- lmer(log.pfu~1+(1|ID), data=phage)
m1 <- lmer(log.pfu~log.cfu+(1|ID), data=phage)
m2 <- lmer(log.pfu~log.cfu*bottleneck+(1|ID), data=phage)
m3 <- lmer(log.pfu~log.cfu*timepoint+(1|ID), data=phage)
m.global <- lmer(log.pfu~log.cfu*bottleneck*timepoint+(1|ID), data=phage)

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

summary(m.global)
R2 <- r.squaredGLMM(m.global)
R2[1]/R2[2]*100

### Model - CFU only ####

m.null <- lmer(log.cfu~1+(1|ID), data=phage)
m1 <- lmer(log.cfu~bottleneck+(1|ID), data=phage)
m2 <- lmer(log.cfu~timepoint+(1|ID), data=phage)
m.global <- lmer(log.cfu~bottleneck*timepoint+(1|ID), data=phage)

plot(m.null)
plot(m1)
plot(m2)
plot(m.global)

AIC(m.null, m1, m2, m.global) %>% compare_AICs()

anova(m.null, m1, m2, m.global, test="F")
drop1(m.global, test="Chisq")

sresid <- resid(m.global, type="pearson")
hist(sresid)

summary(m.global)
R2 <- r.squaredGLMM(m.global)
R2[1]/R2[2]*100

m.global <- lm(log.cfu~bottleneck*timepoint, data=phage)
summary(m.global)
