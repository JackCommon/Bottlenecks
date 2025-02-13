#### Survival analysis ####
## Created 4/6/18 by Jack Common

rm(list=ls(all=TRUE))

#### Dependencies ####
library(survival)
library(rms)
library(car)
library(multcomp)
library(relaimpo)
library(dplyr)
library(magrittr)

#### Functions ####
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

# Experiment 1 ####
phage<-read.csv("./Full_bn_2/original_data/survival.csv", header=T)
phage$bottleneck %<>% as.factor()
phage$replicate %<>% as.factor
attach(phage)
names(phage)

summary(KM<-survfit(Surv(time_to_death,status)~1))
plot(KM, ylab="Survivorship", xlab="Transfer")

# Cox proportional hazards model
m.null <- coxph(Surv(time_to_death, status)~1)
m1 <- coxph(Surv(time_to_death, status)~bottleneck)
m2 <- coxph(Surv(time_to_death, status)~treatment)

AIC(m1, m2) %>% compare_AICs()
anova(m1, m2, test="Chisq")
drop(m1)

sresid <- resid(m1, type="schoenfeld")
hist(sresid)
summary(m.global)

tapply(predict(m.global),bottleneck,mean)

summary(glht(m1, linfct = mcp(bottleneck = "Tukey")))

exp1.tukey <- summary(glht(m1, linfct = mcp(bottleneck = "Tukey")))
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
exp1.HRs
clip = pipe('pbcopy', 'w')
write.table(exp1.HRs, file=clip, sep='\t', row.names = F, col.names = F)
close(clip)

plot(survfit(m.global), xlim=c(0,5))

# Experiment 2
phage<-read.csv("./Phage_bn_exp/original_data/survival.csv", header=T)
phage$bottleneck %<>% as.factor()
phage %<>% filter(bottleneck!="9")
attach(phage)
names(phage)

# Cox proportional hazards model
model3<-coxph(Surv(time_to_death,status)~bottleneck, na.action = na.exclude)
summary(model3)

model3$loglik

anova(model3)
tapply(predict(model3),bottleneck,mean)

exp1.tukey <- summary(glht(model3, linfct = mcp(bottleneck = "Tukey")))
class(exp(exp2.tukey$test$coefficients))
exp2.tukey
exp2.tukey$test$coefficients
exp2.tukey$test$tstat
exp2.tukey$test$pvalues

HRs <- exp(exp2.tukey$test$coefficients)
SEs <- exp(exp2.tukey$test$sigma)
Z <- exp2.tukey$test$tstat
P <- exp2.tukey$test$pvalues

exp2.HRs <- data.frame(HRs, SEs, Z, P)
exp2.HRs
clip = pipe('pbcopy', 'w')
write.table(exp2.HRs, file=clip, sep='\t', row.names = F, col.names = F)
close(clip)

plot(survfit(model3), lty=c(1,2,3),xlim=c(0,5))


# Experiment 3
phage<-read.csv("./Plate_bn_exp/data/exp3.phage.surv.csv", header=T)
phage$bottleneck %<>% relevel(ref="monoculture")
attach(phage)
names(phage)

summary(KM<-survfit(Surv(time_to_death,status)~1))
plot(KM, ylab="Survivorship", xlab="Transfer")

# KM ~ group
summary(KM<-survfit(Surv(time_to_death,status)~bottleneck))

par(mfrow=c(1,1), xpd=TRUE, oma=c(1.5,2.5,1,1), mai=c(1,1,1,1.2), bty="l", pty="s")

plot(survfit(Surv(phage$time_to_death,phage$status)~bottleneck), lty=c(1,3,5), lwd=c(1.5,1.5,1.5), ylab="", xlab="", axes=FALSE, ylim=c(0,1), xlim=c(0,5))

axis(1, tcl=-0.1, pos=0, cex.axis=1)
axis(1, at=2.5, lab="Days post-infection (d.p.i.)", tcl=0, line=2, cex.axis=1.1)

axis(2, tcl=-0.1, pos=-0, cex.axis=1, las=2)
axis(2, at=0.5, lab="Proportion of phage\npopulations surviving", line=2, cex.axis=1.2, tcl=0)

legend(5.1,1, legend=c("1-clone", "5-clone", "50-clone"), bty="o", lty=c(1,3,5), lwd=c(1.5,1.5,1.5), cex=1.1, adj=0)


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
