library(dplyr)
library(magrittr)

#setwd("~/OneDrive/Data/Bottlenecks/")
setwd("~/Documents/OneDrive - University of Exeter/Data/Bottlenecks/")

phage <- read.csv("Phage_bn_exp/data/counts/counts_master_csv.csv")
phage %<>% select(-cfu)
phage %<>% na.exclude()
phage$bottleneck %<>% as.factor
phage$ID %<>% as.factor

library(MASS)
detach("package:MASS")

mod.null <- glm.nb(formula=pfu~1, data=phage, link="log")
1 - pchisq(mod.null$deviance, mod.null$df.residual)

mod.1 <- glm.nb(formula=pfu~bottleneck, data=phage, link="log")
1 - pchisq(mod.1$deviance, mod.1$df.residual)

mod.2 <- glm.nb(pfu~timepoint, data=phage, link = "log")
1- pchisq(mod.2$deviance, mod.2$df.residual)

mod.g <- glm.nb(pfu~bottleneck*timepoint, data=phage, link = "log")
1 - pchisq(mod.g$deviance, mod.g$df.residual)

AIC(mod.null, mod.1, mod.2, mod.g) %>% compare_AICs()
anova(mod.null, mod.1, mod.2, mod.g, test="Chisq")

phage$timepoint %<>% relevel(., ref="t5")
phage$timepoint %<>% relevel(., ref="t4")
phage$timepoint %<>% relevel(., ref="t3")
phage$timepoint %<>% relevel(., ref='t2')
phage$timepoint %<>% relevel(., ref='t1')
phage$timepoint %<>% relevel(., ref='t0')

mod.g <- glm.nb(pfu~bottleneck*timepoint, data=phage, link = "log")

summary(mod.g)
model.tables(aov(mod.g), "mean", se=T)


m1 <- glm(pfu~bottleneck*timepoint, data=phage, family=poisson())
summary(m1)

for (i in seq(41,49,1)){
  print(mod.coef[i,3])
}

p <- predict(mod.g, newdata=phage, type = "response", se.fit = T)
upr <- with(p, fit + (1.96*se.fit))
lwr <- with(p, fit - (1.96*se.fit))
CIs <- data.frame(lwr, upr) %>% distinct()

for (i in seq(23,40,1)){
  print(CIs[i,])
}

CIs
summary(mod.g)

library(scales)

sum <- read.csv("~/OneDrive/Data/Bottlenecks/Phage_bn_exp/data/counts/counts_summary_nbinom.csv")
str(sum)
sum %<>% na.exclude
sum$bottleneck %<>% as.factor
sum <- filter(sum, CI!=0)
sum$mean %<>% +1
sum$lower %<>% +1
sum$upper %<>% +1

qplot(y=log10.mean, x=bottleneck, data=sum )+
  geom_point()+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=log10.lower, ymax=log10.upper), size=0.7, width=0.2)+

  #scale_y_continuous(trans = 'log10',
   #                  breaks = trans_breaks("log10", function(x) 10^x),
    #                 labels = trans_format("log10", math_format(10^.x)),
     #                limits = c(1,5e+05))+
  theme_bw()+
  facet_wrap(~timepoint)


