library(dplyr)
library(magrittr)

setwd("~/OneDrive/Data/Bottlenecks/")

phage <- read.csv("~/OneDrive/Data/Bottlenecks/Phage_bn_exp/data/counts/counts_master_csv.csv")
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

summary(mod.g)
model.tables(aov(mod.g), "mean", se=T)

phage$timepoint %<>% relevel(ref="t3")

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

summary(mod.g)

library(scales)

sum <- read.csv("~/OneDrive/Data/Bottlenecks/Phage_bn_exp/data/counts/counts_summary_nbinom.csv")
str(sum)
sum %<>% na.exclude
sum$bottleneck %<>% as.factor
sum <- filter(sum, CI!=0)

qplot(y=mean, x=bottleneck, data=sum)+
  geom_point()+
  geom_errorbar(aes(ymin=lower, ymax=upper))+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = c(1,1e+07))+
  theme_bw()+
  facet_wrap(~timepoint)
