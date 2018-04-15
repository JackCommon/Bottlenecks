spacer <- read.csv("./Documents/OneDrive - University of Exeter/Data/Bottlenecks/Plate_bn_exp/data/PCR/spacers_poisson.csv")
spacer %<>% select(., -X)
spacer %<>% na.exclude
spacer$Bottleneck %<>% relevel(ref="Monoculture")
str(spacer)

mod.null <- glm(Total.spacers ~ 1, data=spacer, family=poisson(link="identity"))
summary(mod.null)
pseudo_R2(mod.null, mod.null)

mod.1 <- glm(Total.spacers~Bottleneck, data=spacer, family=poisson(link="identity"))
summary(mod.1)
confint(mod.1)

AIC(mod.null, mod.1) %>% compare_AICs()
pseudo_R2(mod.1, mod.null)

model.tables(aov(mod.1), "mean")

spacer_sum <- read.csv("Documents/OneDrive - University of Exeter/Data/Bottlenecks/Plate_bn_exp/data/PCR/PCR_summary_poiss.csv",
                       header=T)
spacer_sum$bottleneck %<>% relevel(ref="monoculture")

qplot(x=bottleneck, y=mean, data=spacer_sum) +
  geom_point()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2, size=0.7)+
  theme_bw()
