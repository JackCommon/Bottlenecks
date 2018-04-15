pheno <- read.csv("./Documents/OneDrive - University of Exeter/Data/Bottlenecks/Plate_bn_exp/data/phenotype/prop_test.csv")
str(pheno)

pheno %<>% na.exclude

library(reshape2)

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

pheno <- melt(pheno, measure.vars = c('CRISPR', 'SM', 'Sensitive'))
pheno$Bottleneck %<>% relevel(ref="Monoculture")

mod.null <- glm(value~1, data=pheno, family=binomial)
mod.1 <- glm(value~Bottleneck, data=pheno, family=binomial)
mod.2 <- glm(value~variable, data=pheno, family=binomial)
mod.g <- glm(value~Bottleneck*variable, data=pheno, family=binomial)

AIC(mod.null, mod.1, mod.2, mod.g) %>% compare_AICs()
pseudo_R2(mod.1, mod.null)
pseudo_R2(mod.2, mod.null)
pseudo_R2(mod.g, mod.null)

anova(mod.g, test="LRT")

x <- summary(mod.g)
confint(mod.g)

model.tables(aov(mod.g), "mean")

pheno$variable %<>% relevel(ref="SM")
pheno$variable %<>% relevel(ref="Sensitive")
pheno$variable %<>% relevel(ref="CRISPR")

pheno$Bottleneck %<>% relevel(ref="Monoculture")
pheno$Bottleneck %<>% relevel(ref="5-clone")
pheno$Bottleneck %<>% relevel(ref="50-clone")

mod.g <- glm(value~Bottleneck*variable, data=pheno, family=binomial)

logit2prob(mod.g$coefficients)
logit2prob(confint(mod.g))


pheno_sum <- read.csv("./Documents/OneDrive - University of Exeter/Data/Bottlenecks/Plate_bn_exp/data/phenotype/phenotype_prob_summary.csv")
pheno_sum$bottleneck %<>% relevel(ref="monoculture")

str(pheno_sum)

qplot(x=bottleneck, y=mean, data=pheno_sum)+
  geom_point(position=pd, aes(colour=variable))+
  geom_errorbar(aes(ymin=lower, ymax=upper, colour=variable), width=0.2, size=0.7, position = pd)+
  theme_bw()
  