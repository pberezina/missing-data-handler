library(VIM)
library(ggplot2)
library(mitools) 
library(tidyr)
library(purrr)
library(mice)
library(dplyr)
library(knitr)
library(reshape2)
library(zoo)

# loading subset of the data 
load('data.Rdata')
wldatasub <- data %>% 
  mutate_each_(funs(factor(.)), c("GENDER", "ETH", "EDU", "EMPLOY")) %>% 
  mutate(FEMALE = recode(GENDER, Female = 1, Male = 0))
subbase <- select(wldatasub, "BL_BMI", "BL_TRI", "BL_CHO", "BL_HDL", "BL_LDL")

# exploring missing data mechanisms and patterns
aggr(subbase, cex.lab=0.8, cex.axis=0.7,xaxs="i", yaxs="i", cex.numbers=0.8, numbers=TRUE)
pbox(subbase, cex.lab=0.8, cex.axis=0.7, interactive = FALSE,)

### 1 
# imputation using MICE predictive mean matching
B <- 100 
lmest.boot <- data.frame(coefficient = rep(c("b0", "b1", "b2", "b3", "b4"), times = B), estimate = rep(NA, B), p.value = rep(NA, B), CI2.5 = rep(NA, B),  CI97.5 = rep(NA, B))
for(b in seq(1, B*5, by=5)){
  imp <- mice(subbase, method="pmm", m=30, maxit=10, donors=3, print=FALSE)
  #plot(imp)
  #library(lattice)
  #stripplot(imp)
  lmfit <- with(imp, lm(BL_BMI ~ BL_TRI+BL_CHO+BL_HDL+BL_LDL))
  lmest <- summary(pool(lmfit), conf.int = TRUE,conf.level = 0.95)
  lmest.boot[b:(b+4), 2] <- lmest$estimate
  lmest.boot[b:(b+4), 3] <- lmest$p.value
  lmest.boot[b:(b+4), 4] <- lmest$`2.5 %`
  lmest.boot[b:(b+4), 5] <- lmest$`97.5 %`
}
bl.completed <- complete(imp)

ggplot(lmest.boot, aes(x=coefficient, y = p.value)) +
  geom_boxplot(fill = "grey", colour = "black") + 
  theme(axis.text=element_text(size=12)) + 
  theme_bw() 

# inference to determine whether lipid profile variables are associated with BMI
lmest.boot %>% 
  group_by(coefficient) %>% 
  summarize(Mean = mean(estimate), '2.5% CI' = mean(CI2.5), '97.5%' = mean(CI97.5)) %>% 
  kable(caption = "Means of coefficient estimates and their 95% CI from MICE-PMM bootsrap", digits = 3)

### 2
wk12 <- select(wldatasub, starts_with("WK12"), -WK12_DAYS_COM)
pbox(wk12, cex.lab=0.8, cex.axis=0.8, interactive = FALSE,)
marginplot(data.frame(Weight=wk12$WK12_WGHT, Cholesterol=wk12$WK12_CHO))
aggr(wk12, numbers = TRUE, cex.axis=0.7, cex.lab=0.8, cex.numbers=0.7,xaxs="i", yaxs="i" )

# imputation using 4 methods: BOCF, LOCF, MICE with normal links, and single imputation
#LOCF, BOCF imputation
wght3 <-  select(wldatasub, ends_with("WGHT"))      #3 cols. all baseline weights were obsereved, no need to use imputed     
wght2 <-  data.frame(BL_WGHT = wldatasub$BL_WGHT, WK12_WGHT = wldatasub$WK12_WGHT) #2 cols

fun_obs_imp <- function(dat){
  n <- nrow(dat)
  m <- ncol(dat)
  for (i in 1:n){
    cz <- zoo(as.numeric(dat[i,]))
    for (j in 1:m){
      imp.vals <- na.locf(cz)
      dat[i,j] <- imp.vals[[j]]
    }}
  dat
}

imp.locf <- fun_obs_imp(wght3)
imp.bocf <- fun_obs_imp(wght2)

# add BMI variable
imp.locf <- mutate(imp.locf, WK12_BMI = WK12_WGHT / (wldatasub$BL_HGT)^2 *10000)
imp.bocf <- mutate(imp.bocf, WK12_BMI = WK12_WGHT / (wldatasub$BL_HGT)^2*10000 )

# single imputation 
imp_mean <- mice(wk12, method="norm.nob", m=1, maxit=1, print=FALSE)
imp.single <- complete(imp_mean)
# add BMI
imp.single <- mutate(imp.single, WK12_BMI = WK12_WGHT / (wldatasub$BL_HGT)^2 *10000 )

# naive, or complete cas imputation
imp.naive <- mutate(wldatasub, WK12_BMI_test = WK12_WGHT/(BL_HGT)^2 * 10000)

lm1 <- lm(WK12_WGHT ~ BL_HGT + FEMALE, wldatasub)
beta.hat <- lm1$coeff
names(beta.hat) <- c("beta0.hat", "beta1.hat", "beta2.hat")
sigma.hat <- sigma(lm1)^2
theta.hat <- c(beta.hat, sigma.hat=sigma.hat)
var.beta.hat <- vcov(lm1)
var.sigma.hat <- 2*sigma.hat^2/(lm1$df.residual) 

D <- 1000
n <- nrow(wldatasub)
wk12.improperMI <- matrix(wldatasub$WK12_WGHT, nrow=D, ncol=n, byrow=TRUE) 
for(i in 1:n){ 
  if (is.na(wldatasub$WK12_WGHT[i])){ 
    beta.hat.d <- MASS::mvrnorm(D, beta.hat, var.beta.hat)
    mean.i <- beta.hat.d[,1] + beta.hat.d[,2]*wldatasub$BL_HGT[i] +beta.hat.d[,3]*wldatasub$FEMALE[i]
    sd.i <- sqrt(pmax(rnorm(D, sigma.hat, sd=sqrt(var.sigma.hat)), 0.1^2))  #a minimum of 0.1
    wk12.improperMI[,i] <- rnorm(D, mean=mean.i, sd=sd.i)
  }}

# analysis for MI imputation completed data
means.mi <- rep(NA,D) 
vars.mi <- rep(NA,D)

for(d in 1:D){
  wk12.BMI <-  wk12.improperMI[d,]/ (wldatasub$BL_HGT)^2 *10000 #convert weight imputed to BMI
  t.mi <- t.test(wldatasub$BL_BMI, wk12.BMI, paired=TRUE)
  means.mi[d] <- t.mi$estimate[[1]]
  vars.mi[d] <- t.mi$stderr ^ 2
}

# 4. Estimate the mean of the posterior distribution of mu0-mu1 using Rubin’s rules. 
ML.postest <- MIcombine(results=as.list(means.mi),
                        variances=as.list(vars.mi))
res.MI <- summary(ML.postest, confint=TRUE)


# inference to estimate the change in BMI
t1 <- t.test(wldatasub$BL_BMI, imp.naive$WK12_BMI_test,paired=TRUE)
t2 <- t.test(wldatasub$BL_BMI, imp.locf$WK12_BMI,paired=TRUE)
t3 <- t.test(wldatasub$BL_BMI, imp.bocf$WK12_BMI,paired=TRUE)
t4 <- t.test(wldatasub$BL_BMI, imp.single$WK12_BMI,paired=TRUE)

data.frame('Avg.change' = c(t1$estimate[[1]], t2$estimate[[1]],t3$estimate[[1]], t4$estimate[[1]], res.MI$results),
           'CI2.5' = c(t1$conf.int[1], t2$conf.int[1],t3$conf.int[1],t4$conf.int[1], res.MI$`(lower`),
           'CI97.5' = c(t1$conf.int[2],t2$conf.int[2],t3$conf.int[2],t4$conf.int[2], res.MI$`upper)`), 
           'Std.error'=c(t1$stderr,t2$stderr,t3$stderr,t4$stderr, res.MI$se),
           'pvalue' = c(t1$p.value,t2$p.value,t3$p.value,t4$p.value, 2*pt(q=res.MI$results/res.MI$se, df=49, lower.tail=FALSE)),
           row.names= c('Naive', 'LOCF', 'BOCF', 'Single', 'Improper imp.')) %>% 
  kable(caption = "Means of coefficient estimates and their 95% CI", digits = 3)