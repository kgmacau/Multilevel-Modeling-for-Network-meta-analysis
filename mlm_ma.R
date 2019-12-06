
library(metafor)
library(lme4)
library(nlme)


################################################
# Two-level model using molly data in 2014
################################################
dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat.molloy2014)
dat[,-c(5:10)]
dat$study <- 1:nrow(dat)


# two-level meta-analytic model using lme funtion
lme0 <- lme(yi ~ 1, random = ~ 1 | study, weights = varFixed(~ vi), 
           data=dat)
lme0


# two-level meta-analytic model using rma funtion
dat$vi2=dat$vi*lme0$sigma^2
rma0 <- rma.mv(yi,vi2, random = ~ 1 |study, data=dat)
rma0


# two-level meta-analytic model using lmer funtion
lmer0 <- lmer(yi ~ 1 + (1 | study), weights = 1/vi, data=dat, 
                 control=lmerControl(check.nobs.vs.nlev="ignore", check.nobs.vs.nRE="ignore"))
lmer0



##### Restrict the level 1 variance to 1
lme1 <- lme(yi ~ 1, random = ~ 1 | study, weights = varFixed(~ vi), 
           control=lmeControl(sigma = 1), data=dat)
lme1


# two-level meta-analytic model using rma funtion
rma1 <- rma.mv(yi,vi, random = ~ 1 |study, data=dat)
rma1




#############################################################
### Three-level data in paper of Konstantopoulos in 2011
#############################################################
kons2011<-dat.konstantopoulos2011


# three-level model
#lme
lme0_kons2011 <- lme(yi~1, random = ~ 1|district/school, 
                    weights=varFixed(~vi),
                    data=kons2011)
lme0_kons2011


# rma.mv
kons2011$vi2=kons2011$vi*lme_kons2011$sigma^2
rma0_kons2011 <- rma.mv(yi, vi2, random = ~ 1 | district/school,
                       data = kons2011)
rma0_kons2011


# lmer
lmer0_kons2011<-lmer(yi ~ 1 + (1 | district/school), weights = 1/vi, data=kons2011, 
                    control=lmerControl(check.nobs.vs.nlev="ignore", check.nobs.vs.nRE="ignore"))
lmer0_kons2011



##### Restrict the level 1 variance to 1
#lme
lme1_kons2011 <- lme(yi~1, random = ~ 1|district/school, 
                    weights=varFixed(~vi), control=lmeControl(sigma = 1),
                    data=kons2011)
lme1_kons2011



# rma.mv
rma1_kons2011 <- rma.mv(yi, vi, random = ~ 1 | district/school,
                       data = kons2011)
rma1_kons2011

