
library(metafor)
library(nlme)
library(lme4)
library(mixmeta)

#############################################################
### Hasselblad multiple treatment comparison data
#############################################################

hasselblad1998<-read.csv("C:/Users/212697818/Desktop/hasselblad1998.csv")

### nlme method
lme0_hasselblad1998 <- lme(yi ~ trt-1, random = ~ 1 | design/study, 
               weights=varFixed(~vi), 
               data=hasselblad1998)
lme0_hasselblad1998


### metafor method
hasselblad1998$vi2=hasselblad1998$vi*lme0_hasselblad1998$sigma^2

rma0_hasselblad1998 <- rma.mv(yi, vi2, mods = ~ trt-1, data = hasselblad1998,
                             random = ~ 1 | design/study)

rma0_hasselblad1998


### lme4 method 
library(lme4)
lmer0_hasselblad1998 <- lmer(yi ~ trt-1 + (1 | design/study), weights = 1/vi, 
                 data=hasselblad1998)
lmer0_hasselblad1998



############################################
### restrict the level-1 variance to 1
############################################
### nlme method
lme1_hasselblad1998 <- lme(yi ~ trt-1, random = ~ 1 | design/study, 
                           weights=varFixed(~vi), control=lmeControl(sigma = 1),
                           data=hasselblad1998)
lme1_hasselblad1998


### metafor method
rma1_hasselblad1998 <- rma.mv(yi, vi, mods = ~ trt-1, data = hasselblad1998,
                             random = ~ 1 | design/study)
rma1_hasselblad1998





# Reference
# http://www.metafor-project.org/doku.php/tips:rma_vs_lm_lme_lmer
# http://www.metafor-project.org/doku.php/analyses:konstantopoulos2011