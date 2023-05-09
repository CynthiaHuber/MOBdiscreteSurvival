#Application example

rm(list=ls())
library(Ecdat)
library(tidyverse)
library(discSurv)
library(model4you)
data("UnempDur")


#Augmented data usingInverse coding of event )(censor4==jobless)

#Current Code of MOB-dS needs ID, trt and CV columns
#Columns are renamed
unemp1 <- UnempDur |>
  transform(
    spell = pmin(spell, 10),
    state = 1 - censor4) |>
  discSurv::dataLong(
    timeColumn = "spell",
    eventColumn = "state",
    timeAsFactor = TRUE) |>
  transform(
    y = factor(y))%>%mutate(ID=NA,CV=NA,trt=NA)%>%dplyr::select(-censor1,-censor2,-censor3,-censor4)
#Column is renamed
colnames(unemp1)[4] <- "time.discrete"



#MOB using the augmented data
tr1 <- glmtree(y ~ timeInt | age + ui + reprate + disrate + logwage + tenure,
               data = unemp1, family = binomial, maxdepth = 4, vcov = "sandwich", cluster = obj)

plot(tr1)

sctest.modelparty(tr1)

#Current Code of MOB-dS needs ID, trt and CV columns in both augmented and short version
#Columns are renamed
UnempDur2 <- UnempDur %>%
  transform(
    spell = pmin(spell, 10)%>%as.factor(),
    state = 1 - censor4)%>%mutate(ID=NA,CV=NA,trt=NA)%>%dplyr::select(-censor1,-censor2,-censor3,-censor4)
colnames(UnempDur2)[1] <- "time.discrete"


#Application of MOB-dS 
source("MOBdS_function.R")


MOBdS2_28_nresamp10000 <- glmtreeMC(y ~ -1+timeInt|age+ui+reprate+disrate+logwage+tenure,
                                    data = unemp1,family=binomial(link=logit),
                                    datenShort=UnempDur2,maxdepth=4,
                                    nresample = 10000,vcov="sandwich")
MOBdS2_28_nresamp10000%>%plot
sctest.modelparty(MOBdS2_28_nresamp10000)


#MOB-sum  with mobster for survival data

unemp2 <- UnempDur %>% transform(
  spell = survival::Surv(pmin(spell, 10), 1 - censor4))
head(unemp2$spell, 10)


ds_augment <- function(y, factor = TRUE) {
  y <- as.matrix(y)
  d <- data.frame(
    id = rep.int(1L:NROW(y), y[, "time"]),
    time = unlist(lapply(y[, "time"], seq.int)),
    event = 0)
  d$event[cumsum(y[, "time"])] <- y[, "status"]
  if(factor) {
    d$time <- factor(d$time)
    d$event <- factor(d$event)
  }
  return(d)
  
}

ds_fit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ...,
                   estfun = FALSE, object = FALSE)
{
  ## augment response
  d <- ds_augment(y)
  
  ## add regressor matrix (if any)
  if(is.null(x) || (NCOL(x) == 1L && colnames(x) == "(Intercept)")) {
    f <- event ~ time
  } else {
    d$x <- x[d$id, , drop = FALSE]
    f <- event ~ time + x
  }
  
  ## fit glm for augmented data
  rval <- glm(f, data = d, family = binomial(link = "logit"))
  
  ## return coefficients, negative log-likelihood, and scores (aggregated within id)
  list(
    coefficients = coef(rval),
    objfun = -as.numeric(logLik(rval)),
    estfun = if(estfun) apply(sandwich::estfun(rval), 2L, rowsum, d$id) else NULL,
    object = if(object) rval else NULL
  )
  
}


mobsum <- mob(spell ~ age + ui + reprate + disrate + logwage + tenure, data = unemp2,
              fit = ds_fit, control = mob_control(maxdepth = 4), vcov="sandwich") 

plot(mobsum)
