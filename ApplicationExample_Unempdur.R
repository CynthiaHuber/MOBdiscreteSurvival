#Application example

rm(list=ls())
library(Ecdat)
library(tidyverse)
library(discSurv)
library(model4you)
data("UnempDur")

#Inverse coding of event )(censor4==jobless)
UnempDur$censor5 <- ifelse(UnempDur$censor4==1,0,1)
UnempDur2 <- UnempDur%>%dplyr::select(-censor1,-censor2,-censor3,-censor4)


#Current Code of MOB-dS needs ID, trt and CV columns
#Columns are renamed
UnempDur2$ID <- 1:length(UnempDur2$spell)
UnempDur2$CV <- NA
UnempDur2$trt <- NA
colnames(UnempDur2)[1] <- "time.discrete"
colnames(UnempDur2)[colnames(UnempDur2)=="censor5"] <- "state"



#truncate time at time=10
UnempDur2 <- UnempDur2%>%mutate(time.discrete=
                                  ifelse(time.discrete>10,10,
                                         time.discrete))
#Generate augmented data
UnempDur2_L <-  dataLong(UnempDur2 ,timeColumn="time.discrete",censColumn="state")



#MOB
mobALL_2 <-glmtree(y ~ timeInt|age+ui+reprate+disrate+logwage+tenure,
                   data = UnempDur2_L ,
                   family=binomial(link=logit),
                   maxdepth=4,vcov="sandwich")


sctest.modelparty(mobALL_2)

mobALL_2%>%plot

source("MOBdS_function.R")


MOBdS2_28_nresamp10000 <- glmtreeMC(y ~ -1+timeInt|age+ui+reprate+disrate+logwage+tenure,
                                    data = UnempDur2_L,family=binomial(link=logit),
                                    datenShort=UnempDur2,maxdepth=4,nresample = 10000,vcov="sandwich")
MOBdS2_28_nresamp10000%>%plot
sctest.modelparty(MOBdS2_28_nresamp10000)



mobALL_2_depth4 <-glmtree(y ~ timeInt|age+ui+reprate+disrate+logwage+tenure,
                          data = UnempDur2_L ,family=binomial(link=logit),maxdepth=4,vcov="sandwich")


sctest.modelparty(mobALL_2_depth4 )

mobALL_2_depth4 %>%plot




save(MOBdS2_28_nresamp10000,mobALL_2,
     file = "~/MOB_discSurv/Anwendung_H0_V2/application_manuscript_lastinterval_Tlarger10_cens4.RData")


load( "~/MOB_discSurv/Anwendung_H0_V2/application_manuscript_lastinterval_Tlarger10_cens4.RData")
