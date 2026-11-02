########## Do years with extreme precipitation anomalies magnify legacy effects in semi-arid grassland ecosystems? #############
# Journal of Ecology
# corresponding author: Kathleen Condon (kathleen.condon@colostate.edu, kcondon@illinois.edu)

# Analysis of discrete variables (ANPP and PRS soil N)

##### Setup =================================================================
library(tidyverse)
library(car)
library(emmeans)
library(lmerTest)
library(rstatix)
library(performance)

plot_theme<-theme_classic()+
  theme(
    strip.background=element_blank(),
    strip.text.x=element_text(size=12,
                              hjust=0,
                              face="bold",colour="black"),
    strip.text.y=element_text(size=12,
                              hjust=0,
                              face="bold",colour="black"),
    panel.background=element_rect(fill="transparent"),
    panel.border=element_rect(fill="transparent",color="black",
                              linewidth=0.5),
    plot.background=element_rect(fill="transparent",color=NA),
    legend.background=element_rect(fill="transparent"),
    legend.position="none",
    # legend.title=element_text(size=16,face="bold"),
    legend.text=element_text(size=12),
    # axis.line.x=element_line(linewidth=1,color="black"),
    axis.line.x=element_blank(),
    axis.line.y=element_blank(),
    axis.ticks=element_line(linewidth=0.5,color="black"),
    axis.ticks.length=unit(-0.05,"in"),
    axis.title = element_text(size=12,
                              colour="black",face="bold"),
    axis.text=element_text(size=10,colour="black"),
    axis.text.x=element_text(angle=0,hjust=0.5,colour="black"),
    axis.text.y=element_text(colour="black"),
    plot.title=element_text(size=12),
    strip.placement="outside")

##### ANPP ====================================================================
anpp<-read.csv("final_data/SGS_EW_ANPP_2022-2024.csv")
str(anpp)

anpp.no.na<-anpp%>%
  # rename variables from original file
  rename(TRT.2022=TRT_22,TRT.2023=TRT_23,C3.GRASSES=C3)%>%
  mutate(TRT.2022=factor(TRT.2022,levels=c("amb","wet")))%>%
  mutate(TRT.2023=factor(TRT.2023,levels=c("amb","wet")))%>%
  # select columns to keep
  select(c(Plot,Year,Clip,TRT.2022,TRT.2023,C3.GRASSES,C4,FORB,ANNUALS,TOTAL.GRASSES,ANPP))%>%
  # drop rows with ANPP NA's
  drop_na(ANPP)%>%
  # create Unique ID for dropping specific observations (outliers)
  mutate(Unique_ID=paste("Plot_",Plot,Clip,sep="-"))

str(anpp.no.na)

anpp.22<-subset(anpp.no.na,Year=="2022")
anpp.23<-subset(anpp.no.na,Year=="2023")
anpp.24<-subset(anpp.no.na,Year=="2024")

#### ANPP outlier checks
#### 2022 ####
dat<-anpp.22
dat$resp<-dat$ANPP

out.plot<-ggplot(dat,aes(x=as.factor(Clip),y=resp,color=TRT.2022))+
  geom_boxplot()+
  scale_color_manual(values=c("darkorange4","steelblue4"))+
  theme_classic()+
  plot_theme+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="right")

out.plot

ext.22<-dat%>%
  group_by(Clip,TRT.2022)%>%
  identify_outliers(resp)%>%
  filter(is.extreme==TRUE)
ext.22

mod.check<-lmer(resp~as.factor(Clip)*TRT.2022+(1|Plot),data=dat)
# performance::check_model(mod.check)

outlierTest(mod.check) # none significant
# view(dat)
dat.22<-dat

#### 2023 ####
dat<-anpp.23%>%
  filter(!Unique_ID %in% c("Plot_-12-Sept_23_clip","Plot_-2-Sept_23_clip","Plot_-31-July_23"))
dat$resp<-dat$ANPP

anpp.23.avs<-dat%>%
  group_by(Year,TRT.2022,TRT.2023)%>%
  dplyr::summarise(
    n.plots=n(),
    av = mean(resp),
    sd = sd(resp),
    SE = sd/sqrt(n.plots))

out.plot<-ggplot(dat,aes(x=as.factor(Clip),y=resp,color=interaction(TRT.2022,TRT.2023)))+
  scale_color_manual(values=c("darkorange4","sandybrown","lightblue","steelblue4"))+
  geom_boxplot()+
  theme_classic()+
  plot_theme+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="right")

out.plot

ext.23<-dat%>%
  group_by(Clip,TRT.2022,TRT.2023)%>%
  identify_outliers(resp)%>%
  filter(is.extreme==TRUE)
ext.23

mod.check<-lmer(resp~as.factor(Clip)*TRT.2022*TRT.2023+(1|Plot),data=dat)
# performance::check_model(mod.check)

outlierTest(mod.check) # outlier @ Plot_-12-Sept_23_clip -- huge (500g/m2) forb skewing ANPP & Plot_-2-Sept_23_clip 230g/m2 forb & Plot_-31-July_23 175g/m2 forb
# view(dat) 
dat.23<-dat

#### 2024 ####
dat<-anpp.24
dat$resp<-dat$ANPP

out.plot<-ggplot(dat,aes(x=as.factor(Clip),y=resp,color=interaction(TRT.2022,TRT.2023)))+
  scale_color_manual(values=c("darkorange4","sandybrown","lightblue","steelblue4"))+
  geom_boxplot()+
  theme_classic()+
  plot_theme+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="right")

out.plot

ext.24<-dat%>%
  group_by(Clip,TRT.2022,TRT.2023)%>%
  identify_outliers(resp)%>%
  filter(is.extreme==TRUE)
ext.24

mod.check<-lm(resp~TRT.2022*TRT.2023,data=dat)
# performance::check_model(mod.check)

outlierTest(mod.check) # no significant outliers
# view(dat) 
dat.24<-dat

#### ANPP end of season
#### 2022 ####
end.22<-subset(dat.22,Clip=="Sept_22_clip")
str(end.22) # 40 observations

## T-test
# check normality
ggplot(end.22, aes(x = resp, fill = TRT.2022)) +
  geom_density(alpha = 0.5)

shapiro.test(end.22$resp[end.22$TRT.2022 == "amb"]) # fine
shapiro.test(end.22$resp[end.22$TRT.2022 == "wet"]) # fine

# if normal, use t.test, if not use transformations or wilcox.test
t.test(resp~TRT.2022,data=end.22) # p < 0.0001

#### 2023 ####
end.23<-subset(dat.23,Clip=="Sept_23_clip")
str(end.23) # 38 observations

#### combined complete ANOVA ####
end.23.comb.mod<-lm(resp ~ TRT.2022*TRT.2023, data=end.23)
plot(end.23.comb.mod,1)
qqnorm(resid(end.23.comb.mod));qqline(resid(end.23.comb.mod))

Anova(end.23.comb.mod,type=3) # interaction NOT significant
pairs(emmeans(end.23.comb.mod,pairwise~TRT.2022|TRT.2023))
p.23.amb<-"Emmeans pairwise p = 0.7828 (t ratio = -0.278, df = 34), interaction main effect n.s."
p.23.wet<-"Emmeans pairwise p = 0.6420 (t ratio = -0.469, df = 34), interaction main effect n.s."

end.23.comb.mod.red<-lm(resp~TRT.2023,data=end.23)
Anova(end.23.comb.mod.red,type=3)
step(end.23.comb.mod)

#### 2024 ####
end.24<-subset(dat.24,Clip=="Sept_24")
str(end.24) # 40 observations

#### combined complete ANOVA ####
end.24.comb.mod<-lm(resp ~ TRT.2022*TRT.2023, data=end.24)
plot(end.24.comb.mod,1)
qqnorm(resid(end.24.comb.mod));qqline(resid(end.24.comb.mod))

Anova(end.24.comb.mod,type=3) # interaction IS significant
pairs(emmeans(end.24.comb.mod,pairwise~TRT.2022|TRT.2023))
p.24.amb<-"Emmeans pairwise p = 0.0262 (t ratio = -2.319, df = 36)"
p.24.wet<-"Emmeans pairwise p = 0.2725 (t ratio = 1.114, df = 36)"
step(end.24.comb.mod) # no reduction

# pairs keeping 2022 treatment the same and looking for differences from 2023
pairs(emmeans(end.24.comb.mod,pairwise~TRT.2023|TRT.2022))
p.24.amb22<-"Emmeans pairwise p = 0.0064 (t ratio = -2.897, df = 36)"
p.24.wet22<-"Emmeans pairwise p = 0.5950 (t ratio = 0.536, df = 36)"


#### ANPP mid season
#### 2022 (watered only - plot only, no stats to do) ####
mid.22<-subset(dat.22,Clip=="July_22")
str(mid.22) # 20 observations

#### 2023 ####
mid.23<-subset(dat.23,Clip=="July_23")
str(mid.23) # 39 observations

#### combined complete ANOVA ####
mid.23.comb.mod<-lm(resp ~ TRT.2022*TRT.2023, data=mid.23)
plot(mid.23.comb.mod,1)
qqnorm(resid(mid.23.comb.mod));qqline(resid(mid.23.comb.mod))

Anova(mid.23.comb.mod,type=3) # interaction NOT significant (0.07, close)
pairs(emmeans(mid.23.comb.mod,pairwise~TRT.2022|TRT.2023))
p.23.amb<-"Emmeans pairwise p = 0.0381 (t ratio = 2.156, df = 35), interaction main effect n.s."
p.23.wet<-"Emmeans pairwise p = 0.6357 (t ratio = -0.478, df = 35), interaction main effect n.s."

step(mid.23.comb.mod) # no reduction

#### ANPP regrowth
#### 2022 (watered only - plot only, no stats) ####
re.22<-subset(dat.22,Clip=="Sept_22_reclip")
str(re.22) # 20 observations

#### 2023 ####
re.23<-subset(dat.23,Clip=="Sept_23_reclip")
str(re.23) # 40 observations

#### combined complete ANOVA ####
re.23.comb.mod<-lm(resp ~ TRT.2022*TRT.2023, data=re.23)
plot(re.23.comb.mod,1)
qqnorm(resid(re.23.comb.mod));qqline(resid(re.23.comb.mod))

Anova(re.23.comb.mod,type=3) # interaction NOT significant
pairs(emmeans(re.23.comb.mod,pairwise~TRT.2022|TRT.2023))
p.23.amb<-"Emmeans pairwise p = 0.2402 (t ratio = -1.194, df = 36), interaction main effect n.s."
p.23.wet<-"Emmeans pairwise p = 0.7524 (t ratio = -0.318, df = 36), interaction main effect n.s."

##### PRS ====================================================================
prs.dat<-read.csv("final_data/SGS_EW_PRS.csv")
str(prs.dat) # 600 observations

prs.no.na<-prs.dat%>%
  # rename variables from original file
  rename(TRT.2022=TRT_2022,TRT.2023=TRT_2023)%>%
  mutate(TRT.2022=factor(TRT.2022,levels=c("amb","wet")))%>%
  mutate(TRT.2023=factor(TRT.2023,levels=c("amb","wet")))%>%
  # select columns to keep
  select(c(Plot,Burial,TRT.2022,TRT.2023,Nutrient,rate,below_MDL,Unique_ID,Notes))%>%
  # drop rows with NA's & below MDL
  drop_na(rate)%>%
  filter(below_MDL == "no")

str(prs.no.na) # 487 observations

prs.no3<-subset(prs.no.na,Nutrient=="NO3_N")
str(prs.no3) # 40 observations


##### PRS outlier check 
dat<-prs.no3

dat$resp<-dat$rate
str(dat)

out.plot<-ggplot(dat,aes(x=as.factor(Burial),y=resp,color=interaction(TRT.2022,TRT.2023)))+
  scale_color_manual(values=c("darkorange4","sandybrown","lightblue","steelblue4"))+
  geom_boxplot()+
  theme_classic()+
  plot_theme+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="right")

out.plot

ext.23<-dat%>%
  group_by(Burial,TRT.2022,TRT.2023)%>%
  identify_outliers(resp)%>%
  filter(is.extreme==TRUE)
ext.23

mod.check<-lmer(resp~as.factor(Burial)*TRT.2022*TRT.2023+(1|Plot),data=dat)
# performance::check_model(mod.check)

outlierTest(mod.check) # no significant outliers

dat.23<-dat

##### PRS analysis =======================================================
#### 2023 ####
str(dat.23) # 40 observations

## Analysis
periods.23<-dat.23%>%
  mutate(period=case_when(
    Burial == 1 ~ "early-mid",
    Burial == 2 ~ "mid-late"))%>%
  mutate(period=factor(period,levels=c("early-mid","mid-late")))
str(periods.23)

wet.23<-subset(periods.23,TRT.2023=="wet")
amb.23<-subset(periods.23,TRT.2023=="amb")

#### full combined ANOVA ####
prs.comb.mod<-lmer(log(resp) ~ period*TRT.2022*TRT.2023 + (1|Plot), data=periods.23)

plot(prs.comb.mod)
qqnorm(resid(prs.comb.mod));qqline(resid(prs.comb.mod))

anova(prs.comb.mod,ddf="Kenward-Roger") # three way interaction NOT significant
pairs(emmeans(prs.comb.mod,pairwise~TRT.2022|TRT.2023|period))

p.23.amb_early<-"Emmeans pairwise p = 0.0070 (t ratio = 2.946, df = 24.1), interaction main effect n.s."
p.23.wet_early<-"Emmeans pairwise p = 0.2481 (t ratio = 1.184, df = 24.1), interaction main effect n.s."
p.23.amb_late<-"Emmeans pairwise p = 0.0367 (t ratio = 2.213, df = 24.1), interaction main effect n.s."
p.23.wet_late<-"Emmeans pairwise p = 0.8710 (t ratio = 0.164, df = 24.1), interaction main effect n.s."

step(prs.comb.mod)

prs.comb.mod.red<-lmer(log(resp)~period+TRT.2022+TRT.2023+(1|Plot),data=periods.23)
anova(prs.comb.mod.red,ddf="Kenward-Roger")

#### across growing season ####
prs.all<-lmer(log(resp)~TRT.2022*TRT.2023+(1|Plot),data=periods.23)
plot(prs.all)
qqnorm(resid(prs.all));qqline(resid(prs.all))

anova(prs.all,ddf="Kenward-Roger")

prs.all.reduced<-lmer(log(resp)~TRT.2022+TRT.2023+(1|Plot),data=periods.23)

anova(prs.all.reduced,ddf="Kenward-Roger")

