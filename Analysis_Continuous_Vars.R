########## Do years with extreme precipitation anomalies magnify legacy effects in semi-arid grassland ecosystems? #############
# Journal of Ecology
# corresponding author: Kathleen Condon (kathleen.condon@colostate.edu, kcondon@illinois.edu)

# Analysis of continuous variables (VWC, GCC, Soil Respiration/Efflux)

##### Setup =================================================================
library(tidyverse)
library(lubridate)
library(car)
library(emmeans)
library(lmerTest)
library(rstatix)

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

##### VWC =====================================================================
#### VWC Full data ####
vwc.20<-read.csv("final_data/SGS_EW_20cm_VWC.csv")

vwc.20.no.na<-vwc.20%>%
  # set date formats
  mutate(new.Date=as.Date(Full.Date,format="%m/%d/%y"),yr=year(new.Date),mon=month(new.Date))%>%
  mutate(short.Date=paste("2024",mon,day(new.Date),sep="-"))%>%
  mutate(short.Date=as.Date(short.Date,format="%Y-%m-%d"))%>%
  # rename variables from original file
  rename(TRT.2022=Trt_2022)%>%
  mutate(TRT.2022=factor(TRT.2022,levels=c("amb","wet")))%>%
  rename(TRT.2023=Trt_2023)%>%
  mutate(TRT.2023=factor(TRT.2023,levels=c("amb","wet")))%>%
  # select columns to keep
  select(c(new.Date,short.Date,yr,mon,Plot,TRT.2022,TRT.2023,VWC,Notes))%>%
  # drop rows with VWC NA's
  drop_na(VWC)%>%
  # create Unique_ID for dropping specific observations (outliers)
  mutate(Unique_ID=paste(new.Date,"_Plot",Plot,sep=""))

str(vwc.20) # 1900 observations
str(vwc.20.no.na) # 1900 observations

vwc.20cm.2024<-subset(vwc.20.no.na,yr=="2024")
vwc.20cm.2023<-subset(vwc.20.no.na,yr=="2023")
vwc.20cm.2022<-subset(vwc.20.no.na,yr=="2022")

##### VWC Daily
#### 2022 ####
dat<-vwc.20cm.2022
dat$resp<-dat$VWC

# check outliers
out.plot<-ggplot(dat,aes(x=as.factor(new.Date),y=resp,color=TRT.2022))+
  geom_boxplot()+
  scale_color_manual(values=c("darkorange4","steelblue4"))+
  theme_classic()+
  plot_theme+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="right")
out.plot

ext.22<-dat%>%
  group_by(new.Date,TRT.2022)%>%
  identify_outliers(resp)%>%
  filter(is.extreme==TRUE)
ext.22

mod.check<-lmer(resp~as.factor(new.Date)*TRT.2022+(1|Plot),data=dat)

outlierTest(mod.check) # One outlier 2022-07-05_Plot23 - not totally unreasonable in boxplot so leave in
# view(dat)
dat.22<-dat

#### 2022 mod ####
mod.check<-lmer(log(resp)~as.factor(new.Date)*TRT.2022+(1|Plot),data=dat)
plot(mod.check)
qqnorm(resid(mod.check));qqline(resid(mod.check))

anova(mod.check,ddf="Kenward-Roger") # date-trt interaction p<2.2.e-16

mod.pairs<-pairs(emmeans(mod.check,pairwise~TRT.2022|new.Date))
pairs.df<-as.data.frame(summary(mod.pairs))
vwc.22.sig.pairs<-subset(pairs.df,p.value<0.05)

#### 2023 ####
dat<-vwc.20cm.2023%>%
  filter(!Unique_ID %in% c("2023-08-29_Plot3"))
dat$resp<-dat$VWC

# check outliers
out.plot<-ggplot(dat,aes(x=as.factor(new.Date),y=resp,color=interaction(TRT.2022,TRT.2023)))+
  # scale_color_manual(values=c("orange","navy"))+
  scale_color_manual(values=c("darkorange4","sandybrown","lightblue","steelblue4"))+
  geom_boxplot()+
  theme_classic()+
  plot_theme+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="right")
out.plot

ext.23<-dat%>%
  group_by(new.Date,TRT.2022,TRT.2023)%>%
  identify_outliers(resp)%>%
  filter(is.extreme==TRUE )
ext.23

mod.check<-lmer(resp~as.factor(new.Date)*TRT.2022*TRT.2023+(1|Plot),data=dat)

outlierTest(mod.check) # two outliers 2023-08-29_Plot3 & 2023-08-29_Plot8 (both in the same trt at opposite ends) - removing high outlier (Plot 3) also stops second from being outlier
# view(dat)
dat.23<-dat

#### 2023 mod ####
mod.check<-lmer(resp~as.factor(new.Date)*TRT.2022*TRT.2023+(1|Plot),data=dat)
plot(mod.check)
qqnorm(resid(mod.check));qqline(resid(mod.check))

anova(mod.check,ddf="Kenward-Roger")
# step(mod.check) # no reduction

mod.pairs<-pairs(emmeans(mod.check,~TRT.2022|TRT.2023*new.Date))
pairs.df<-as.data.frame(summary(mod.pairs))
vwc.23.sig.pairs<-subset(pairs.df,p.value<0.05)

#### 2024 ####
dat<-vwc.20cm.2024
dat$resp<-dat$VWC

# check outliers
out.plot<-ggplot(dat,aes(x=as.factor(new.Date),y=resp,color=interaction(TRT.2022,TRT.2023)))+
  # scale_color_manual(values=c("orange","navy"))+
  scale_color_manual(values=c("darkorange4","sandybrown","lightblue","steelblue4"))+
  geom_boxplot()+
  theme_classic()+
  plot_theme+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="right")
out.plot

ext.24<-dat%>%
  group_by(new.Date,TRT.2022,TRT.2023)%>%
  identify_outliers(resp)%>%
  filter(is.extreme==TRUE )
ext.24

mod.check<-lmer(resp~as.factor(new.Date)*TRT.2022*TRT.2023+(1|Plot),data=dat)

outlierTest(mod.check) # no outliers
# view(dat)
dat.24<-dat

#### 2024 mod ####
mod.check<-lmer(log(resp)~as.factor(new.Date)*TRT.2022*TRT.2023+(1|Plot),data=dat)
plot(mod.check)
qqnorm(resid(mod.check));qqline(resid(mod.check))

anova(mod.check,ddf="Kenward-Roger")

step(mod.check)

mod.red<-lmer(log(resp)~as.factor(new.Date)+TRT.2022+TRT.2023+as.factor(new.Date)*TRT.2023+(1|Plot),data=dat)
plot(mod.red)
qqnorm(resid(mod.red));qqline(resid(mod.red))
anova(mod.red,ddf="Kenward-Roger")

mod.pairs<-pairs(emmeans(mod.red,~TRT.2023|new.Date))
pairs.df<-as.data.frame(summary(mod.pairs))
vwc.24.sig.pairs<-subset(pairs.df,p.value<0.05)


#####  VWC comparisons for growing seasons
#### 2022 growing season - Fig 2a ####
str(dat.22) # 800 observations

## Analysis
gs.avs.22<-dat.22%>%
  group_by(Plot,TRT.2022)%>%
  dplyr::summarise(
    n.samp.dates = n(),
    season.av = mean(resp),
    season.sd = sd(resp),
    season.SE = season.sd/sqrt(n.samp.dates))

## T-test
ggplot(gs.avs.22, aes(x = season.av, fill = TRT.2022)) +
  geom_density(alpha = 0.5)

shapiro.test(gs.avs.22$season.av[gs.avs.22$TRT.2022 == "amb"]) # fine
shapiro.test(gs.avs.22$season.av[gs.avs.22$TRT.2022 == "wet"]) # fine

t.pval<-t.test(season.av~TRT.2022,data=gs.avs.22) # p < 0.0001
vwc.t.22.pval<-t.pval$p.value

#### 2023 growing season ####
str(dat.23) # 739 observations

## Analysis
gs.avs.23<-dat.23%>%
  group_by(Plot,TRT.2022,TRT.2023)%>%
  dplyr::summarise(
    n.samp.dates = n(),
    season.av = mean(resp),
    season.sd = sd(resp),
    season.SE = season.sd/sqrt(n.samp.dates))

## ANOVA
gs.23<-lm(season.av~TRT.2022*TRT.2023,data=gs.avs.23)
plot(gs.23,1)
qqnorm(resid(gs.23));qqline(resid(gs.23))

Anova(gs.23,type=3) # interaction not significant
pairs(emmeans(gs.23,pairwise~TRT.2022|TRT.2023)) # amb 2023 p = 0.0580 for legacy, wet 2023 p = 0.0113 for legacy

mod.pairs<-pairs(emmeans(gs.23,~TRT.2022|TRT.2023))
vwc.gs.23.pairs<-as.data.frame(summary(mod.pairs))

gs.23.ni<-lm(season.av~TRT.2022+TRT.2023,data=gs.avs.23)
Anova(gs.23.ni,type=3) # both significant

#### 2024 growing season - Figure 2c ####
str(dat.24) # 360 observations

## Analysis
gs.avs.24<-dat.24%>%
  group_by(Plot,TRT.2022,TRT.2023)%>%
  dplyr::summarise(
    n.samp.dates = n(),
    season.av = mean(resp),
    season.sd = sd(resp),
    season.SE = season.sd/sqrt(n.samp.dates))

## ANOVA
gs.24<-lm(season.av~TRT.2022*TRT.2023,data=gs.avs.24)
qqnorm(resid(gs.24));qqline(resid(gs.24)) # fine
plot(gs.24,1)

Anova(gs.24,type=3) # interaction not significant, main effects not significant
mod.pairs<-pairs(emmeans(gs.24,~TRT.2022|TRT.2023)) # amb 2023 0.1972 leg., wet 2023 0.0742 leg.
vwc.gs.24.og.pairs<-as.data.frame(summary(mod.pairs))

# add rearranged bars pairwise for revisions
mod.pairs<-pairs(emmeans(gs.24,~TRT.2023|TRT.2022))
vwc.gs.24.pairs<-as.data.frame(summary(mod.pairs))

# check reduced model
step(gs.24)
gs.24.ni<-lm(season.av~TRT.2022+TRT.2023,data=gs.avs.24)
Anova(gs.24.ni,type=3) # TRT.2022 significant only

gs.24.red<-lm(season.av~TRT.2022,data=gs.avs.24)
Anova(gs.24.red,type=3)

##### VWC mid and late periods
#### 2022 ####
str(dat.22) # 800 observations

## Analysis
periods.22<-dat.22%>%
  mutate(period=case_when(
    new.Date <= ymd("2022-07-06") ~ "early-mid",
    new.Date > ymd("2022-07-06") ~ "mid-late"))%>%
  mutate(period=factor(period,levels=c("early-mid","mid-late")))
str(periods.22)

## full combined seasonal ANOVAs
log.per.22<-lmer(log(resp)~TRT.2022*period+(1|Plot/period),data=periods.22)
plot(log.per.22)
qqnorm(resid(log.per.22));qqline(resid(log.per.22))

anova(log.per.22,ddf="Kenward-Roger") # drop interaction
pairs(emmeans(log.per.22,pairwise~TRT.2022|period)) # both < 0.0001

step(log.per.22,reduce.random=FALSE) # keep 2022 trt only 

log.per.22.red<-lmer(log(resp)~TRT.2022+(1|Plot/period),data=periods.22)
anova(log.per.22.red,ddf="Kenward-Roger")

anova(log.per.22,log.per.22.red) # prefer reduced model

# check pairwise comps of best model
pairs(emmeans(log.per.22.red,pairwise~TRT.2022)) # sig. different (<0.0001) averaged across periods
# check pairwise comps of full
pairs(emmeans(log.per.22,pairwise~TRT.2022|period)) # both are < 0.0001 alone 

#### 2023 - Figure 2b ####
str(dat.23) # 739 observations

## Analysis
periods.23<-dat.23%>%
  mutate(period=case_when(
    new.Date <= ymd("2023-07-06") ~ "early-mid",
    new.Date > ymd("2023-07-06") ~ "mid-late"))%>%
  mutate(period=factor(period,levels=c("early-mid","mid-late")))
str(periods.23)

period.avs.23<-periods.23%>%
  group_by(period,Plot,TRT.2022,TRT.2023)%>%
  dplyr::summarise(
    n.samp.dates= n(),
    season.av = mean(resp),
    season.sd = sd(resp),
    season.SE = season.sd/sqrt(n.samp.dates))

## full combined seasonal ANOVAs
season.23.avs.mod<-lmer(season.av ~ period*TRT.2022*TRT.2023+(1|Plot),data=period.avs.23)
plot(season.23.avs.mod)
qqnorm(resid(season.23.avs.mod));qqline(resid(season.23.avs.mod))
anova(season.23.avs.mod,ddf="Kenward-Roger") # three way interaction NOT significant
mod.pairs<-pairs(emmeans(season.23.avs.mod,~TRT.2022|TRT.2023*period))
vwc.period.23.pairs<-as.data.frame(summary(mod.pairs))

step(season.23.avs.mod)
season.23.avs.mod.red<-lmer(season.av ~ period+TRT.2022+TRT.2023+(1|Plot),data=period.avs.23)
plot(season.23.avs.mod.red)
qqnorm(resid(season.23.avs.mod.red));qqline(resid(season.23.avs.mod.red))
anova(season.23.avs.mod.red,ddf="Kenward-Roger")

#### 2024 ####
str(dat.24) # 360 observations

## Analysis
periods.24<-dat.24%>%
  mutate(period=case_when(
    new.Date <= ymd("2024-07-06") ~ "early-mid",
    new.Date > ymd("2024-07-06") ~ "mid-late"))%>%
  mutate(period=factor(period,levels=c("early-mid","mid-late")))
str(periods.24)

period.avs.24<-periods.24%>%
  group_by(period,Plot,TRT.2022,TRT.2023)%>%
  dplyr::summarise(
    n.samp.dates = n(),
    season.av = mean(resp),
    season.sd = sd(resp),
    season.SE = season.sd/sqrt(n.samp.dates))

## full combined season ANOVAs
season.24.avs.mod<-lmer(season.av ~ period*TRT.2022*TRT.2023+(1|Plot),data=period.avs.24)
plot(season.24.avs.mod)
qqnorm(resid(season.24.avs.mod));qqline(resid(season.24.avs.mod))
anova(season.24.avs.mod,ddf="Kenward-Roger") # only trt 2022 main effect is significant
pairs(emmeans(season.24.avs.mod,~TRT.2022|TRT.2023*period))

step(season.24.avs.mod)
season.24.avs.mod.red<-lmer(season.av ~ TRT.2022 + (1|Plot), data = period.avs.24)
plot(season.24.avs.mod.red)
qqnorm(resid(season.24.avs.mod.red));qqline(resid(season.24.avs.mod.red))
anova(season.24.avs.mod.red,ddf="Kenward-Roger")

##### GCC =====================================================================
#### GCC Full data ####
gcc<-read.csv("final_data/SGS_EW_GCC_2022-2024.csv")
str(gcc)

gcc.no.na<-gcc%>%
  # set date formats
  mutate(new.Date=as.Date(Date,format="%m/%d/%y"),yr=year(new.Date),mon=month(new.Date))%>%
  mutate(short.Date=paste("2024",mon,day(new.Date),sep="-"))%>%
  mutate(short.Date=as.Date(short.Date,format="%Y-%m-%d"))%>%
  # rename variables from original file
  rename(TRT.2022=TRT_2022)%>%
  mutate(TRT.2022=factor(TRT.2022,levels=c("amb","wet")))%>%
  rename(TRT.2023=TRT_2023)%>%
  mutate(TRT.2023=factor(TRT.2023,levels=c("amb","wet")))%>%
  # select columns to keep
  rename(Notes=Notes.)%>%
  select(c(new.Date,short.Date,yr,mon,Plot,TRT.2022,TRT.2023,GCC,Notes))%>%
  # drop rows with VWC NA's
  drop_na(GCC)%>%
  # create Unique_ID for dropping specific observations (outliers)
  mutate(Unique_ID=paste(new.Date,"_Plot",Plot,sep=""))

str(gcc) # 1880 observations
str(gcc.no.na) # 19870 observations

gcc.2024<-subset(gcc.no.na,yr=="2024")
gcc.2023<-subset(gcc.no.na,yr=="2023")
gcc.2022<-subset(gcc.no.na,yr=="2022")

##### GCC Daily
#### 2022 ####
dat<-gcc.2022%>%
  filter(!Unique_ID %in% c("2022-08-11_Plot9"))
dat$resp<-dat$GCC

out.plot<-ggplot(dat,aes(x=as.factor(new.Date),y=resp,color=TRT.2022))+
  geom_boxplot()+
  scale_color_manual(values=c("darkorange4","steelblue4"))+
  theme_classic()+
  plot_theme+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="right")

out.plot

ext.22<-dat%>%
  group_by(new.Date,TRT.2022)%>%
  identify_outliers(resp)%>%
  filter(is.extreme==TRUE)
ext.22

mod.check<-lmer(resp~as.factor(new.Date)*TRT.2022+(1|Plot),data=dat)
plot(mod.check)
qqnorm(resid(mod.check));qqline(resid(mod.check))

outlierTest(mod.check) # One outlier 2022-08-11_Plot9 (dropped). Second check outlier 2022-08-11_Plot2 (not dropped) and 2022-08-11_Plot36 (not dropped)
# view(dat)
dat.22<-dat

#### 2022 mod ####
mod.check<-lmer(resp~as.factor(new.Date)*TRT.2022+(1|Plot),data=dat)
plot(mod.check)
qqnorm(resid(mod.check));qqline(resid(mod.check))

anova(mod.check,ddf="Kenward-Roger") # all sig.

mod.pairs<-pairs(emmeans(mod.check,pairwise~TRT.2022|new.Date))
pairs.df<-as.data.frame(summary(mod.pairs))
gcc.22.sig.pairs<-subset(pairs.df,p.value<0.05)

#### 2023 ####
dat<-gcc.2023%>%
  filter(!Unique_ID %in% c(
    "2023-06-28_Plot2",
    "2023-07-11_Plot2",
    "2023-07-11_Plot12",
    "2023-06-28_Plot12",
    "2023-07-04_Plot2",
    "2023-07-04_Plot12",
    "2023-07-18_Plot2"
  ))
str(dat)

dat$resp<-dat$GCC

out.plot<-ggplot(dat,aes(x=as.factor(new.Date),y=resp,color=interaction(TRT.2022,TRT.2023)))+
  # scale_color_manual(values=c("orange","navy"))+
  scale_color_manual(values=c("darkorange4","sandybrown","lightblue","steelblue4"))+
  geom_boxplot()+
  theme_classic()+
  plot_theme+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="right")
out.plot

ext.23<-dat%>%
  group_by(new.Date,TRT.2022,TRT.2023)%>%
  identify_outliers(resp)%>%
  filter(is.extreme==TRUE )
ext.23

mod.check<-lmer(resp~as.factor(new.Date)*TRT.2022*TRT.2023+(1|Plot),data=dat)
plot(mod.check)
qqnorm(resid(mod.check));qqline(resid(mod.check))

outlierTest(mod.check) # first outlier 2023-06-28_Plot2, second 2023-07-11_Plot2,  third 2023-07-11_Plot12, 2023-06-28_Plot12, 2023-07-04_Plot2 --- all removed (see above)... 2023-06-19_Plot2 also outlier but kept b/c out.plot shows not as unreasonable 

# view(dat)
dat.23<-dat

#### 2023 mod ####
mod.check<-lmer(resp~as.factor(new.Date)*TRT.2022*TRT.2023+(1|Plot),data=dat)
plot(mod.check)
qqnorm(resid(mod.check));qqline(resid(mod.check))

anova(mod.check,ddf="Kenward-Roger") # three way is sig.
# step(mod.check) # no reduction

mod.pairs<-pairs(emmeans(mod.check,~TRT.2022|TRT.2023*new.Date))
pairs.df<-as.data.frame(summary(mod.pairs))
gcc.23.sig.pairs<-subset(pairs.df,p.value<0.05)

#### 2024 ####
dat<-gcc.2024
dat$resp<-dat$GCC

out.plot<-ggplot(dat,aes(x=as.factor(new.Date),y=resp,color=interaction(TRT.2022,TRT.2023)))+
  # scale_color_manual(values=c("orange","navy"))+
  scale_color_manual(values=c("darkorange4","sandybrown","lightblue","steelblue4"))+
  geom_boxplot()+
  theme_classic()+
  plot_theme+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="right")
out.plot

ext.24<-dat%>%
  group_by(new.Date,TRT.2022,TRT.2023)%>%
  identify_outliers(resp)%>%
  filter(is.extreme==TRUE )
ext.24

mod.check<-lmer(resp~as.factor(new.Date)*TRT.2022*TRT.2023+(1|Plot),data=dat)
plot(mod.check)
qqnorm(resid(mod.check));qqline(resid(mod.check))

outlierTest(mod.check) # no outliers
# view(dat)
dat.24<-dat

#### 2024 mod ####
mod.check<-lmer(resp~as.factor(new.Date)*TRT.2022*TRT.2023+(1|Plot),data=dat)
plot(mod.check)
qqnorm(resid(mod.check));qqline(resid(mod.check))

anova(mod.check,ddf="Kenward-Roger") # only date significant

##### GCC comparisons for growing seasons
#### 2022 growing season - Figure 3d ####
str(dat.22) # 793 observations

## Analysis
gs.avs.22<-dat.22%>%
  group_by(Plot,TRT.2022)%>%
  dplyr::summarise(
    n.samp.dates = n(),
    season.av = mean(resp),
    season.sd = sd(resp),
    season.SE = season.sd/sqrt(n.samp.dates))

## T-test
ggplot(gs.avs.22,aes(x = season.av, fill = TRT.2022))+
  geom_density(alpha = 0.5)

shapiro.test(gs.avs.22$season.av[gs.avs.22$TRT.2022 == "amb"]) # fine
shapiro.test(gs.avs.22$season.av[gs.avs.22$TRT.2022 == "wet"]) # fine

t.pval<-t.test(season.av ~ TRT.2022, data = gs.avs.22) # p < 0.0001
gcc.t.22.pval<-t.pval$p.value

#### 2023 growing season ####
str(dat.23) # 749 observations

## Analysis
gs.avs.23<-dat.23%>%
  group_by(Plot,TRT.2022,TRT.2023)%>%
  dplyr::summarise(
    n.samp.dates = n(),
    season.av = mean(resp),
    season.sd = sd(resp),
    season.SE = season.sd/sqrt(n.samp.dates))

## ANOVA
gs.23<-lm(season.av~TRT.2022*TRT.2023,data=gs.avs.23)
plot(gs.23,1)
qqnorm(resid(gs.23));qqline(resid(gs.23))

Anova(gs.23,type=3) # interaction not significant
mod.pairs<-pairs(emmeans(gs.23,pairwise~TRT.2022|TRT.2023)) # amb legacy = 0.2285 n.s., wet = 0.1121 n.s. 
gcc.gs.23.pairs<-as.data.frame(summary(mod.pairs))

gs.23.ni<-lm(season.av~TRT.2022+TRT.2023,data=gs.avs.23)
Anova(gs.23.ni,type=3) # both significant

#### 2024 growing season - Figure 3f ####
str(dat.24) # 360 observations

## Analysis
gs.avs.24<-dat.24%>%
  group_by(Plot,TRT.2022,TRT.2023)%>%
  dplyr::summarise(
    n.samp.dates = n(),
    season.av = mean(resp),
    season.sd = sd(resp),
    season.SE = season.sd/sqrt(n.samp.dates))

## ANOVA
gs.24<-lm(season.av~TRT.2022*TRT.2023,data=gs.avs.24)
qqnorm(resid(gs.24));qqline(resid(gs.24))
plot(gs.24,1)

Anova(gs.24,type=3) # interaction not significant, main effects not significant
mod.pairs<-pairs(emmeans(gs.24,pairwise~TRT.2022|TRT.2023)) # p = 0.5884 amb, 0.9227 wet legacies
gcc.gs.24.og.pairs<-as.data.frame(summary(mod.pairs))

# add rearranged bars pairwise
mod.pairs<-pairs(emmeans(gs.24,~TRT.2023|TRT.2022))
gcc.gs.24.pairs<-as.data.frame(summary(mod.pairs))

step(gs.24) # nothing
gs.24.ni<-lm(season.av~TRT.2022+TRT.2023,data=gs.avs.24)
Anova(gs.24.ni,type=3) # none significant

##### GCC mid and late periods
#### 2022 ####
str(dat.22) # 793 observations

## Analysis
periods.22<-dat.22%>%
  mutate(period=case_when(
    new.Date <= ymd("2022-07-06") ~ "early-mid",
    new.Date > ymd("2022-07-06") ~ "mid-late"))%>%
  mutate(period=factor(period,levels=c("early-mid","mid-late")))
str(periods.22)

## full combined seasonal ANOVA ##
log.per.22<-lmer(log(resp)~TRT.2022*period+(1|Plot/period),data=periods.22)
plot(log.per.22)
qqnorm(resid(log.per.22));qqline(resid(log.per.22))

anova(log.per.22,ddf="Kenward-Roger") # all significant
step(log.per.22,reduce.random=FALSE) # keep all

# check pairwise comps of best model
pairs(emmeans(log.per.22,pairwise~TRT.2022|period)) # sig. different for both periods (< 0.0001)

#### 2023 - Figure 3e ####
str(dat.23) # 749 observations

## Analysis
periods.23<-dat.23%>%
  mutate(period=case_when(
    new.Date <= ymd("2023-07-06") ~ "early-mid",
    new.Date > ymd("2023-07-06") ~ "mid-late"))%>%
  mutate(period=factor(period,levels=c("early-mid","mid-late")))
str(periods.23)

period.avs.23<-periods.23%>%
  group_by(period,Plot,TRT.2022,TRT.2023)%>%
  dplyr::summarise(
    n.samp.dates= n(),
    season.av = mean(resp),
    season.sd = sd(resp),
    season.SE = season.sd/sqrt(n.samp.dates))

## full combined seasonal ANOVAs
season.23.avs.mod<-lmer(season.av ~ period*TRT.2022*TRT.2023+(1|Plot),data=period.avs.23)
plot(season.23.avs.mod)
qqnorm(resid(season.23.avs.mod));qqline(resid(season.23.avs.mod))
anova(season.23.avs.mod,ddf="Kenward-Roger") # three way interaction IS significant
mod.pairs<-pairs(emmeans(season.23.avs.mod,~TRT.2022|TRT.2023*period))
gcc.period.23.pairs<-as.data.frame(mod.pairs)

step(season.23.avs.mod) # no reduction

#### 2024 ####
str(dat.24) # 320 observations

## Analysis
periods.24<-dat.24%>%
  mutate(period=case_when(
    new.Date <= ymd("2024-07-06") ~ "early-mid",
    new.Date > ymd("2024-07-06") ~ "mid-late"))%>%
  mutate(period=factor(period,levels=c("early-mid","mid-late")))
str(periods.24)

period.avs.24<-periods.24%>%
  group_by(period,Plot,TRT.2022,TRT.2023)%>%
  dplyr::summarise(
    n.samp.dates = n(),
    season.av = mean(resp),
    season.sd = sd(resp),
    season.SE = season.sd/sqrt(n.samp.dates))

## full combined season ANOVAs
season.24.avs.mod<-lmer(season.av ~ period*TRT.2022*TRT.2023+(1|Plot),data=period.avs.24)
plot(season.24.avs.mod)
qqnorm(resid(season.24.avs.mod));qqline(resid(season.24.avs.mod))
anova(season.24.avs.mod,ddf="Kenward-Roger") # three way interaction NOT sig.
pairs(emmeans(season.24.avs.mod,~TRT.2022|TRT.2023*period))

step(season.24.avs.mod) # keep only period
season.24.avs.mod.red<-lmer(season.av ~ period+(1|Plot),data=period.avs.24)
plot(season.24.avs.mod.red)
qqnorm(resid(season.24.avs.mod.red));qqline(resid(season.24.avs.mod.red))
anova(season.24.avs.mod.red,ddf="Kenward-Roger")

##### Respiration/Efflux ======================================================
#### Efflux Full data ####
resp<-read.csv("final_data/SGS_EW_Efflux_2022-2023.csv")
str(resp)

resp[resp=="ambient"]<-"amb"

resp.no.na<-resp%>%
  # set date formats
  mutate(new.Date=as.Date(Full.Date,format="%m/%d/%y"),yr=year(new.Date),mon=month(new.Date))%>%
  mutate(short.Date=paste("2024",mon,day(new.Date),sep="-"))%>%
  mutate(short.Date=as.Date(short.Date,format="%Y-%m-%d"))%>%
  # rename variables from original file
  rename(TRT.2022=TRT_2022)%>%
  mutate(TRT.2022=factor(TRT.2022,levels=c("amb","wet")))%>%
  rename(TRT.2023=TRT_2023)%>%
  mutate(TRT.2023=factor(TRT.2023,levels=c("amb","wet")))%>%
  # select columns to keep
  select(c(new.Date,short.Date,yr,mon,Plot,TRT.2022,TRT.2023,Soil.temp,Efflux,Notes))%>%
  # drop rows with VWC NA's
  drop_na(Efflux)%>%
  # create Unique_ID for dropping specific observations (outliers)
  mutate(Unique_ID=paste(new.Date,"_Plot",Plot,sep=""))

str(resp) # 721 observations
str(resp.no.na) # 712 observations

resp.2023<-subset(resp.no.na,yr=="2023")
resp.2022<-subset(resp.no.na,yr=="2022")

##### Resp Daily
#### 2022 ####
dat<-resp.2022
dat$resp<-dat$Efflux

# check outliers
out.plot<-ggplot(dat,aes(x=as.factor(new.Date),y=resp,color=TRT.2022))+
  geom_boxplot()+
  scale_color_manual(values=c("darkorange4","steelblue4"))+
  theme_classic()+
  plot_theme+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="right")

out.plot

ext.22<-dat%>%
  group_by(new.Date,TRT.2022)%>%
  identify_outliers(resp)%>%
  filter(is.extreme==TRUE)
ext.22

mod.check<-lmer(resp~as.factor(new.Date)*TRT.2022+(1|Plot),data=dat)
plot(mod.check)
qqnorm(resid(mod.check));qqline(resid(mod.check))

outlierTest(mod.check) # two outliers but not on plot, so leaving both in (one has note that ran again in field to check high and was high again so OK)

dat.22<-dat

#### 2022 mod ####
mod.check<-lmer(sqrt(resp)~as.factor(new.Date)*TRT.2022+(1|Plot),data=dat)
plot(mod.check)
qqnorm(resid(mod.check));qqline(resid(mod.check))

anova(mod.check,ddf="Kenward-Roger") # interaction sig

mod.pairs<-pairs(emmeans(mod.check,pairwise~TRT.2022|new.Date))
pairs.df<-as.data.frame(summary(mod.pairs))
resp.22.sig.pairs<-subset(pairs.df,p.value<0.05)

#### 2023 ####
dat<-resp.2023%>%
  filter(!Unique_ID %in% c("2023-06-28_Plot35"))
str(dat)

dat$resp<-dat$Efflux

# check outliers
out.plot<-ggplot(dat,aes(x=as.factor(new.Date),y=resp,color=interaction(TRT.2022,TRT.2023)))+
  scale_color_manual(values=c("darkorange4","sandybrown","lightblue","steelblue4"))+
  geom_boxplot()+
  theme_classic()+
  plot_theme+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="right")
out.plot

ext.23<-dat%>%
  group_by(new.Date,TRT.2022,TRT.2023)%>%
  identify_outliers(resp)%>%
  filter(is.extreme==TRUE )
ext.23

mod.check<-lmer(resp~as.factor(new.Date)*TRT.2022*TRT.2023+(1|Plot),data=dat)
plot(mod.check)
qqnorm(resid(mod.check));qqline(resid(mod.check))

outlierTest(mod.check) # outlier on 6/28 had note about coller looking extremely wet from watering day before still, second reading was also high BUT dropped b/c seems like it wasn't drying same as other plots for some reason that day                                                                        
# view(dat)
dat.23<-dat

#### 2023 mod ####
mod.check<-lmer(sqrt(resp)~as.factor(new.Date)*TRT.2022*TRT.2023+(1|Plot),data=dat)
plot(mod.check)
qqnorm(resid(mod.check));qqline(resid(mod.check)) 

anova(mod.check,ddf="Kenward-Roger") # both trt's with date interaction sig, but not three way interaction

mod.pairs<-pairs(emmeans(mod.check,~TRT.2022|TRT.2023*new.Date))
pairs.df<-as.data.frame(summary(mod.pairs))
resp.23.sig.pairs<-subset(pairs.df,p.value<0.05)

step(mod.check)

mod.red<-lmer(sqrt(resp)~as.factor(new.Date)+TRT.2022+TRT.2023+as.factor(new.Date)*TRT.2022+as.factor(new.Date)*TRT.2023+(1|Plot),data=dat)
plot(mod.red)
qqnorm(resid(mod.red));qqline(resid(mod.red))

anova(mod.red,ddf="Kenward-Roger")
mod.pairs<-pairs(emmeans(mod.red,~TRT.2022|new.Date))
pairs.df<-as.data.frame(summary(mod.pairs))
resp.23.sig.pairs.across23trt<-subset(pairs.df,p.value<0.05)

##### Resp comparisons for growing seasons
#### 2022 growing season ####
str(dat.22) # 373 observations

## Analysis
gs.avs.22<-dat.22%>%
  group_by(Plot,TRT.2022)%>%
  dplyr::summarise(
    n.samp.dates = n(),
    season.av = mean(resp),
    season.sd = sd(resp),
    season.SE = season.sd/sqrt(n.samp.dates))

## T-test
ggplot(gs.avs.22, aes(x = season.av, fill = TRT.2022)) +
  geom_density(alpha = 0.5)

shapiro.test(gs.avs.22$season.av[gs.avs.22$TRT.2022 == "amb"]) # fine
shapiro.test(gs.avs.22$season.av[gs.avs.22$TRT.2022 == "wet"]) # fine

t.pval<-t.test(season.av~TRT.2022,data=gs.avs.22) # p < 0.0001
resp.t.22.pval<-t.pval$p.value

#### 2023 growing season ####
str(dat.23) # 338 observations

## Analysis
gs.avs.23<-dat.23%>%
  group_by(Plot,TRT.2022,TRT.2023)%>%
  dplyr::summarise(
    n.samp.dates = n(),
    season.av = mean(resp),
    season.sd = sd(resp),
    season.SE = season.sd/sqrt(n.samp.dates))

## ANOVA
gs.23<-lm(season.av~TRT.2022*TRT.2023,data=gs.avs.23)
plot(gs.23,1)
qqnorm(resid(gs.23));qqline(resid(gs.23))

Anova(gs.23,type=3) # interaction not significant
mod.pairs<-pairs(emmeans(gs.23,pairwise~TRT.2022|TRT.2023)) # amb 2023 legacy is 0.2135, wet 2023 legacy is 0.0580
resp.gs.23.pairs<-as.data.frame(summary(mod.pairs))

gs.23.ni<-lm(season.av~TRT.2022+TRT.2023,data=gs.avs.23)
Anova(gs.23.ni,type=3) # both significant (2023 0.066)
pairs(emmeans(gs.23.ni,pairwise~TRT.2022)) # p = 0.0275 legacy effect across 2023 trts

step(gs.23) # use no interaction model

##### Resp mid and late periods
#### 2022 ####
str(dat.22) # 373 observations

## Analysis
periods.22<-dat.22%>%
  mutate(period=case_when(
    new.Date <= ymd("2022-07-06") ~ "early-mid",
    new.Date > ymd("2022-07-06") ~ "mid-late"))%>%
  mutate(period=factor(period,levels=c("early-mid","mid-late")))
str(periods.22)

## full combined seasonal ANOVAs
per.22<-lmer(sqrt(resp)~TRT.2022*period+(1|Plot/period),data=periods.22)
plot(per.22)
qqnorm(resid(per.22));qqline(resid(per.22))

anova(per.22,ddf="Kenward-Roger") # interaction not sig.
pairs(emmeans(per.22,pairwise~TRT.2022|period)) # both p < 0.0001

step(per.22,reduce.random=FALSE) # keep additive

per.22.ni<-lmer(sqrt(resp)~TRT.2022+period+(1|Plot/period),data=periods.22)
plot(per.22.ni)
qqnorm(resid(per.22.ni));qqline(resid(per.22.ni))
anova(per.22.ni,ddf="Kenward-Roger") # both significant
pairs(emmeans(per.22.ni,pairwise~TRT.2022)) # p < 0.0001 across trts

# check preferred mod
anova(per.22,per.22.ni) # prefer reduced model

#### 2023 - Figure 5a ####
str(dat.23) # 749 observations

## Analysis
periods.23<-dat.23%>%
  mutate(period=case_when(
    new.Date <= ymd("2023-07-06") ~ "early-mid",
    new.Date > ymd("2023-07-06") ~ "mid-late"))%>%
  mutate(period=factor(period,levels=c("early-mid","mid-late")))
str(periods.23)

period.avs.23<-periods.23%>%
  group_by(period,Plot,TRT.2022,TRT.2023)%>%
  dplyr::summarise(
    n.samp.dates= n(),
    season.av = mean(resp),
    season.sd = sd(resp),
    season.SE = season.sd/sqrt(n.samp.dates))

## full combined seasonal ANOVAs
season.23.avs.mod<-lmer(season.av ~ period*TRT.2022*TRT.2023+(1|Plot),data=period.avs.23)
plot(season.23.avs.mod)
qqnorm(resid(season.23.avs.mod));qqline(resid(season.23.avs.mod))
anova(season.23.avs.mod,ddf="Kenward-Roger") # three way interaction NOT significant
mod.pairs<-pairs(emmeans(season.23.avs.mod,~TRT.2022|TRT.2023*period))
resp.period.23.pairs<-as.data.frame(summary(mod.pairs))

step(season.23.avs.mod)
season.23.avs.mod.red<-lmer(season.av ~ period + TRT.2022 + (1|Plot), data=period.avs.23)
anova(season.23.avs.mod.red,ddf="Kenward-Roger")
plot(season.23.avs.mod.red)
qqnorm(resid(season.23.avs.mod.red));qqline(resid(season.23.avs.mod.red))

