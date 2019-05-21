###Phenological Variation in Ringed Salamanders###

#set working directory
setwd("C:/Users/Tom/Documents/GitRepos/AnnulatumPhenologicalVariation/")

#clear workspace
rm(list=ls())

#load necessary packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)
library(vegan)
library(cowplot)
library(nlme)
library(emmeans)
library(lme4)

#Make function to change dates
extrdate<-function(dat,name){
  dat$Year<-as.numeric(format(as.Date(as.character(dat[,name]),"%m/%d/%Y"),"%Y"))
  dat$Month<-as.numeric(format(as.Date(as.character(dat[,name]),"%m/%d/%Y"),"%m"))
  dat$Day<-as.numeric(format(as.Date(as.character(dat[,name]),"%m/%d/%Y"),"%d"))
  dat$Jdate<-as.numeric(format(as.Date(as.character(dat[,name]),"%m/%d/%Y"),"%j"))
  return(dat)
}

#load data
meta_size<-read.csv("Data/metamorph_size.csv",na.strings = c("-"))
larv_size <- read.csv("Data/larval_size.csv")
tank_key<-read.csv("Data/tank_key.csv",na.strings = c("-"))
larv_surv<-read.csv("Data/larval_survival.csv",na.strings = c("-"))
chla<-read.csv("Data/chla.csv")
zoop<-read.csv("Data/zooplankton.csv")

#Combine and extract relevant
meta_size<-merge(meta_size,tank_key,all=T,by="Tank")
larv_surv<-merge(larv_surv,tank_key,by="Tank")
chla<-merge(chla,tank_key,by.x="Sample",by.y="Tank")

#extract month, day, year information and add Julian date column
meta_size<-extrdate(meta_size,name="Date")
larv_size<-extrdate(larv_size,name="Date")
larv_surv<-extrdate(larv_surv,name="Date")
chla<-extrdate(chla,name="DateCollected")
zoop<-extrdate(zoop,name="Date_Sampled")

#Attach larval counts at end of experiment to rest of data to calculate survival
meta_size<-bind_rows(meta_size,larv_surv)

#Calculate quantity of metamorphs 
meta_size$Quantity<-ifelse(is.na(meta_size$Quantity),1,meta_size$Quantity)

# Annulatum Larval Results ------------------------------------------------

#Take the average, sd and cv of each tank for larval size
larvdat<-larv_size%>%
  group_by(Tank,Treatment,Jdate,Date)%>%
  summarise(meanHW=mean(HW,na.rm=T),sdHW=sd(HW,na.rm=T),cvHW=sdHW/meanHW)

#take the treatment average, sd of each quantity above to plot means +/- SD
larv_plotdat<-larvdat%>%
  group_by(Treatment,Jdate,Date)%>%
  summarise(mean_mean=mean(meanHW,na.rm=T),mean_sd=mean(sdHW,na.rm=T),mean_cv=mean(cvHW,na.rm=T),
            sd_mean=sd(meanHW,na.rm=T),sd_sd=sd(sdHW,na.rm=T),sd_cv=sd(cvHW,na.rm=T))

#generate error bars for both mean and CV of head width
lim.mean<-aes(ymax=larv_plotdat$mean_mean+larv_plotdat$sd_mean,ymin=larv_plotdat$mean_mean-larv_plotdat$sd_mean)
lim.cv<-aes(ymax=larv_plotdat$mean_cv+larv_plotdat$sd_cv,ymin=larv_plotdat$mean_cv-larv_plotdat$sd_cv)

#Plot mean larval head width
mean.pl<-ggplot(larv_plotdat,aes(Date,mean_mean,shape=Treatment,group=Treatment,color=Treatment,shape=Treatment))+
  geom_point(size=3)+geom_errorbar(lim.mean,width=0.75)+
  geom_line()+labs(y="Larval Head Width (mm)",x="")+
  scale_x_discrete(limits=c("3/21/2018","4/9/2018","4/26/2018"))+
  scale_color_manual(values=c("black","orange","skyblue"),breaks=c("one_date","three_dates","six_dates"),labels=c("Low","Medium","High"))+
  theme(legend.position=c(0.1,0.8),legend.text=element_text(size=10),legend.title=element_text(size=10))+
  scale_shape_manual(values=c(15,17,19),breaks=c("one_date","three_dates","six_dates"),labels=c("Low","Medium","High"))
mean.pl

#Plot CV of larval head width
cv.pl<-ggplot(larv_plotdat,aes(Date,mean_cv,shape=Treatment,group=Treatment,color=Treatment))+
  geom_point(size=3)+
  geom_errorbar(lim.cv,width=0.75)+
  geom_line()+
  scale_x_discrete(limits=c("3/21/2018","4/9/2018","4/26/2018"))+
  labs(y="CV of Larval HW",x="Julian Date")+
  theme(legend.position = "none")+
  scale_color_manual(values=c("black","orange","skyblue"))+
  scale_shape_manual(values=c(15,17,19),breaks=c("one_date","three_dates","six_dates"),labels=c("Low","Medium","High"))

#save plot as .tiff file
tiff("Results/Fig1.tiff",height=7,width=3.5,units="in",res=600,compression=c("lzw"))
plot_grid(mean.pl,cv.pl,ncol=1,labels=c("A","B"),align="hv")
dev.off()

#Analyze Growth Rates
growthmod<-lme(meanHW~Treatment*Jdate,random=~Jdate|Tank,data=larvdat,na.action=na.omit)
summary(growthmod)
Anova(growthmod,type =3, test.statistic="F")

cvmod<-lme(cvHW~Treatment+Jdate,random=~Jdate|Tank,data=larvdat,na.action=na.omit)
summary(cvmod)
Anova(cvmod,test.statistic="F")


# Annulatum Metamorph Results ---------------------------------------------

#Aggregate annulatum metamorphosis data
aman_dat<-meta_size%>%
  filter(Species=="AMAN",Treatment!="Control")%>%
  group_by(Treatment,Tank)%>%
  summarise(SVL_mean=mean(SVL*10,na.rm=T),Mass_mean=mean(Mass,na.rm=T),Day_mean=mean(Jdate,na.rm=T),
            SVL_sd=sd(SVL*10,na.rm=T),Mass_sd=sd(Mass,na.rm=T),Day_sd=sd(Jdate,na.rm=T),
            SVL_cv=SVL_sd/SVL_mean,Mass_cv=Mass_sd/Mass_mean,Day_cv=Day_sd/Day_mean,
            Quantity=sum(Quantity,na.rm=T))%>%
  complete(Tank)
aman_dat$Quantity<-ifelse(is.na(aman_dat$Quantity),0,aman_dat$Quantity)

#plot annulatum data by treatment
aa.svl.pl<-ggplot(aman_dat,aes(Treatment,SVL_mean))+
  geom_boxplot()+labs(y="SVL (cm)",x="")+
  scale_x_discrete(labels=c("","",""));aa.svl.pl
aa.mass.pl<-ggplot(aman_dat,aes(Treatment,Mass_mean))+
  geom_boxplot()+
  labs(y="Mass (g)",x="")+
  scale_x_discrete(labels=c("","",""));aa.mass.pl
aa.svlcv.pl<-ggplot(aman_dat,aes(Treatment,SVL_cv))+
  geom_boxplot()+
  labs(y="CV of SVL",x="")+
  scale_x_discrete(labels=c("","",""))+
  annotate(geom="text",x=c(1,2,3),y=0.12,label=c("A","B","AB"));aa.svlcv.pl
aa.masscv.pl<-ggplot(aman_dat,aes(Treatment,Mass_cv))+
  geom_boxplot()+
  labs(y="CV of Mass",x="")+
  scale_x_discrete(labels=c("","",""));aa.masscv.pl
aa.day.pl<-ggplot(aman_dat,aes(Treatment,Day_cv))+
  geom_boxplot()+
  labs(y="CV of Metamorphosis Date",x="Phenology Treatment")+
  scale_x_discrete(labels=c("Low","Medium","High"));aa.day.pl
aa.surv.pl<-ggplot(aman_dat,aes(Treatment,Quantity/36))+
  geom_boxplot()+lims(y=c(0,1))+
  labs(y="Percent Survival",x="Phenology Treatment")+
  scale_x_discrete(labels=c("Low","Medium","High"));aa.surv.pl

#save plot as .tiff file
tiff("Results/Fig2.tiff",res=600,height=10.5,width=7,units="in",compression=c("lzw"))
plot_grid(aa.svl.pl,aa.svlcv.pl,aa.mass.pl,aa.masscv.pl,aa.surv.pl,aa.day.pl,
          ncol=2,align="hv",labels=LETTERS[1:6])
dev.off()

#Annulatum analysis
svl.mod<-lm(SVL_mean~Treatment,data=aman_dat)
svlcv.mod<-lm(SVL_cv~Treatment,data=aman_dat)
mass.mod<-lm(Mass_mean~Treatment,data=aman_dat)
masscv.mod<-lm(Mass_cv~Treatment,data=aman_dat)
daycv.mod<-lm(Day_cv~Treatment,data=aman_dat)
surv.mod<-glmer(cbind(Quantity,36-Quantity)~Treatment+(1|Tank),data=aman_dat,binomial)
Anova(svl.mod)
Anova(svlcv.mod)
Anova(mass.mod)
Anova(masscv.mod)
Anova(daycv.mod)
Anova(surv.mod)
summary(svl.mod)
summary(svlcv.mod)
summary(mass.mod)
summary(masscv.mod)
summary(daycv.mod)
summary(surv.mod)

# Prey Results ------------------------------------------------------------

#Calculate prey survival and diversity
surv_wide<-meta_size%>%
  filter(Treatment!="Control")%>%
  group_by(Tank,Species)%>%
  summarise(Quantity=sum(Quantity,na.rm=T))%>%
  spread(key=Species,value = Quantity,fill = "0")
surv_wide<-sapply(as.data.frame(surv_wide),as.numeric)
AMSP<-surv_wide[,"AMMA"]+surv_wide[,"AMTE"]
surv_wide<-as.data.frame(cbind(surv_wide,AMSP))
surv_wide$Diversity<-diversity(surv_wide[,c("BUAM","PSFE","RASP","HYLA","AMSP")])
surv_wide$Richness<-apply(surv_wide[,c("BUAM","PSFE","RASP","HYLA","AMSP")],1,function(x){sum(x>0)})

#merge with larval and metamorph annulatum data
all_wide<-merge(surv_wide,tank_key,by="Tank")
all_wide<-merge(all_wide,aman_dat[,c('Tank','Treatment','SVL_mean','Day_mean')],by=c("Tank","Treatment"),all=T)
all_wide<-merge(all_wide,larvdat[larvdat$Jdate==80,c("Tank","meanHW")],by="Tank",all=T)

#calculate survival
all_wide$TotalSurvival<-rowSums(all_wide[,c("AMSP","BUAM","PSFE","HYLA","RASP")],na.rm=T)
all_wide$TotalInitial<-rowSums(all_wide[,c("InitialAMMA","InitialBUAM","InitialPSFE","InitialHYLA","InitialRASP")],na.rm=T)
all_wide$PercentSurvival<-all_wide$TotalSurvival/all_wide$TotalInitial

#Analyze total survival of amphiban prey
totsurv.mod<-glmer(cbind(TotalSurvival,TotalInitial-TotalSurvival)~Treatment*meanHW+AMAN+(1|Tank),data=all_wide,family=binomial)
Anova(totsurv.mod)
pairs(emmeans::emmeans(totsurv.mod,specs="Treatment"))

#re-run model ignoring potential outlier HW value
summary(update(totsurv.mod,subset=c(meanHW>5)))
Anova(update(totsurv.mod,subset=c(meanHW>5)))

#Analyze survival and diversity with predator size
div.size<-lm(Diversity~Treatment+meanHW,data=all_wide)
Anova(div.size)
summary(div.size)

#re-run model ignoring potential outlier HW value
summary(update(div.size,subset=c(meanHW>5)))

#Plot overall survival and diversity relationships
surv.trt.pl<-ggplot(all_wide,aes(Treatment,PercentSurvival))+
  geom_boxplot()+
  labs(y="Total Prey Percent Survival",x="")+
  scale_x_discrete(labels=c("Low","Medium","High"),breaks=c("1 date","3 dates","6 dates"))
surv.hw.pl<-ggplot(all_wide,aes(meanHW,PercentSurvival,group=Treatment))+
  geom_smooth(method = "glm", method.args = list(family = "binomial"),aes(color=Treatment),size=2,se=F)+
  geom_point(aes(shape=Treatment,color=Treatment),size=3)+
  labs(x="Mean Larval HW (mm)",y="")+
  scale_shape_manual(name="Treatment",label=c("Low","Medium","High"),values=c(15,17,19))+
  scale_color_manual(name="Treatment",label=c("Low","Medium","High"),values=c("black","orange","skyblue"))+
  theme(legend.position=c(0.7,0.8),legend.text=element_text(size=10),legend.title=element_text(size=10))
surv.aman.pl<-ggplot(all_wide,aes(AMAN,PercentSurvival))+
  geom_point(aes(shape=Treatment,color=Treatment),size=3)+
  geom_smooth(method = "glm", method.args = list(family = "binomial"),color="black",size=2,se=F)+
  labs(x=expression("Number of "~italic(A.~annulatum)),y="")+
  scale_shape_manual(name="Treatment",label=c("Low","Medium","High"),values=c(15,17,19))+
  scale_color_manual(values=c("black","orange","skyblue"),name="Treatment",label=c("Low","Medium","High"))+
  theme(legend.position="none")
div.trt.pl<-ggplot(all_wide,aes(Treatment,Diversity))+
  geom_boxplot()+
  labs(y="Shannon Diversity",x="Treatment")+
  scale_x_discrete(labels=c("Low","Medium","High"),breaks=c("1 date","3 dates","6 dates"))+
  lims(y=c(0,1.5))
div.hw.pl<-ggplot(all_wide,aes(meanHW,Diversity))+
  geom_point(aes(shape=Treatment,color=Treatment),size=3)+
  geom_smooth(method = "lm",color="black",size=2,se=F)+
  scale_shape_manual(values=c(15,17,19))+
  scale_color_manual(values=c("black","orange","skyblue"))+
  labs(x="Mean Larval HW (mm)",y="")+
  lims(y=c(0,1.5))+theme(legend.position="none")

#generate survival summary data for each species for plotting
surv_dat<-all_wide%>%
  mutate(PerSurv_BUAM=BUAM/InitialBUAM,
         PerSurv_PSFE=PSFE/InitialPSFE,
         PerSurv_HYLA=HYLA/InitialHYLA,
         PerSurv_RASP=RASP/InitialRASP,
         PerSurv_AMSP=AMSP/InitialAMMA)%>%
  dplyr::select(Tank,Treatment,PerSurv_BUAM:PerSurv_AMSP)%>%
  gather(Species,Survival,PerSurv_BUAM:PerSurv_AMSP)

#generate color palette
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
preysurv<-ggplot(surv_dat,aes(Treatment,Survival,fill=Species))+
  geom_boxplot()+
  labs(y="Percent Survival",x="Treatment")+
  lims(y=c(0,0.75))+
  scale_fill_manual(values=cbbPalette,name="",labels=c("AMSP","ANAM","HYLA","PSFE","RASP"))+
  theme(legend.position = c(x=0.01,y=0.9),legend.direction = "horizontal",legend.text=element_text(size=10),legend.title=element_text(size=10))+
  scale_x_discrete(labels=c("Low","Medium","High"),limits=c("1 date","3 dates","6 dates"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

#save plot as .tiff file
tiff("Results/Fig3.tiff",height=7,width=9,units="in",compression=c("lzw"),res=600)
plot_grid(surv.trt.pl,surv.hw.pl,surv.aman.pl,preysurv,div.trt.pl,div.hw.pl,
          align="hv",labels=c(LETTERS[1:6]),ncol=3,hjust=-2)
dev.off()

#Ambystoma analysis
amsp_dat<-meta_size%>% 
  filter(Species%in%c("AMMA","AMTE") & Treatment !="Control")%>%
  group_by(Treatment,Tank)%>%
  summarise(SVL_mean=mean(SVL*10,na.rm=T),Mass_mean=mean(Mass,na.rm=T),Day_mean=mean(Jdate,na.rm=T),
            SVL_sd=sd(SVL*10,na.rm=T),Mass_sd=sd(Mass,na.rm=T),Day_sd=sd(Jdate,na.rm=T),
            SVL_cv=SVL_sd/SVL_mean,Mass_cv=Mass_sd/Mass_mean,Day_cv=Day_sd/Day_mean,
            Quantity=sum(Quantity,na.rm=T))%>%
  complete(Tank)

#combine with annulatum data
amsp_dat<-merge(amsp_dat,all_wide[,c("Tank","AMAN","meanHW")],by="Tank")

#Analyze Ambystoma sp. responses
amsp.svl<-lm(SVL_mean~Treatment,data=amsp_dat)
amsp.mass<-lm(Mass_mean~Treatment,data=amsp_dat)
amsp.day<-lm(Day_mean~Treatment,data=amsp_dat)
amsp.surv<-glmer(cbind(AMSP,InitialAMMA-AMSP)~Treatment*meanHW+AMAN+(1|Tank),data=all_wide,family=binomial)
rbind(Anova(amsp.svl),Anova(amsp.mass),Anova(amsp.day))
Anova(amsp.svl)
Anova(amsp.mass)
Anova(amsp.day)
Anova(amsp.surv)

summary(amsp.svl)
summary(amsp.mass)
summary(amsp.day)
summary(amsp.surv)

#Plot Ambystoma responses in relation to phenology treatments
am.svl.pl<-ggplot(amsp_dat,aes(Treatment,SVL_mean))+
  geom_boxplot()+labs(y="SVL (cm)",x="")+
  scale_x_discrete(breaks=c("1 date","3 dates","6 dates"),labels=c("Low","Medium","High"));am.svl.pl
am.mass.pl<-ggplot(amsp_dat,aes(Treatment,Mass_mean))+
  geom_boxplot()+
  labs(y="Mass (g)",x="")+
  scale_x_discrete(breaks=c("1 date","3 dates","6 dates"),labels=c("Low","Medium","High"));am.mass.pl
am.day.pl<-ggplot(amsp_dat,aes(Treatment,Day_mean))+
  geom_boxplot()+
  labs(y="Date of Metamorphosis",x="")+
  scale_x_discrete(breaks=c("1 date","3 dates","6 dates"),labels=c("Low","Medium","High"));am.day.pl
am.surv.pl<-ggplot(amsp_dat,aes(Treatment,Quantity/45))+
  geom_boxplot()+
  labs(y="Percent Survival")+
  scale_x_discrete(breaks=c("1 date","3 dates","6 dates"),labels=c("Low","Medium","High"));am.surv.pl

#save plot as .tiff file
tiff("Results/FigA1.tiff",res=600,height=7,width=7,units="in",compression=c("lzw"))
plot_grid(am.svl.pl,am.mass.pl,am.day.pl,am.surv.pl,align="hv",labels=LETTERS[1:4])
dev.off()

#Plot Ambystoma responses in relation to annulatum density and size
fig4a<-ggplot(amsp_dat,aes(AMAN,Quantity/45))+
  geom_point(aes(shape=Treatment,color=Treatment),size=3)+lims(y=c(0,1))+
  geom_smooth(method="glm",method.args = list(family = "binomial"),se=F,color="black")+
  labs(y="Percent Survival",x=expression("Number of surviving"~italic(A.~annulatum)))+
  scale_shape_manual(values=c(15,17,19),breaks=c("1 date","3 dates","6 dates"),labels=c("Low","Medium","High"))+
  scale_color_manual(values=c("black","orange","skyblue"),name="Treatment",label=c("Low","Medium","High"),breaks=c("1 date","3 dates","6 dates"))+
  theme(legend.position=c(0.7,0.8),legend.text=element_text(size=10),legend.title=element_text(size=10))
fig4b<-ggplot(amsp_dat,aes(meanHW,Quantity/45))+
  geom_point(aes(shape=Treatment,color=Treatment),size=3)+lims(y=c(0,1))+
  geom_smooth(method="glm",method.args = list(family = "binomial"),se=F,aes(color=Treatment))+
  labs(y="",x="Mean Larval HW (mm)")+
  scale_shape_manual(values=c(15,17,19),breaks=c("1 date","3 dates","6 dates"),labels=c("Low","Medium","High"))+
  scale_color_manual(values=c("black","orange","skyblue"),name="Treatment",label=c("Low","Medium","High"),breaks=c("1 date","3 dates","6 dates"))+
  theme(legend.position="none")

#save plot as .tiff file
tiff("Results/Fig4.tiff",res=600,height=3.5,width=7,units="in",compression=c("lzw"))
plot_grid(fig4a,fig4b,ncol=2,labels=c("A","B"))
dev.off()

#RASP data aggregation
rasp_dat<-meta_size%>% 
  filter(Species=="RASP")%>%
  group_by(Treatment,Tank)%>%
  summarise(SVL_mean=mean(SVL*10,na.rm=T),Mass_mean=mean(Mass,na.rm=T),Day_mean=mean(Jdate,na.rm=T),
            SVL_sd=sd(SVL*10,na.rm=T),Mass_sd=sd(Mass,na.rm=T),Day_sd=sd(Jdate,na.rm=T),
            SVL_cv=SVL_sd/SVL_mean,Mass_cv=Mass_sd/Mass_mean,Day_cv=Day_sd/Day_mean,
            Quantity=sum(Quantity,na.rm=T))%>%
  complete(Tank)

#merge with annulatum data
rasp_dat<-merge(rasp_dat,all_wide[,c("Tank","AMAN","meanHW")],by="Tank")

#Plot RASP data
rs.svl.pl<-ggplot(rasp_dat,aes(Treatment,SVL_mean))+
  geom_boxplot()+
  labs(y="SVL (mm)",x="")+
  scale_x_discrete(breaks=c("1 date","3 dates","6 dates"),labels=c("Low","Medium","High"));rs.svl.pl
rs.mass.pl<-ggplot(rasp_dat,aes(Treatment,Mass_mean))+
  geom_boxplot()+
  labs(y="Mass (g)",x="")+
  scale_x_discrete(breaks=c("1 date","3 dates","6 dates"),labels=c("Low","Medium","High"));rs.mass.pl
rs.day.pl<-ggplot(rasp_dat,aes(Treatment,Day_mean))+
  geom_boxplot()+
  labs(y="Date of Metamorphosis",x="Treatment")+
  scale_x_discrete(breaks=c("1 date","3 dates","6 dates"),labels=c("Low","Medium","High"));rs.day.pl
rs.surv.pl<-ggplot(rasp_dat,aes(Treatment,Quantity/200))+
  geom_boxplot()+
  labs(y="Percent Survival")+
  scale_x_discrete(breaks=c("1 date","3 dates","6 dates"),labels=c("Low","Medium","High"));rs.surv.pl

#save plot as .tiff file
tiff("Results/FigA2.tiff",res=600,height=7,width=7,units="in",compression=c("lzw"))
plot_grid(rs.svl.pl,rs.mass.pl,rs.day.pl,rs.surv.pl,align="hv",labels=LETTERS[1:4])
dev.off()

#Plot RASP survival in relation to annulatum body size
rasp.aman.pl<-ggplot(all_wide,aes(meanHW,RASP/200))+
  geom_point(aes(shape=Treatment,color=Treatment),size=3)+
  geom_smooth(method="glm",method.args = list(family = "binomial"),se=F,color="black")+
  labs(y="Percent Survival",x="Mean Larval HW (mm)")+
  scale_shape_manual(values=c(15,17,19),breaks=c("1 date","3 dates","6 dates"),labels=c("Low","Medium","High"))+
  scale_color_manual(values=c("black","orange","skyblue"),name="Treatment",label=c("Low","Medium","High"),breaks=c("1 date","3 dates","6 dates"))+
  theme(legend.position=c(0.7,0.8),legend.text=element_text(size=10),legend.title=element_text(size=10))

#save plot as .tiff file
tiff("Results/FigA3.tiff",res=600,height=3.5,width=3.5,units="in",compression=c("lzw"))
rasp.aman.pl
dev.off()

#Analyze RASP responses
rasp.svl<-lm(SVL_mean~Treatment,data=rasp_dat)
rasp.mass<-lm(Mass_mean~Treatment,data=rasp_dat)
rasp.smi<-lm(SMI_mean~Treatment,data=rasp_dat)
rasp.day<-lm(Day_mean~Treatment,data=rasp_dat)
rasp.surv<-glmer(cbind(RASP,InitialRASP-RASP)~Treatment+meanHW+(1|Tank),data=all_wide,family=binomial)
Anova(rasp.svl)
Anova(rasp.mass)
Anova(rasp.day)
Anova(rasp.surv)
#re-run survival to see if HW outlier matters
Anova(update(rasp.surv,subset=c(meanHW>5)))

summary(rasp.svl)
summary(rasp.mass)
summary(rasp.day)
summary(rasp.surv)
#re-run survival to see if HW outlier matters
summary(update(rasp.surv,subset=c(meanHW>5)))

#analyze PSFE survival
psfe.surv<-glmer(cbind(PSFE,InitialPSFE-PSFE)~Treatment+(1|Tank),data=all_wide,family=binomial)
Anova(psfe.surv)

# Plankton Results --------------------------------------------------------

#Phytoplankton plots and analysis
chla.pl<-ggplot(filter(chla),aes(Treatment,log(TotChla_mgL)))+
  geom_boxplot()+
  scale_x_discrete(limits=c("1 date","3 dates","6 dates"),labels=c("Low","Medium","High"))+
  labs(x="",y="log-Chl a (mg/L)");may_chla.pl

#combine with anuran and salamander data
surv_wide$TotAnurans<-rowSums(surv_wide[,c("BUAM","HYLA","PSFE","RASP")])
chla<-merge(surv_wide,chla,by.x=c("Tank"),by.y=c("Sample"),all=T)

#models
chla.mod<-lm(log(TotChla_mgL)~Treatment+AMAN+TotAnurans,data=chla,subset=c(DateCollected=="5/4/2018"))
Anova(chla.mod)
summary(chla.mod)

#combine zoolankton data with salamander data
zoop<-merge(surv_wide,zoop,by="Tank",all=T)

totzoop.pl<-ggplot(zoop,aes(Treatment,log1p(Total.L)))+
  geom_boxplot()+
  scale_x_discrete(limits=c("1 date","3 dates","6 dates"),labels=c("Low","Medium","High"))+
  labs(x="Treatment",y="log-Total Zooplankton (per L)");may_totzoop.pl

#save plot as .tiff file
tiff("Results/FigA5.tiff",res=600,height=7,width=3.5,compression=c("lzw"),units="in")
plot_grid(chla.pl,totzoop.pl,ncol=1,labels=c("A","B"))
dev.off()

#Plankton analysis
tot.logmod<-lm(log1p(Total.L)~Treatment,data=zoop)
Anova(tot.logmod)
summary(tot.logmod)

tot.poismod<-glmer(Total.L~Treatment+(1|Tank),data=zoop,family="poisson")
Anova(tot.poismod)
tot.negbinmod<-glmer.nb(Total.L~Treatment+(1|Tank),data=zoop)
Anova(tot.negbinmod)

# Temporal Overlap --------------------------------------------------------

#Temporal overlap metrics
par(mfrow=c(3,2),mar=c(3,4,0,0))
hist(aman_dat$Day_mean-64)
rasp.ol<-mean(aman_dat$Day_mean-64,na.rm=T)

hist(aman_dat$Day_mean-70)
psfe.ol<-mean(aman_dat$Day_mean-70,na.rm=T)

hist(aman_dat$Day_mean-93)
amsp.ol<-mean(aman_dat$Day_mean-93,na.rm=T)

hist(aman_dat$Day_mean-97)
anam.ol<-mean(aman_dat$Day_mean-97,na.rm=T)

hist(aman_dat$Day_mean-129)
hyla.ol<-mean(aman_dat$Day_mean-129,na.rm=T)
