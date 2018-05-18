##----------------------------------
## BLW 2018 glmer
## Zhengpei Ye 20.03.2018 Innsbruck
##-----------------------------------

setwd("C:/Users/c7701076/1_Zhengpei_YE/R_analysis/BLW 2018")

library(MASS)
library(effects)
library(nlme)
library(lme4)
library(reshape2)

D1<-read.csv(file="Sample 2012 BLW 2018.csv",header=T)

str(D1)
names(D1)
unique(paste(D1$field,D1$plot,D1$fertilization))


D1$hyperparasitoid<-ifelse(rowSums(D1[,c(22:27)])>=1,1,0)


te<-table(D1$Hden,D1$Rins)
chisq.test(te)

D2<-subset(D1,Save==1)
ts<-table(D2$Hden,D2$Rins)
chisq.test(ts)

D3<-subset(D1,Mdir==1)
tm<-table(D3$Hden,D3$Rins)
chisq.test(tm)

D4<-subset(D1,Rpad==1)
tr<-table(D4$Hden,D4$Rins)
chisq.test(tr)


Data<-D1[,c(1:13,28)]

names(Data)
str(Data)

Data[,1]<-as.character(Data[,1])

for(i in 2:6){
  
  Data[,i]<-as.factor(Data[,i])
  
}

Data$ParaSt<-ifelse(Data$AphidMummy=="aph",ifelse(Data$parasitoid==1,"1P",0),"2M")
Data$ParaSt<-as.factor(Data$ParaSt)

Data$Aphid<-ifelse(Data$Mdir==0,ifelse(Data$Save==0,"Rp","Sa"),"Md")
Data$Aphid<-as.factor(Data$Aphid)


names(Data)

dim(Data[Data$AphidMummy=="aph",])
dim(Data[Data$AphidMummy=="mum",])

dim(Data[Data$AphidMummy=="aph"&Data$Mdir==1,])
dim(Data[Data$AphidMummy=="aph"&Data$Save==1,])
dim(Data[Data$AphidMummy=="aph"&Data$Rpad==1,])

dim(Data[Data$AphidMummy=="mum"&Data$Mdir==1,])
dim(Data[Data$AphidMummy=="mum"&Data$Save==1,])
dim(Data[Data$AphidMummy=="mum"&Data$Rpad==1,])

dim(Data[Data$AphidMummy=="aph"&Data$EsymL=="H0",])/dim(Data[Data$AphidMummy=="aph",])
dim(Data[Data$AphidMummy=="aph"&Data$EsymL=="R0",])/dim(Data[Data$AphidMummy=="aph",])
dim(Data[Data$AphidMummy=="aph"&Data$EsymL=="RH",])/dim(Data[Data$AphidMummy=="aph",])

dim(Data[Data$AphidMummy=="mum"&Data$EsymL=="H0",])/dim(Data[Data$AphidMummy=="mum",])
dim(Data[Data$AphidMummy=="mum"&Data$EsymL=="R0",])/dim(Data[Data$AphidMummy=="mum",])
dim(Data[Data$AphidMummy=="mum"&Data$EsymL=="RH",])/dim(Data[Data$AphidMummy=="mum",])

dim(Data[Data$AphidMummy=="aph"&Data$Mdir==1&Data$endosymbiont==1,])/dim(Data[Data$AphidMummy=="aph"&Data$Mdir==1,])
dim(Data[Data$AphidMummy=="aph"&Data$Rpad==1&Data$endosymbiont==1,])/dim(Data[Data$AphidMummy=="aph"&Data$Rpad==1,])
dim(Data[Data$AphidMummy=="aph"&Data$Save==1&Data$endosymbiont==1,])/dim(Data[Data$AphidMummy=="aph"&Data$Save==1,])

dim(Data[Data$AphidMummy=="mum"&Data$Mdir==1&Data$endosymbiont==1,])/dim(Data[Data$AphidMummy=="mum"&Data$Mdir==1,])
dim(Data[Data$AphidMummy=="mum"&Data$Rpad==1&Data$endosymbiont==1,])/dim(Data[Data$AphidMummy=="mum"&Data$Rpad==1,])
dim(Data[Data$AphidMummy=="mum"&Data$Save==1&Data$endosymbiont==1,])/dim(Data[Data$AphidMummy=="mum"&Data$Save==1,])

dim(Data[Data$Mdir==1&Data$endosymbiont==1,])/dim(Data[Data$Mdir==1,])
dim(Data[Data$Rpad==1&Data$endosymbiont==1,])/dim(Data[Data$Rpad==1,])
dim(Data[Data$Save==1&Data$endosymbiont==1,])/dim(Data[Data$Save==1,])


m1<-glmer(endosymbiont~ParaSt+fertilization+Aphid+(1|field/plot/tube), Data, family="binomial",
          control=glmerControl(optCtrl=list(maxfun=2e5), optimizer="bobyqa"))

m2<-glmer(endosymbiont~ParaSt*fertilization+Aphid+(1|field/plot/tube), Data, family="binomial",
          control=glmerControl(optCtrl=list(maxfun=2e5), optimizer="bobyqa"))
m3<-glmer(endosymbiont~ParaSt*fertilization+Aphid*fertilization+(1|field/plot/tube), Data, family="binomial",
          control=glmerControl(optCtrl=list(maxfun=2e5), optimizer="bobyqa"))

m4<-glmer(endosymbiont~ParaSt*Aphid+fertilization+(1|field/plot/tube), Data, family="binomial",
          control=glmerControl(optCtrl=list(maxfun=2e5), optimizer="bobyqa"))

m5<-glmer(endosymbiont~ParaSt*Aphid+fertilization*Aphid+(1|field/plot/tube), Data, family="binomial",
          control=glmerControl(optCtrl=list(maxfun=2e5), optimizer="bobyqa"))

m6<-glmer(endosymbiont~ParaSt*fertilization+ParaSt*Aphid+(1|field/plot/tube), Data, family="binomial",
          control=glmerControl(optCtrl=list(maxfun=2e5), optimizer="bobyqa"))

m7<-glmer(endosymbiont~ParaSt+fertilization*Aphid+(1|field/plot/tube), Data, family="binomial",
          control=glmerControl(optCtrl=list(maxfun=2e5), optimizer="bobyqa"))

m8<-glmer(endosymbiont~ParaSt*fertilization*Aphid+(1|field/plot/tube), Data, family="binomial",
          control=glmerControl(optCtrl=list(maxfun=2e5), optimizer="bobyqa"))

AIC(m1,m2,m3,m4,m5,m6,m7,m8)

drop1(m1,test="Chi")
drop1(m4,test="Chi")

ms1<-glmer(endosymbiont~ParaSt+Aphid+(1|field/plot/tube), Data, family="binomial",
          control=glmerControl(optCtrl=list(maxfun=2e5), optimizer="bobyqa"))
ms2<-glmer(endosymbiont~ParaSt*Aphid+(1|field/plot/tube), Data, family="binomial",
          control=glmerControl(optCtrl=list(maxfun=2e5), optimizer="bobyqa"))

AIC(ms1,ms2,m1,m4)

drop1(ms1,test="Chi")

plot(ms2)
library(arm)
x<-predict(ms2)
y<-resid(ms2)
binnedplot(x,y)

plot(allEffects(ms2))
summary(ms2)
anova(ms2)

Data$IntPsA<-with(Data, interaction(Aphid, ParaSt))

ms2.int<-glmer(endosymbiont~IntPsA+(1|field/plot/tube), Data, family="binomial",
           control=glmerControl(optCtrl=list(maxfun=2e5), optimizer="bobyqa"))

plot(allEffects(ms2.int))

library(multcomp)
PcomPsA<-glht(ms2.int, linfct=mcp(IntPsA="Tukey"))
summary(PcomPsA, test=adjusted("fdr"))


e<-as.data.frame(effect("ParaSt:Aphid",ms2) )
e75<-as.data.frame(effect("ParaSt:Aphid",ms2,confidence.level=0.75))
e$lower75<-e75$lower
e$upper75<-e75$upper
e$ParaSt<-factor(e$ParaSt,labels=c("living unparasitised aphids",
                                   "living parasitised aphids",
                                   "mummified aphids"))

e$sig<-c(" a"," ab"," c"," d"," cd"," cd"," b"," b"," cd")

library(ggplot2)
fp1<-ggplot(mapping=aes(middle=fit*100,
                        ymax=upper*100,ymin=lower*100,
                        upper=upper75*100,lower=lower75*100))

fp2<-fp1+aes(label=sig, y=upper75*100)

fp3<-fp2%+%e+aes(x=Aphid, fill=ParaSt)

fp4<-fp3+
  geom_boxplot(stat="identity", colour="black", position=position_dodge(0.93))+
  geom_text(aes(y=upper75*100+1.5),position=position_dodge(0.93), size=8, colour="black", vjust=0, hjust=0)+
  labs(y="Endosymbiont infection (%)", fill="Parasitism types")+
  scale_y_continuous(limits=c(0,105),expand=c(0,0))+
  scale_fill_manual(values=c("#C6DBEF", "#4292C6","#08306B"))+
  scale_x_discrete(labels=c(expression(italic("M. dirhodum")),
                            expression(italic("R. padi")),
                            expression(italic("S. avenae"))))


fp4+theme_bw()+theme(axis.text.y=element_text(size=18,colour="black"),                     
                     axis.title.y=element_text(size=20),
                     axis.title.x=element_blank(),
                     axis.text.x=element_text(size=20,colour="black"),
                     legend.justification=c(1,1),
                     legend.position=c(0.6,0.99),
                     legend.title=element_text(size=18),
                     legend.text=element_text(size=16))



## --------------------------------------------------------------------------------

names(Data)
Data$EsymL<-ifelse(Data$Hden==1,ifelse(Data$Rins==1,"RH","H0"),ifelse(Data$Rins==1,"R0",0))
Data$EsymL<-as.factor(Data$EsymL)

Dp<-subset(Data,AphidMummy=="aph")

names(Dp)
Dp<-subset(Data,AphidMummy=="mum")
Da<-Dp[,c(3,4,16)]
dim(Da)
library(reshape2)
names(Da)
Da$counter<-1
str(Da)
m1<-melt(Da)

dcast(m1, field+plot~variable, sum)

nrow(Dp)
nrow(Dp[Dp$EsymL=="RH",c(1,13)])
# Dp<-subset(Dp,EsymL!="RH")

mp1<-glmer(parasitoid~EsymL+fertilization+Aphid+(1|field), Dp, family="binomial")

mp2<-glmer(parasitoid~EsymL*fertilization+Aphid+(1|field), Dp, family="binomial")
mp3<-glmer(parasitoid~EsymL*fertilization+Aphid*fertilization+(1|field), Dp, family="binomial")

mp4<-glmer(parasitoid~EsymL*Aphid+fertilization+(1|field), Dp, family="binomial")
mp5<-glmer(parasitoid~EsymL*Aphid+fertilization*Aphid+(1|field), Dp, family="binomial")

mp6<-glmer(parasitoid~EsymL*fertilization+EsymL*Aphid+(1|field), Dp, family="binomial")

mp7<-glmer(parasitoid~EsymL+fertilization*Aphid+(1|field), Dp, family="binomial")

mp8<-glmer(parasitoid~EsymL*fertilization*Aphid+(1|field), Dp, family="binomial")

AIC(mp1,mp2,mp3,mp4,mp5,mp6,mp7,mp8)

drop1(mp1,test="Chi")
drop1(mp7,test="Chi")

plot(allEffects(mp2))

mps1<-glmer(parasitoid~fertilization+Aphid+(1|field), Dp, family="binomial")
mps2<-glmer(parasitoid~fertilization*Aphid+(1|field), Dp, family="binomial")

AIC(mps1,mps2,mp1,mp7)

drop1(mps1,test="Chi")
drop1(mps2,test="Chi")

mpss1<-glmer(parasitoid~Aphid+(1|field), Dp, family="binomial")
mpss2<-glmer(parasitoid~fertilization+(1|field), Dp, family="binomial")
mpss3<-glmer(parasitoid~1+(1|field), Dp, family="binomial")

AIC(mps2,mp7,mpss1,mpss2,mpss3)

plot(mps2)
library(arm)
x<-predict(mps2)
y<-resid(mps2)
binnedplot(x,y)

plot(allEffects(mpss1))
summary(mps2)
anova(mps2)

Dp$IntfA<-with(Dp, interaction(Aphid, fertilization))

summary(mps2)

mps2.int<-glmer(parasitoid~IntfA+(1|field), Dp, family="binomial")

plot(allEffects(mps2.int))

library(multcomp)
PcomfA<-glht(mps2.int, linfct=mcp(IntfA="Tukey"))
summary(PcomfA, test=adjusted("fdr"))

library(multcomp)
PcomA<-glht(mpss1, linfct=mcp(Aphid="Tukey"))
summary(PcomA, test=adjusted("fdr"))


## --------------------------------------------------------------------------------
names(Data)
dim(Data[Data$hyperparasitoid==1,])/dim(Data)

Dh<-subset(Data, parasitoid==1)
dim(Dh)

mh1<-glmer(hyperparasitoid~EsymL+fertilization+Aphid+(1|field), Dh, family="binomial")

mh2<-glmer(hyperparasitoid~EsymL*fertilization+Aphid+(1|field), Dh, family="binomial")
mh3<-glmer(hyperparasitoid~EsymL*fertilization+Aphid*fertilization+(1|field), Dh, family="binomial")

mh4<-glmer(hyperparasitoid~EsymL*Aphid+fertilization+(1|field), Dh, family="binomial")
mh5<-glmer(hyperparasitoid~EsymL*Aphid+fertilization*Aphid+(1|field), Dh, family="binomial")

mh6<-glmer(hyperparasitoid~EsymL*fertilization+EsymL*Aphid+(1|field), Dh, family="binomial")

mh7<-glmer(hyperparasitoid~EsymL+fertilization*Aphid+(1|field), Dh, family="binomial")

mh8<-glmer(hyperparasitoid~EsymL*fertilization*Aphid+(1|field), Dh, family="binomial")

AIC(mh1,mh2,mh3,mh4,mh5,mh6,mh7,mh8)

drop1(mh1,test="Chi")
drop1(mh7,test="Chi")

mhs1<-glmer(hyperparasitoid~EsymL+Aphid+(1|field), Dh, family="binomial")
mhs2<-glmer(hyperparasitoid~EsymL*Aphid+(1|field), Dh, family="binomial")

mhs3<-glmer(hyperparasitoid~EsymL+fertilization+(1|field), Dh, family="binomial")
mhs4<-glmer(hyperparasitoid~EsymL*fertilization+(1|field), Dh, family="binomial")

AIC(mhs1,mhs2,mhs3,mhs4)

drop1(mhs1,test="Chi")
drop1(mhs2,test="Chi")
drop1(mhs3,test="Chi")
drop1(mhs4,test="Chi")

mhss1<-glmer(hyperparasitoid~EsymL+(1|field), Dh, family="binomial")
mhss2<-glmer(hyperparasitoid~1+(1|field), Dh, family="binomial")

AIC(mhss1,mhss2)

AIC(mh7,mhs1,mhss1)

plot(allEffects(mhss1))

summary(mhss1)
anova(mhss1)

library(multcomp)
Pcomh<-glht(mhss1, linfct=mcp(EsymL="Tukey"))
summary(Pcomh, test=adjusted("fdr"))

plot(mhss1)
library(arm)
x<-predict(mhss1)
y<-resid(mhss1)
binnedplot(x,y)

h<-as.data.frame(effect("EsymL",mhss1) )
h75<-as.data.frame(effect("EsymL",mhss1,confidence.level=0.75))
h$lower75<-h75$lower
h$upper75<-h75$upper

library(ggplot2)
h$sig<-c(" a"," ab"," b", " ab")
fp1<-ggplot(mapping=aes(middle=fit*100,
                        ymax=upper*100,ymin=lower*100,
                        upper=upper75*100,lower=lower75*100))

fp2<-fp1+aes(label=sig, y=upper75*100)

fp3<-fp2%+%h+aes(x=EsymL)

fp4<-fp3+
  geom_boxplot(stat="identity", colour="black", fill="#1F78B4")+
  geom_text(aes(y=upper75*100+1), size=8, colour="black", vjust=0, hjust=0)+
  labs(y="Hyperparasitism rate (%)")+
  scale_y_continuous(limits=c(5,70),expand=c(0,0))+
  scale_x_discrete(labels=c("uninfected",
                            expression(paste(italic("H. defensa"),"-infected",sep="")),
                            expression(paste(italic("R. insecticola"),"-infected",sep="")),
                            "superinfected"))


fp4+theme_bw()+theme(axis.text.y=element_text(size=18,colour="black"),                     
                     axis.title.y=element_text(size=20),
                     axis.title.x=element_blank(),
                     axis.text.x=element_text(size=20,colour="black"))


