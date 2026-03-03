##############Lizard_traits fiter############
###########Zhong-wen Jiang 2021/12/28########
library(car) #for the leveneTest
library(nlme)
library(emmeans)
library(FSA)
library(effsize)
library(lme4)
library(lmerTest)


#####Morphology/Behaviour/Physiology####
##Body mass
#Read in data (body mass was measured when measuring CTmax, Tsel and food passage time)
BM<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/BM_all.csv")

mod_BMSVL<-lm(BM~SVL,data = BM)
summary(mod_BMSVL)#R²=0.74,p<0.001

####SVL as the covariates and the effect of elevation on body mass is nonsignificant
####morphological data of lizard (Phrynocephalus vlangalii) from literatures about different elevation 
BM_litera<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/BM_litera.csv")
BM_litera$Elevation<-as.factor(BM_litera$Elevation)
leveneTest(dat=BM_litera,BM~Elevation)
mod_litera<-lm(BM~Elevation+SVL,data = BM_litera)###Elevation:t=0.34,p=0.745
summary(mod_litera)
shapiro.test(x=residuals(object=mod_litera)) ###p=0.300

#SVL>5 at sexual maturity both (Delingha: elevation 2910m, Daotanghe: 3367m; Maduo: 4257m)
#Li J, Zhou R, Liu N. Life-history variation among three populations of the toad-headed lizard Phrynocephalus vlangalii along an elevation gradient on the northeastern Tibetan Plateau[J]. The Herpetological Journal, 2014, 24(1): 17-23.
dat_BM<-BM[BM$SVL>5,]
table(dat_BM$Elevation)


#Run Levene's test (check the homogeneity of variances)
leveneTest(dat=dat_BM,BM~Elevation) #p<0.001, there is significant difference between the group variances in Elevation
#Fit Linear Model Using Generalized Least Squares due to unequal variance; SVl as the covariates
mod_BM<-gls(BM~SVL+Elevation,weights=varIdent(form=~1|Elevation),data=dat_BM)
summary(mod_BM)
anova(mod_BM)#F=1.756,p=0.177
#Run Shapiro-Wilk test (check the normality of residuals)
shapiro.test(x=residuals(object=mod_BM)) #p=0.004, reject the null hypothesis that the residuals are normally distributed
hist(residuals(object=mod_BM))
qqnorm(residuals(object=mod_BM))

emmeans(mod_BM, list(pairwise~Elevation), adjust = "tukey")

#Elevation emmean     SE   df lower.CL upper.CL
#2600m       6.42 0.180 50.5     6.06     6.78
#3400m       6.13 0.117 72.0     5.89     6.36
#4200m       6.60 0.285 24.5     6.02     7.19
#1             estimate    SE  df t.ratio p.value
#2600m - 3400m   0.295 0.226 92.2  1.303  0.3972  
#2600m - 4200m  -0.181 0.344 46.8 -0.525  0.8593  
#3400m - 4200m  -0.476 0.305 31.5 -1.563  0.2763  

mean(dat_BM[dat_BM$Elevation=="2600m",]$BM) #6.125g
mean(dat_BM[dat_BM$Elevation=="3400m",]$BM) #5.112g
mean(dat_BM[dat_BM$Elevation=="4200m",]$BM) #5.150g
mean(dat_BM$BM)#6.31g

#####SVL
leveneTest(dat=dat_BM,SVL~Elevation) #p=0.075
mod_SVL<-aov(SVL~Elevation,data=dat_BM)
summary(mod_SVL)
anova(mod_SVL)
emmeans(mod_SVL, list(pairwise~Elevation), adjust = "tukey")



##Selected body temperature
#Read in data
Tsel<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/Tsel.csv")[,c(1,2,5,8)]
colnames(Tsel)<-c("Elevation","SVL","Gravid","Tsel")

mod_Tsel<-lm(data=Tsel,Tsel~SVL)
summary(mod_Tsel)#R²=0.018,p=0.305


dat_Tsel<-Tsel[Tsel$Gravid=="NO",c(1,4)]
#Run Levene's test (check the homogeneity of variances)
#library(car) #for the leveneTest
leveneTest(dat=Tsel,Tsel~Elevation) #p=0.073
mod_Tsel<-aov(data=dat_Tsel,Tsel~Elevation)
summary(mod_Tsel)
shapiro.test(resid(mod_Tsel))#p-value = 0.217
emmeans(mod_Tsel, list(pairwise~Elevation), adjust = "tukey")

# 1             estimate    SE df t.ratio p.value
#2600m - 3400m   -0.836 0.421 32 -1.984  0.1324 
#2600m - 4200m    0.975 0.317 50  3.078  0.0093 
#3400m - 4200m    1.621 0.312 50  5.192  <.0001 



mean(dat_Tsel[dat_Tsel$Elevation=="4200m",]$Tsel) #34.66
mean(dat_Tsel[dat_Tsel$Elevation!="4200m",]$Tsel) #35.97  alter in the second round




########Body temperatures#####
#Read in data
dat_Tb<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/Tbody.csv")
dat_NTb<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/Burrow.csv")

mod_TbSVL<-lm(data=dat_Tb,Tbody~Elevation+SVL)
summary(mod_TbSVL)#p=0.226

mod_NTb<-lm(data=dat_NTb,Tb~Elevation+SVL)
summary(mod_NTb)#p=0.503


#####The SVL females at sexual maturity 
Gravid_Tb<-dat_Tb[dat_Tb$Gravid=="YES",c(3,6)]
Gravid_NTb<-dat_NTb[dat_NTb$Gravid=="YES",c(5,7)]
Gravid<-na.omit(rbind(Gravid_Tb,Gravid_NTb))
Gravid_2600m<-Gravid[Gravid$Elevation=="2600m",]
Gravid_3400m<-Gravid[Gravid$Elevation=="3400m",]
Gravid_4200m<-Gravid[Gravid$Elevation=="4200m",]
range(Gravid_2600m$SVL)#4.9-7.2
range(Gravid_3400m$SVL)#4.1-6,5
range(Gravid_4200m$SVL)#4.7-6.5

#Remove data from gravid female (with lower their body temperature)
dat_Tb<-dat_Tb[dat_Tb$Gravid=="NO",c(3,4,5,6,7)]#dat_Tb$SVL>5 & 
dat_NTb<-dat_NTb[dat_NTb$Gravid=="NO",c(5,6,7,8)]#dat_NTb$SVL>5

table(dat_NTb$Elevation)
table(dat_Tb$Elevation)

#Keep the central 90% of data for each population
#####according to the Tsel, we combined the Tbody data of 2600m and 3400m.
i<-"4200m"
dat_Tb_90<-as.data.frame(matrix(NA,ncol=3,nrow=0))

  Tb_High<-dat_Tb[dat_Tb$Elevation==i,]
  Tb_High<-Tb_High[Tb_High$Tbody>=quantile(Tb_High$Tbody,0.05)&Tb_High$Tbody<=quantile(Tb_High$Tbody,0.95),]
  dat_Tb_90<-rbind(dat_Tb_90,Tb_High)
  
  Tb_Low<-dat_Tb[dat_Tb$Elevation!=i,]
  Tb_Low<-Tb_Low[Tb_Low$Tbody>=quantile(Tb_Low$Tbody,0.05)&Tb_Low$Tbody<=quantile(Tb_Low$Tbody,0.95),]
  dat_Tb_90<-rbind(dat_Tb_90,Tb_Low)
  
dat_NTb_90<-as.data.frame(matrix(NA,ncol=3,nrow=0))

NTb_High<-dat_NTb[dat_NTb$Elevation==i,]
NTb_High<-NTb_High[NTb_High$Tb>=quantile(NTb_High$Tb,0.05)&NTb_High$Tb<=quantile(NTb_High$Tb,0.95),]
dat_NTb_90<-rbind(dat_NTb_90,NTb_High)

NTb_Low<-dat_NTb[dat_NTb$Elevation!=i,]
NTb_Low<-NTb_Low[NTb_Low$Tb>=quantile(NTb_Low$Tb,0.05)&NTb_Low$Tb<=quantile(NTb_Low$Tb,0.95),]
dat_NTb_90<-rbind(dat_NTb_90,NTb_Low)


#Keep the central 99% of data for each population
dat_Tb_99<-as.data.frame(matrix(NA,ncol=3,nrow=0))

  Tb_High<-dat_Tb[dat_Tb$Elevation==i,]
  Tb_High<-Tb_High[Tb_High$Tbody>=quantile(Tb_High$Tbody,0.005)&Tb_High$Tbody<=quantile(Tb_High$Tbody,0.995),]
  dat_Tb_99<-rbind(dat_Tb_99,Tb_High)
  
  Tb_Low<-dat_Tb[dat_Tb$Elevation!=i,]
  Tb_Low<-Tb_Low[Tb_Low$Tbody>=quantile(Tb_Low$Tbody,0.005)&Tb_Low$Tbody<=quantile(Tb_Low$Tbody,0.995),]
  dat_Tb_99<-rbind(dat_Tb_99,Tb_Low)
  
  
dat_NTb_99<-as.data.frame(matrix(NA,ncol=3,nrow=0))
NTb_High<-dat_NTb[dat_NTb$Elevation==i,]
NTb_High<-NTb_High[NTb_High$Tb>=quantile(NTb_High$Tb,0.005)&NTb_High$Tb<=quantile(NTb_High$Tb,0.995),]
dat_NTb_99<-rbind(dat_NTb_99,NTb_High)

NTb_Low<-dat_NTb[dat_NTb$Elevation!=i,]
NTb_Low<-NTb_Low[NTb_Low$Tb>=quantile(NTb_Low$Tb,0.005)&NTb_Low$Tb<=quantile(NTb_Low$Tb,0.995),]
dat_NTb_99<-rbind(dat_NTb_99,NTb_Low)

range(dat_Tb_90$Tbody[dat_Tb_90$Elevation=="4200m"])
range(dat_Tb_99$Tbody[dat_Tb_99$Elevation=="4200m"])
range(dat_NTb_90$Tb[dat_NTb_90$Elevation=="4200m"])
range(dat_NTb_99$Tb[dat_NTb_99$Elevation=="4200m"])

range(dat_Tb_90$Tbody[dat_Tb_90$Elevation!="4200m"])
range(dat_Tb_99$Tbody[dat_Tb_99$Elevation!="4200m"])
range(dat_NTb_90$Tb[dat_NTb_90$Elevation!="4200m"])
range(dat_NTb_99$Tb[dat_NTb_99$Elevation!="4200m"])


#Cold tolerance
#Read in data
CTmin<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/CTm.csv")[,c(1,3,5)]

mod_SVLCTmin<-lm(data=CTmin,CTmin~SVL)
summary(mod_SVLCTmin)#R²=0.005,p=0.598

colnames(CTmin)<-c("Elevation","SVL","CTmin")
#dat_CTmin<-CTmin[CTmin$SVL>4.9 ,c(1:3)]
##、Run Levene's test (check the homogeneity of variances)
leveneTest(dat=CTmin,CTmin~Elevation) #p=0.048

aov_CTmin<-aov(data=CTmin,CTmin~Elevation)
summary(aov_CTmin)
#Run Shapiro-Wilk test (check the normality of residuals)
shapiro.test(x=residuals(object=aov_CTmin)) #p=0.365
emmeans(aov_CTmin, list(pairwise~Elevation), adjust = "tukey")
#1             estimate    SE df t.ratio p.value
#2600m - 3400m   0.0625 0.357 61 0.175   0.9832 
#2600m - 4200m   0.8875 0.357 61 2.488   0.0408  
#3400m - 4200m   0.8250 0.373 61 2.214   0.0768 
mean(CTmin[CTmin$Elevation=="4200m",]$CTmin) #3.9
mean(CTmin[CTmin$Elevation!="2600m",]$CTmin) #4.3





#Heatt tolerance
#Read in data
CTmax<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/CTm.csv")[,c(1,3,6)]
colnames(CTmax)<-c("Elevation","SVL","CTmax")

mod_SVLCTmax<-lm(data=CTmax,CTmax~SVL)
summary(mod_SVLCTmax)#R²=0.11,P=0.008

dat_CTmax<-CTmax[CTmax$SVL>5 ,c(1:3)]
##、Run Levene's test (check the homogeneity of variances)
leveneTest(dat=dat_CTmax,CTmax~Elevation) #p=0.83
#Compute one-way ANOVA test
aov_CTmax<-aov(data=dat_CTmax,CTmax~SVL+Elevation)#F2,27=0.768,p=0.4738
summary(aov_CTmax)

#Run Shapiro-Wilk test (check the normality of residuals)
shapiro.test(x=residuals(object=aov_CTmax)) #p=0.273
emmeans(aov_CTmax, list(pairwise~Elevation), adjust = "tukey")
# 1             estimate    SE df t.ratio p.value
#2600m - 3400m    0.154 0.291 28  0.527  0.8588   
#2600m - 4200m   -0.166 0.281 28 -0.591  0.8263
#3400m - 4200m   -0.319 0.319 28 -1.000  0.5830 

#library(effects)
#####remove the effects of covariates
#mean_CTmax<-mean(effect("Elevation",aov_CTmax)[["fit"]])
mean(dat_CTmax$CTmax)#p=45.89


#####metabolic rate
#####Model built according to Andrews, R. M., & Pough, F. H. (1985). Metabolism of squamate reptiles: allometric and ecological relationships. Physiological Zoology, 58(2), 214-231.
####MR from field data collection

MR<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/MRtest.csv")#### the best and final data is MRtest 
#MR_ak<-read.csv("/Users/jiangzhongwen/Desktop/2021/QT/Zhong-wen/Data/Aksai/MR_ak.csv")

MR$No.<-as.factor(MR$No.)
MR$Tb<-as.factor(MR$Tb)
#MR$BM<-as.factor(MR$BM)
#MR<-subset(MR,No.!=1&No.!=6&No.!=7&No.!=10)
#MR<-subset(MR,No.!=7)
MR2600m<-subset(MR,Elevation=="2600m")
MR3400m<-subset(MR,Elevation=="3400m")
MR4200m<-subset(MR,Elevation=="4200m")

#MR4200m<-subset(MR4200m,BM!="8.76"&BM!="3.64"&BM!="7.22")
leveneTest(log10(MR$V)~MR$Elevation)#p=0.492
#MR_mod<-aov(log10(V) ~log10(BM)+Tb*Elevation+Error(No.),data = MR)
leveneTest(MR$V..BM~MR$Elevation)#0.681
#MR_mod<-lmer(log10(V..BM)~Elevation+(1|No.),data = MR)
MR_mod<-lme(log10(V..BM)~BM+Elevation,random = ~1|No./Tb,data = MR)


summary(MR_mod)
anova(MR_mod)

shapiro.test(residuals(object=MR_mod))#0.01017
qqnorm(residuals(object=MR_mod))
hist(residuals(object=MR_mod))

library(emmeans)
emmeans(MR_mod, list(pairwise~Elevation), adjust = "tukey")
#$`pairwise differences of Elevation`
#1             estimate     SE df t.ratio p.value
#2600m - 3400m  0.0596 0.0766 26  0.778  0.7195 
#2600m - 4200m  -0.2052 0.0744 26 -2.760  0.0273 
#3400m - 4200m  -0.2648 0.0742 26 -3.567  0.0040 



library(plyr)
library(dplyr)

MR_sum<-ddply(MR,c("Elevation"),summarise,N=sum(!is.na(V..BM)),
                 mean=mean(V..BM,na.rm=T),
                 sd=sd(V..BM,na.rm=T),
                 se=sd/sqrt(N))




#2600m
MR2600m_mod<-lm(log10(V)~log10(BM)+Tb,data = MR2600m)
summary(MR2600m_mod)
#R²=0.7521, lgV=0.059479lgBM+0.040919Tb-1.372556; V=0.043*BM^0.059*10^0.041*Tb

#3400m
MR3400m_mod<-lm(log10(V)~log10(BM)+Tb,data = MR3400m)
summary(MR3400m_mod)
#R²=0.8017, lgV=0.463171lgBM+0.038788Tb-1.670339; V=0.021*BM^0.46*10^0.039*Tb

#####2600m&3400m
MR_low<-subset(MR,Elevation!="4200m")
MR_low_mod<-lm(log10(V)~log10(BM)+Tb,data = MR_low)
summary(MR_low_mod)
#R²=0.77, lgV=0.32lgBM+0.040Tb-1.574;V=0.027*BM^0.32*10^0.04*Tb


#4200m
MR4200m_mod<-lm(log10(V)~log10(BM)+Tb,data = MR4200m)
summary(MR4200m_mod)
#R²=0.808, lgV=1.020382lgBM+0.031665Tb-1.599509; V=0.025*BM^1.02*10^0.032*Tb







########Figures of behaviourial and physiological data########
library(ggplot2)
library(plyr)


#####Figure 2a Active body temperatures 
#dat_Tb<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/Tbody.csv")
dat_Tb_Low<-dat_Tb_99[dat_Tb_99$Elevation!="4200m",]
dat_Tb_High<-dat_Tb_99[dat_Tb_99$Elevation=="4200m",]
dat_Tb_Low["pop"]<-"Low"
dat_Tb_High["pop"]<-"High"
pop_Tb<-rbind(dat_Tb_High,dat_Tb_Low)
#pop_Tb<-pop_Tb[pop_Tb$Tbody>21,]
table(dat_Tb_99$Elevation)

Tb_sum<-ddply(pop_Tb,c("pop","Time"),summarise,N=sum(!is.na(Tbody)),
                 mean=mean(Tbody,na.rm=T),
                 sd=sd(Tbody,na.rm=T),
                 se=sd/sqrt(N))

Figure2a<-ggplot() +
  geom_line(data=Tb_sum,aes(x = Time, y =mean,group = pop,color=pop),size=1,position=position_dodge(0.3))+
  geom_point(data=Tb_sum,aes(x = Time, y =mean,group = pop,shape=pop,color=pop),size=6,position=position_dodge(0.3))+
  geom_errorbar(data=Tb_sum,size=1,aes(x = Time,ymin=mean-se, ymax=mean+se,color=pop), width=1,
                position=position_dodge(0.3))+
  #scale_color_manual(values=c("black","black"))+
  #scale_shape_manual(values = c(16,2))+
  scale_color_manual(values=c("#4169E1","#DC143C"))+
  scale_x_continuous(breaks = seq(8,20,1))+
  #scale_y_continuous(limits=c(20,40),breaks = seq(20.0,40.0,5.0))+
  labs(x="Hour", y = expression(bold(paste("Body temperature ("^"o","C)"))))+
  theme_classic()+
  theme(legend.title=element_blank(),
        #legend.position=c(0.1,0.85),
        legend.position="right",
        legend.text=element_text(size = rel(1),face=c("bold")),
        axis.text = element_text(size = rel(1),color="black",face=c("bold")),
        axis.title = element_text(size = rel(1),face=c("bold")),
        axis.line = element_line(colour = 'black', size = 0.8),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(size = 0.8))+
  theme(strip.text = element_text(colour = "black",family = "serif",face = "bold", size = 14), strip.background = element_rect(linetype = 0 ))+
  theme(plot.subtitle = element_text(family = "serif", 
                                     size = 19, face = "bold", colour = "black", 
                                     hjust = 0.5), plot.caption = element_text(family = "serif", 
                                                                               face = "bold", hjust = 0.5), axis.title = element_text(family = "serif"), 
        axis.text = element_text(family = "serif"), 
        axis.text.x = element_text(family = "serif"), 
        axis.text.y = element_text(family = "serif"), 
        plot.title = element_text(family = "serif"), 
        legend.text = element_text(family = "serif")) + 
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 24), 
        plot.title = element_text(size = 24))+ 
  theme(legend.text = element_text(size = 22), 
        legend.background = element_rect(fill = NA))+ 
  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0.025, 
                                  vjust = -5))+ 
  #theme(legend.position ="none")+
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 24), 
        legend.text = element_text(size = 24)) +
  theme(plot.title = element_text(size = 24))

  ggsave(Figure2a,file="/Volumes/My Passport/QT/Result/Figure2a.pdf",device = "pdf",height = 6,width = 9,dpi = 300)

  
  #####Figure 2b CTmin#######
  CTmin<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/CTm.csv")[,c(1,3,5)]
  colnames(CTmin)<-c("Elevation","SVL","CTmin")
  Low_CTmin<-CTmin[CTmin$Elevation!="4200m",c(2,3)]
  High_CTmin<-CTmin[CTmin$Elevation=="4200m",c(2,3)]
  Low_CTmin["pop"]<-"Low"
  High_CTmin["pop"]<-"High"
  pop_CTmin<-rbind(Low_CTmin,High_CTmin)


  Figure2b<-ggplot(pop_CTmin, aes(x=pop, y=CTmin, color=pop)) +
    labs(x="Elevation")+
    geom_boxplot(size=1)+
    scale_color_manual(values=c("#4169E1","#DC143C"))+
    labs(x="Elevation", y = "CTmin")+
    theme_classic()+
    theme(legend.title=element_blank(),
          #legend.position=c(0.1,0.85),
          legend.position="right",
          legend.text=element_text(size = rel(1),face=c("bold")),
          axis.text = element_text(size = rel(1),color="black",face=c("bold")),
          axis.title = element_text(size = rel(1),face=c("bold")),
          axis.line = element_line(colour = 'black', size = 0.8),
          axis.ticks.length=unit(.12, "cm"),
          axis.ticks=element_line(size = 0.8))+
    theme(strip.text = element_text(colour = "black",family = "serif",face = "bold", size = 14), strip.background = element_rect(linetype = 0 ))+
    theme(plot.subtitle = element_text(family = "serif", 
                                       size = 19, face = "bold", colour = "black", 
                                       hjust = 0.5), plot.caption = element_text(family = "serif", 
                                                                                 face = "bold", hjust = 0.5), axis.title = element_text(family = "serif"), 
          axis.text = element_text(family = "serif"), 
          axis.text.x = element_text(family = "serif"), 
          axis.text.y = element_text(family = "serif"), 
          plot.title = element_text(family = "serif"), 
          legend.text = element_text(family = "serif")) + 
    theme(axis.title = element_text(size = 24), 
          axis.text = element_text(size = 24), 
          plot.title = element_text(size = 24))+ 
    theme(legend.text = element_text(size = 22), 
          legend.background = element_rect(fill = NA))+ 
    theme(plot.title = element_text(face = "bold")) +
    theme(plot.title = element_text(hjust = 0.025, 
                                    vjust = -5))+ 
    #theme(legend.position ="none")+
    theme(axis.title = element_text(size = 24), 
          axis.text = element_text(size = 24), 
          legend.text = element_text(size = 24)) +
    theme(plot.title = element_text(size = 24))
  
  ggsave(Figure2b,file="/Volumes/My Passport/QT/Result/Figure2b.pdf",device = "pdf",height = 6,width = 9,dpi = 300)
  
  
  #####Figure 2c metabolic rate #######
  MR<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/MRtest.csv")[,c(1,2,4,7,14)]#### the best and final data is MRtest 
  colnames(MR)<-c("ID","Ele","Tb","BM","MR")
  Low_MR<-MR[MR$Ele!="4200m",]
  High_MR<-MR[MR$Ele=="4200m",]
  Low_MR["pop"]<-"Low"
  High_MR["pop"]<-"High"
  pop_MR<-rbind(Low_MR,High_MR)
  
  popMR_sum<-ddply(pop_MR,c("Tb","pop"),summarise,N=sum(!is.na(MR)),
                mean=mean(MR,na.rm=T),
                sd=sd(MR,na.rm=T),
                se=sd/sqrt(N))
  
  Figure2c<-ggplot() +
    geom_line(data=popMR_sum,aes(x = Tb, y =mean,group = pop,color=pop),size=1,position=position_dodge(0.3))+
    geom_point(data=popMR_sum,aes(x = Tb, y =mean,group = pop,shape=pop,color=pop),size=6,position=position_dodge(0.3))+
    geom_errorbar(data=popMR_sum,size=1,aes(x = Tb,ymin=mean-se, ymax=mean+se,color=pop), width=2.5,
                  position=position_dodge(0.3))+
    #scale_color_manual(values=c("black","black"))+
    #scale_shape_manual(values = c(16,2))+
    scale_color_manual(values=c("#4169E1","#DC143C"))+
    scale_x_continuous(breaks = seq(10,40,5))+
    #scale_y_continuous(limits=c(20,40),breaks = seq(20.0,40.0,5.0))+
    labs(x="Body temperature", y = "MR(ml/g/h)")+
    theme_classic()+
    theme(legend.title=element_blank(),
          #legend.position=c(0.1,0.85),
          legend.position="right",
          legend.text=element_text(size = rel(1),face=c("bold")),
          axis.text = element_text(size = rel(1),color="black",face=c("bold")),
          axis.title = element_text(size = rel(1),face=c("bold")),
          axis.line = element_line(colour = 'black', size = 0.8),
          axis.ticks.length=unit(.12, "cm"),
          axis.ticks=element_line(size = 0.8))+
    theme(strip.text = element_text(colour = "black",family = "serif",face = "bold", size = 14), strip.background = element_rect(linetype = 0 ))+
    theme(plot.subtitle = element_text(family = "serif", 
                                       size = 19, face = "bold", colour = "black", 
                                       hjust = 0.5), plot.caption = element_text(family = "serif", 
                                                                                 face = "bold", hjust = 0.5), axis.title = element_text(family = "serif"), 
          axis.text = element_text(family = "serif"), 
          axis.text.x = element_text(family = "serif"), 
          axis.text.y = element_text(family = "serif"), 
          plot.title = element_text(family = "serif"), 
          legend.text = element_text(family = "serif")) + 
    theme(axis.title = element_text(size = 24), 
          axis.text = element_text(size = 24), 
          plot.title = element_text(size = 24))+ 
    theme(legend.text = element_text(size = 22), 
          legend.background = element_rect(fill = NA))+ 
    theme(plot.title = element_text(face = "bold")) +
    theme(plot.title = element_text(hjust = 0.025, 
                                    vjust = -5))+ 
    #theme(legend.position ="none")+
    theme(axis.title = element_text(size = 24), 
          axis.text = element_text(size = 24), 
          legend.text = element_text(size = 24)) +
    theme(plot.title = element_text(size = 24))
  
  ggsave(Figure2c,file="/Volumes/My Passport/QT/Result/Figure2c.pdf",device = "pdf",height = 6,width = 9,dpi = 300)
  
  #####Figure 2d Active body temperatures and preferred body temperatures
  Tsel<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/Tsel.csv")[,c(1,2,5,8)]
  colnames(Tsel)<-c("Elevation","SVL","Gravid","Tsel")
  dat_Tsel<-Tsel[Tsel$Gravid=="NO",c(1,4)]
  Low_Tsel<-dat_Tsel[dat_Tsel$Elevation!="4200m",]
  High_Tsel<-dat_Tsel[dat_Tsel$Elevation=="4200m",]
  Low_Tsel["pop"]<-"Low"
  High_Tsel["pop"]<-"High"
  pop_Tsel<-rbind(Low_Tsel,High_Tsel)
  
  
  Figure2d<-ggplot(pop_Tsel, aes(x=pop, y=Tsel, color=pop)) +
    labs(x="Elevation")+
    geom_boxplot(size=1)+
    scale_color_manual(values=c("#4169E1","#DC143C"))+
    labs(x="Elevation", y = "Tsel")+
    theme_classic()+
    theme(legend.title=element_blank(),
          #legend.position=c(0.1,0.85),
          legend.position="right",
          legend.text=element_text(size = rel(1),face=c("bold")),
          axis.text = element_text(size = rel(1),color="black",face=c("bold")),
          axis.title = element_text(size = rel(1),face=c("bold")),
          axis.line = element_line(colour = 'black', size = 0.8),
          axis.ticks.length=unit(.12, "cm"),
          axis.ticks=element_line(size = 0.8))+
    theme(strip.text = element_text(colour = "black",family = "serif",face = "bold", size = 14), strip.background = element_rect(linetype = 0 ))+
    theme(plot.subtitle = element_text(family = "serif", 
                                       size = 19, face = "bold", colour = "black", 
                                       hjust = 0.5), plot.caption = element_text(family = "serif", 
                                                                                 face = "bold", hjust = 0.5), axis.title = element_text(family = "serif"), 
          axis.text = element_text(family = "serif"), 
          axis.text.x = element_text(family = "serif"), 
          axis.text.y = element_text(family = "serif"), 
          plot.title = element_text(family = "serif"), 
          legend.text = element_text(family = "serif")) + 
    theme(axis.title = element_text(size = 24), 
          axis.text = element_text(size = 24), 
          plot.title = element_text(size = 24))+ 
    theme(legend.text = element_text(size = 22), 
          legend.background = element_rect(fill = NA))+ 
    theme(plot.title = element_text(face = "bold")) +
    theme(plot.title = element_text(hjust = 0.025, 
                                    vjust = -5))+ 
    #theme(legend.position ="none")+
    theme(axis.title = element_text(size = 24), 
          axis.text = element_text(size = 24), 
          legend.text = element_text(size = 24)) +
    theme(plot.title = element_text(size = 24))
  
  ggsave(Figure2d,file="/Volumes/My Passport/QT/Result/Figure2d.pdf",device = "pdf",height = 6,width = 9,dpi = 300)
  
  
  
  
  
  
  
  
  
  
  
  #########translating the dorsal skin temperature to cloocal temperatures
  library("ggpmisc")
  library(ggplot2)
  
  #Correct the infrared gun in 2600m
  Correct_2600m<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/Correct_2600m.csv")
  
  formula <- y ~ x     #####For smooth   require(ggpmisc)
  gg_gun2600m = ggplot(Correct_2600m, aes(Ts,Tc)) +
    geom_point(alpha=.8,size=3)+
    ggtitle("2600m")+
    geom_smooth(method=lm,formula = formula,col="black")+  #####添加拟合， se=FALSE 可以去掉置信区间线,aes(fill=GainLoss)，用于分组
    xlab("Skin temperatures")+
    ylab("Cloacal temperatures")+
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), ####添加公式和R2
                 label.x.npc = "left",# label.y.npc = 0.15,
                 formula = formula, parse = TRUE, size = 3.5,rr.digits = 2)+   ####不同颜色的公式
    ####出总公式的话添加，color="black"
    # geom_abline(intercept = 0,slope = 1,lwd=2,lty=2,col="blue")+
    theme_light()+
    ggtitle("2600m")+
    theme(plot.title = element_text(hjust = 0.5))     
  
  
  #Correct the infrared gun in 3400m
  Correct_3400m<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/Correct_3400m.csv")
  mod_correct_3400m<-lm(Tc~Ts,data = Correct_3400m)
  summary(mod_correct_3400m)#c1 y=0.84863x+5.25413 R²>0.9    #c2 y=0.76906x+7.91610 R²=0.8428
  formula <- y ~ x     #####For smooth   require(ggpmisc)
  gg_gun3400m = ggplot(Correct_3400m, aes(Ts,Tc)) +
    geom_point(alpha=.8,size=3)+
    geom_smooth(method=lm,formula = formula,col="black")+  #####添加拟合， se=FALSE 可以去掉置信区间线,aes(fill=GainLoss)，用于分组
    xlab("Skin temperatures")+
    ylab("Cloacal temperatures")+
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), ####添加公式和R2
                 label.x.npc = "left",# label.y.npc = 0.15,
                 formula = formula, parse = TRUE, size = 3.5,rr.digits = 2)+   ####不同颜色的公式
    ####出总公式的话添加，color="black"
    #geom_abline(intercept = 0,slope = 1,lwd=2,lty=2,col="blue")+
    theme_light()+
    ggtitle("3400m")+
    theme(plot.title = element_text(hjust = 0.5))   
  
  #Correct the infrared gun in 4200m
  Correct_4200m<-read.csv("/Volumes/My Passport/QT/NicheMapR/Data/Correct_4200m.csv")
  mod_correct_4200m<-lm(Tc~Ts,data = Correct_4200m)
  summary(mod_correct_4200m)
  
  formula <- y ~ x     #####For smooth   require(ggpmisc)
  gg_gun4200m = ggplot(Correct_4200m, aes(Ts,Tc)) +
    geom_point(alpha=.8,size=3)+
    geom_smooth(method=lm,formula = formula,col="black")+  #####添加拟合， se=FALSE 可以去掉置信区间线,aes(fill=GainLoss)，用于分组
    xlab("Skin temperatures")+
    ylab("Cloacal temperatures")+
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), ####添加公式和R2
                 label.x.npc = "left",# label.y.npc = 0.15,
                 formula = formula, parse = TRUE, size = 3.5,rr.digits = 2)+   ####不同颜色的公式
    ####出总公式的话添加，color="black"
    #geom_abline(intercept = 0,slope = 1,lwd=2,lty=2,col="blue")+
    theme_light()+
    ggtitle("4200m")+
    theme(plot.title = element_text(hjust = 0.5))   
  
  library(ggpubr)
  pdf("/Volumes/My Passport/QT/Result/gg_gun.pdf",height = 4,width = 12)
  
  ggarrange(gg_gun2600m,gg_gun3400m,gg_gun4200m,ncol = 3,nrow = 1,labels = c("a","b","c"))
  dev.off()
  
  
