#####EWL validation for QT_MS####
####zhong-wen Jiang Apr2-2022####
rm(list = ls())
gc()
library(NicheMapR)


EWL_Pv<-read.csv("/Users/jiangzhongwen/Desktop/2021/QT/QT_MS/validation_data/EWL.csv")

data_pct_wet<-data.frame()

for (i in 1:nrow(EWL_Pv)) {
  #i=1
  ID_test<-EWL_Pv[i,]
  
  wet<-as.data.frame(get_p_wet(T_b = ID_test$Tb,mass.g = ID_test$BM,E_tot = ID_test$EWL/3600,RHin = 15) )
  ID_test["ID"]<-i
  ID_test["pct_wet"]<-wet$p_wet
  
  data_pct_wet<-rbind(data_pct_wet,ID_test)
  
}

mean(data_pct_wet$pct_wet)
#RHin=0,pct_wetmean<-0.009479852
#RHin=5,pct_wetmean<-0.009978829
#RHin=10,pct_wetmean<-0.01053338
#RHin=15,pct_wetmean<-0.01115332
#mean=0.01028635


micro<-micro_global(run.gads = 2,timeinterval = 12)
EWL_Niche<-data.frame()



for (i in 1:nrow(EWL_Pv)) {
  
  ID_test<-EWL_Pv[i,]
  
  Temp=ID_test$Tb
  
  EWL_Niche_hum<-data.frame()
  
 # for (hum in c(0,5,10,15)) {
    hum=15
    metoutdata<-as.data.frame(micro$metout[1:24,])
    metoutdata$TALOC=Temp#air temperature (°C) at local height (specified by 'Usrhyt' variable)
    metoutdata$TAREF=Temp#air temperature (°C) at reference height (specified by 'Refhyt', 1.2m default)
    metoutdata$RHLOC=hum#relative humidity (%) at local height (specified by 'Usrhyt' variable)
    metoutdata$RH=hum   #relative humidity (%) at reference height (specified by 'Refhyt', 1.2m default)
    metoutdata$VLOC=0.00329 #wind speed (m/s) at local height (specified by 'Usrhyt' variable)
    metoutdata$VREF=0.00329 #wind speed (m/s) at reference height (specified by 'Refhyt', 1.2m default)
    metoutdata$SNOWMELT=0#snowmelt (mm)
    metoutdata$POOLDEP=0#water pooling on surface (mm)
    metoutdata$PCTWET=0 #soil surface wetness (%)
    metoutdata$ZEN=0 #zenith angle of sun (degrees - 90 = below the horizon)
    metoutdata$SOLR=0 #solar radiation (W/m2) (unshaded, horizontal plane)
    metoutdata$TSKYC=0#Temp #sky radiant temperature (°C)
    metoutdata$DEW=0 # dew presence (0 or 1)
    metoutdata$FROST=0 #frost presence (0 or 1)
    metoutdata$SNOWFALL=0#snow predicted to have fallen (cm)
    metoutdata$SNOWDEP=0 #predicted snow depth (cm)
    metoutdata$SNOWDENS=0 #snow density (g/cm3)
    metout<-as.matrix(metoutdata)
    shadmet<-as.matrix(metoutdata)
    
    soildata =as.data.frame(micro$soil[1:24,]) 
    soildata$D0cm=Temp
    soildata$D2.5cm=Temp
    soildata$D5cm=Temp
    soildata$D10cm=Temp
    soildata$D15cm=Temp
    soildata$D20cm=Temp
    soildata$D30cm=Temp
    soildata$D50cm=Temp
    soildata$D100cm=Temp
    soildata$D200cm=Temp
    soil<-as.matrix(soildata)
    
    soilmoistdata =as.data.frame(micro$soilmoist[1:24,]) 
    soildata$D0cm=0
    soildata$D2.5cm=0
    soildata$D5cm=0
    soildata$D10cm=0
    soildata$D15cm=0
    soildata$D20cm=0
    soildata$D30cm=0
    soildata$D50cm=0
    soildata$D100cm=0
    soildata$D200cm=0
    soilmoist<-as.matrix(soilmoistdata)
    
    humiddata =as.data.frame(micro$humid[1:24,]) 
    humiddata$D0cm=hum*0.01
    humiddata$D2.5cm=hum*0.01
    humiddata$D5cm=hum*0.01
    humiddata$D10cm=hum*0.01
    humiddata$D15cm=hum*0.01
    humiddata$D20cm=hum*0.01
    humiddata$D30cm=hum*0.01
    humiddata$D50cm=hum*0.01
    humiddata$D100cm=hum*0.01
    humiddata$D200cm=hum*0.01
    humid<-as.matrix(humiddata)
    
    lizard_loc<-ectotherm(
      live=1,
      Ww_g	= ID_test$BM, #default=40; Wet weight of animal (g), note this model is 'steady state' so no lags in heating/cooling due to mass
      M_1	=0.025, #trait_list$M_1, #Metabolic rate parameter 1 V_O2=M_1*M^M_2*10^(M_3*Tb), in ml O2 / h, default parameters for lizards based on Eq. 2 from Andrews & Pough 1985. Physiol. Zool. 58:214-231
      M_2	=1.02, #trait_list$M_2, #Metabolic rate parameter 2
      M_3	=0.032, #trait_list$M_3, #Metabolic rate parameter 3
      
      minshades=0,
      maxshades=0,
      pct_wet=1,#%of surface area acting as a free-water exchanger, for computing cutaneous water loss
      #pct_eyes=0.03,#% of surface area taken up by open eyes, for computing ocular water loss (only when active)
      #pct_mouth=5,#% of surface area taken up by open mouth, for computing panting water loss
      #pantmax=1,# maximum multiplier on breathing rate, for respiratory water loss via panting
      #F_O2=20,#% oxygen extraction efficiency, for respiratory water loss
      #delta_air=0.1,#Temperature difference (??C) between expired and inspired air, for computing respiratory water loss
      
      #Details
      #Parameters controling how the model runs:
      
      metout = metout, #Microclimate model output for above ground, minimum shade conditions
      shadmet = shadmet, #Microclimate model output for above ground, maximum shade conditions
      soil = soil, #Microclimate model output for soil temperature, minimum shade conditions
      shadsoil = soil, #Microclimate model output for soil temperature, maximum shade conditions
      soilmoist = soilmoist, #Microclimate model output for soil moisture, minimum shade conditions
      shadmoist = soilmoist, #Microclimate model output for soil moisture, maximum shade conditions
      humid = humid, #Microclimate model output for soil humidity, minimum shade conditions
      shadhumid = humid, #Microclimate model output for soil humidity, maximum shade conditions
      #soilpot = 0, #Microclimate model output for soil water potential, minimum shade conditions
      #shadpot = 0, #Microclimate model output for soil water potential, maximum shade conditions
      rainfall = 0, #Vector of daily rainfall (mm)
      rainhr = 0, #Vector of hourly rainfall (mm), overwrites rainfall if not negative
    )
    
    lizard_environ_out<-data.frame(lizard_loc[["environ"]])
    lizard_masbal_out<-data.frame(lizard_loc[["masbal"]])[1,15:17]
    lizard_masbal_out["TA"]<-Temp
    #lizard_masbal_out["Humidy"]<-hum
    lizard_masbal_out["TC"]<-lizard_environ_out$TC[1]
    lizard_masbal_out["ID"]<-i
    lizard_masbal_out["EWL_model"]<-lizard_masbal_out$H2OResp_g+lizard_masbal_out$H2OCut_g+lizard_masbal_out$H2OEye_g
    
    EWL_Niche_hum<-rbind(EWL_Niche_hum,lizard_masbal_out)
    
  
  EWL_Niche_hum<-cbind(EWL_Niche_hum[1,c(4,7)],TB=mean(EWL_Niche_hum$TC),EWL_model=mean(EWL_Niche_hum$EWL_model))
  EWL_Niche<-rbind(EWL_Niche,EWL_Niche_hum)
  
  
  
}






EWL_Niche["EWL_PV"]<-EWL_Pv$EWL

write.csv(EWL_Niche,file = "/Users/jiangzhongwen/Desktop/2021/QT/QT_MS/EWL_validation.csv",row.names = FALSE)

mod<-lm(data = EWL_Niche,EWL_model~EWL_PV)
summary(mod)



library(ModelMetrics)
library(ggplot2)

EWL_Niche<-read.csv("/Users/jiangzhongwen/Desktop/2021/QT/QT_MS/validation_data/EWL_validation.csv")
cor_pb<-cor.test(EWL_Niche$EWL_model,EWL_Niche$EWL_PV)
rmse(EWL_Niche$EWL_model,EWL_Niche$EWL_PV)

library(plyr)

EWLmodel_sum<-ddply(EWL_Niche,c("TA"),summarise,N=sum(!is.na(EWL_model)),
                    mean_mod=mean(EWL_model,na.rm=T),
                    sd_mod=sd(EWL_model,na.rm=T),
                    se_mod=sd_mod/sqrt(N))

EWLPV_sum<-ddply(EWL_Niche,c("TA"),summarise,N=sum(!is.na(EWL_PV)),
                    mean_obs=mean(EWL_PV,na.rm=T),
                    sd_obs=sd(EWL_PV,na.rm=T),
                    se_obs=sd_obs/sqrt(N))

cor_pb<-cor.test(EWLmodel_sum$mean_mod,EWLPV_sum$mean_obs)


Fig_validation_EWL<-ggplot()+
  #geom_line(data=EWLPV_sum,aes(x = TA, y =mean_obs),size=1,position=position_dodge(0.3),color="blue")+
  geom_point(data=EWLPV_sum,aes(x = TA, y =mean_obs),size=6,position=position_dodge(0.3),color="blue")+
  geom_errorbar(data=EWLPV_sum,size=1,aes(x = TA,ymin=mean_obs-sd_obs, ymax=mean_obs+sd_obs),color="blue", width=1,
                position=position_dodge(0.3))+
geom_line(data=EWLmodel_sum,aes(x = TA, y =mean_mod),size=1,position=position_dodge(0.3),color="red")+
  geom_point(data=EWLmodel_sum,aes(x = TA, y =mean_mod),size=6,position=position_dodge(0.3),color="red")+
  geom_errorbar(data=EWLmodel_sum,size=1,aes(x = TA,ymin=mean_mod-sd_mod, ymax=mean_mod+sd_mod),color="red", width=1,
               position=position_dodge(0.3))+

  theme_classic()+
  #ylim(0,0.04)+
 # scale_x_continuous(limits=c(25,40),breaks = seq(20.0,45.0,5.0))+
  ylab("Evaporated water loss (g/h)")+
  xlab("Air temperatures (°C)")+
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
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 18), 
        plot.title = element_text(size = 18))+ 
  theme(legend.text = element_text(size = 18), 
        legend.background = element_rect(fill = NA))+ 
  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0.025, 
                                  vjust = -5))+ 
  #theme(legend.position ="none")+
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 18)) +
  theme(plot.title = element_text(size = 18))




ggsave(Fig_validation_EWL,file="/Users/jiangzhongwen/Desktop/2021/QT/QT_MS/Figure2.pdf",device = "pdf",height = 6,width = 6,dpi = 300)

#°C

