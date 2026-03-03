#####Response plot#####

library(ggplot2)
library(usdm)
library(ecospat)
library(dplyr)
library(raster)
library(reshape2)
library(rgeos)
library(rgdal)

limit="3600"

for (limit in c("3600","3700")) {
  


path_rds<-paste0("E:/QT/SDM/",limit,"m/SDM_Group/")
species="Phrynocephalus vlangalii"
varnames<-c("Bio3","Bio4","Bio14","Bio15","OC","WL")

#pop="High"

#####Figure S1 Response curves####
for (pop in c("Low","High")) {
  
  
  mod<-readRDS(paste(path_rds,"/",pop,"/",species,"_",pop,"_mod.rds",sep=""))
  mod.EM<-readRDS(paste0(path_rds,pop,"/",species,"_",pop,"_mod.EM.rds"))
  names(mod.EM$weights)<-sub("BIOMOD","PA1",names(mod.EM$weights))
  
  
  
  pdf(paste0(path_rds,pop,"_Fig_S1.pdf"),height=6,width=12)
  dat_rc<-ecospat.ESM.responsePlot(mod.EM,mod,fixed.var.metric = "median")
  dev.off()
  write.csv(dat_rc,paste0(path_rds,"dat_rc_",pop,".csv"),row.names = F)
  
  
}

#####Figure S2 (variable importance)####
dat_varimp_mean<-matrix(NA,ncol=3,nrow=0)
colnames(dat_varimp_mean)<-c("Population","Variable","Value")
for(pop in c("Low","High")){
  
  dat<-read.csv(paste0(path_rds,pop,"/Phrynocephalus vlangalii_",pop,"_EM_Contrib.csv"))
  varimp_i<-do.call(cbind.data.frame,list(rep(pop,6),varnames,dat$ENS))
  colnames(varimp_i)<-c("Population","Variable","Value")
  dat_varimp_mean<-rbind.data.frame(dat_varimp_mean,varimp_i)
}

dat_varimp_mean$Variable <- factor(dat_varimp_mean$Variable, levels=varnames)
write.csv(dat_varimp_mean,paste0(path_rds,"dat_varimp_mean.csv"),row.names=F)

pdf(paste0(path_rds,"Fig_S2.pdf"),height=6,width=12)
ggplot(dat_varimp_mean, aes(x = Variable, y = Value, fill = Population)) + 
  scale_fill_manual(values = c("#4169E1","#DC143C"))+
  geom_bar(stat = "identity")+
  ylab("Variable importance")+
  theme_light()+
  theme(strip.text = element_text(face = "italic",size = 11),
        axis.text.x = element_text(size = 11, angle = 45,margin = margin(t = 17, r = 0, b = 0, l = 0)),
        axis.title = element_text(size = 13),
        legend.position = "none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_grid(Population~ .)
dev.off()




dat_varimp_mean<-read.csv(paste0(path_rds,"dat_varimp_mean.csv"))
var_clim<-varnames[1:2]
mech_land<-varnames[3:4]
aov_list<-vector()
for(pop in c("Low","High")){
  
  dat_spe<-dat_varimp_mean[dat_varimp_mean$Population==pop,]
  dat_spe$category<-c(rep("clim",4),rep("mech",2))
  # Compute the analysis of variance
  res.aov <- aov(Value ~ category, data = dat_spe)
  # Summary of the analysis
  print(summary(res.aov))
}







#####Figure S3 The contribution of variables to the Lost and Colonized of suitable habitats####

####----Prepare data structure for figures----####

ssp126_change<-data.frame()
ssp585_change<-data.frame()


Vlangalii_IUCN<-readOGR("E:/QT/SDM/P. vlangalii/data_P.shp")
spe<-"Phrynocephalus.vlangalii"
IUCN_buffer<-gBuffer(Vlangalii_IUCN,width=0.5)
#Elevation
elmat = raster("E:/QT/SDM/wc2.1_10m_elev.tif")
elmat=crop(elmat,c(79,110,24,43))

#pop="Low"

all_changeresult<-data.frame()

for(pop in c("Low","High")){
  
  
  path="E:/QT/SDM/United/Niche_var/"
  env_1970_2000<-stack(paste(path,"Eco_",pop,"_convar_1970_2000_10m.tif",sep=""))
  varname<-read.csv("E:/QT/SDM/varname1.csv")$ ﻿V1
  names(env_1970_2000)<-varname
  env_1970_2000_sub<-env_1970_2000[[c(3,4,14,15,23,26)]]
  #env_1970_2000_sub.NoCor <- vifcor(env_1970_2000_sub, th = 0.7)
  env_1970_2000_sub.NoCor <- vifstep(env_1970_2000_sub, th = 5)
  env_1970_2000.Reduced <- exclude(env_1970_2000_sub, env_1970_2000_sub.NoCor)
  env_2081_2100_ssp126<-stack(paste(path,"Eco_",pop,"_ssp126_convar_2081_2100_10m.tif",sep=""))
  env_2081_2100_ssp585<-stack(paste(path,"Eco_",pop,"_ssp585_convar_2081_2100_10m.tif",sep=""))
  names(env_2081_2100_ssp126)<-varname
  names(env_2081_2100_ssp585)<-varname
  env_2081_2100_ssp126_Reduce<-exclude(env_2081_2100_ssp126,env_1970_2000_sub.NoCor)
  env_2081_2100_ssp585_Reduce<-exclude(env_2081_2100_ssp585,env_1970_2000_sub.NoCor)
  
  
  
  
  ####Lost/Colonized suitable habitat
  
  current_pre<-raster(paste0("E:/QT/SDM/",limit,"m/SDM_Group/",pop,"/TIFresult/","binary_current_",pop,".tif"))
  current_crop<-crop(current_pre,IUCN_buffer)
  current_mask<-mask(current_crop,IUCN_buffer)
  
  dat_current_mask<-as.data.frame(rasterToPoints(current_mask))
  
  
  dat_current_mask["Elevation"]<-extract(elmat,dat_current_mask[,c(1,2)])
  colnames(dat_current_mask)<-c("Longitude","Latitude","Preabs","Elevation") #Preabs: presence (1) vs absence (0);
  
  if(pop=="Low"){  dat_current_mask<-dat_current_mask[dat_current_mask$Elevation<limit,] 
  }else{  dat_current_mask<-dat_current_mask[dat_current_mask$Elevation>=limit,]}
  
  current_raster<-rasterFromXYZ(dat_current_mask[dat_current_mask$Preabs==1,]) 
  current_polygen<-rasterToPolygons(current_raster)
  future_buffer<-gBuffer(current_polygen,width=0.5)######the future dispersed boundary limited by predicted current with a 50km buffer n 
  
  
  future_ssp126_pre<-raster(paste0("E:/QT/SDM/",limit,"m/SDM_Group/",pop,"/TIFresult/","binary_ssp126_wc_",pop,".tif"))
  future_ssp585_pre<-raster(paste0("E:/QT/SDM/",limit,"m/SDM_Group/",pop,"/TIFresult/","binary_ssp585_wc_",pop,".tif"))
  
  
  current<-mask(current_pre,future_buffer)
  future_ssp126<-mask(future_ssp126_pre,future_buffer)
  future_ssp585<-mask(future_ssp585_pre,future_buffer)

  Current<-as.data.frame(rasterToPoints(current))
  future_ssp126<-as.data.frame(rasterToPoints(future_ssp126))
  future_ssp585<-as.data.frame(rasterToPoints(future_ssp585))
  
  
  Current<-na.omit(Current)
  future_ssp126<-na.omit(future_ssp126)
  future_ssp585<-na.omit(future_ssp585)
  
  Current["Elevation"]<-extract(elmat,Current[,c(1,2)])
  colnames(Current)<-c("Longitude","Latitude","EF","Elevation")
  colnames(future_ssp126)<-c("Longitude","Latitude","EF")
  colnames(future_ssp585)<-c("Longitude","Latitude","EF")
  
  if(pop=="Low"){Current[Current$EF==1&Current$Elevation>=limit,]$EF<-0}else{
    Current[Current$EF==1&Current$Elevation<limit,]$EF<-0}
  
  
  #Lost suitable habitat
  Current_SH<-Current[Current$EF==1,1:2]
  future_ssp126_Lost<-Current[Current$EF==1&future_ssp126$EF!=1,1:2]
  future_ssp585_Lost<-Current[Current$EF==1&future_ssp585$EF!=1,1:2]
  
  #Colonizeded suitable habitat
  future_ssp126_Colonized<-Current[Current$EF!=1&future_ssp126$EF==1,1:2]
  future_ssp585_Colonized<-Current[Current$EF!=1&future_ssp585$EF==1,1:2]
  
  
  dat_resp<-read.csv(paste0(path_rds,"dat_rc_",pop,".csv"))
  
  ####Occurrence probability induced by each variable at Lost/Colonizeded suitable habitat
  ##Occurrence probability induced by each ssp126 no adaptation assuming all other variables are constant(mean value) at Lost/Colonizeded suitable habitat (induced by only LUC)
  future_ssp126_Lost_past<-future_ssp126_Lost
  future_ssp126_Lost_future<-future_ssp126_Lost
  future_ssp126_Colonized_past<-future_ssp126_Colonized
  future_ssp126_Colonized_future<-future_ssp126_Colonized
  
  
  ssp126change_Lost<-vector()
  ssp126change_Colonized_Low<-vector()
  ssp126change_Colonized_High<-vector()
  #n=4
  
  ssp126_Lost_result<-data.frame()
  ssp126_Colonized_result<-data.frame()
  
  for(n in 1:6){
    dat_resp_n<-dat_resp[,c((n*5-4),n*5)]
    resp_fun<-approxfun(dat_resp_n[,1],dat_resp_n[,2]) #Linear interpolation
    
    ssp126_Lost_past<-raster::extract(env_1970_2000.Reduced[[n]],future_ssp126_Lost[,1:2], sp = T) #mean variable value at Lost habitat in the past
    ssp126_Lost_future<-raster::extract(env_2081_2100_ssp126_Reduce[[n]],future_ssp126_Lost[,1:2], sp = T) #mean variable value at Lost habitat in the future
    ssp126change_Lost[n]<-ifelse((mean(na.omit(ssp126_Lost_future))-mean(na.omit(ssp126_Lost_past)))>0,1,-1) #increase=1, decrease=-1
    
    future_ssp126_Lost["current"]<-raster::extract(env_1970_2000.Reduced[[n]],future_ssp126_Lost[,1:2], sp = T) #mean variable value at Lost habitat in the past
    future_ssp126_Lost["ssp126"]<-raster::extract(env_2081_2100_ssp126_Reduce[[n]],future_ssp126_Lost[,1:2], sp = T) #mean variable value at Lost habitat in the future
    future_ssp126_Lost["Traits"]<-varnames[n]
    ssp126_Lost<-melt(future_ssp126_Lost,id.vars = c("Longitude","Latitude","Traits"))
    ssp126_Lost_result<-rbind(ssp126_Lost_result,ssp126_Lost)
    
    ###replace the na by the last value
    test_Lost_past<-unlist(lapply(ssp126_Lost_past,resp_fun))
    test_Lost_past[is.na(test_Lost_past)]<- dat_resp_n[order(dat_resp_n[,1]),][nrow(dat_resp_n),2]
    test_Lost_future<-unlist(lapply(ssp126_Lost_future,resp_fun))
    test_Lost_future[is.na(test_Lost_future)]<- dat_resp_n[order(dat_resp_n[,1]),][nrow(dat_resp_n),2]
    
    future_ssp126_Lost_past<-cbind.data.frame(future_ssp126_Lost_past,test_Lost_past)
    future_ssp126_Lost_future<-cbind.data.frame(future_ssp126_Lost_future,test_Lost_future)
    
    
    Ele_colonized<-future_ssp126_Colonized
    Ele_colonized["Ele"]<-extract(elmat,Ele_colonized[,1:2])
    Colonized_low_ssp126<-Ele_colonized[Ele_colonized$Ele<limit,]
    Colonized_high_ssp126<-Ele_colonized[Ele_colonized$Ele>=limit,]
    
    ssp126_Colonized_past<-raster::extract(env_1970_2000.Reduced[[n]],future_ssp126_Colonized[,1:2], sp = T) #mean variable value at Colonizeded habitat in the past
    ssp126_Colonized_future<-raster::extract(env_2081_2100_ssp126_Reduce[[n]],future_ssp126_Colonized[,1:2], sp = T) #mean variable value at Colonizeded habitat in the future
 
    
    ssp126_Colonized_past_Low<-raster::extract(env_1970_2000.Reduced[[n]],Colonized_low_ssp126[,1:2], sp = T) #mean variable value at Colonizeded habitat in the past
    ssp126_Colonized_future_Low<-raster::extract(env_2081_2100_ssp126_Reduce[[n]],Colonized_low_ssp126[,1:2], sp = T) #mean variable value at Colonizeded habitat in the future
    ssp126change_Colonized_Low[n]<-ifelse((mean(na.omit(ssp126_Colonized_future_Low))-mean(na.omit(ssp126_Colonized_past_Low)))>0,1,-1) #increase=1, decrease=-1
    
    ssp126_Colonized_past_High<-raster::extract(env_1970_2000.Reduced[[n]],Colonized_high_ssp126[,1:2], sp = T) #mean variable value at Colonizeded habitat in the past
    ssp126_Colonized_future_High<-raster::extract(env_2081_2100_ssp126_Reduce[[n]],Colonized_high_ssp126[,1:2], sp = T) #mean variable value at Colonizeded habitat in the future
    ssp126change_Colonized_High[n]<-ifelse((mean(na.omit(ssp126_Colonized_future_High))-mean(na.omit(ssp126_Colonized_past_High)))>0,1,-1) #increase=1, decrease=-1
    
    
    
    future_ssp126_Colonized["current"]<-raster::extract(env_1970_2000.Reduced[[n]],future_ssp126_Colonized[,1:2], sp = T) #mean variable value at Colonized habitat in the past
    future_ssp126_Colonized["ssp126"]<-raster::extract(env_2081_2100_ssp126_Reduce[[n]],future_ssp126_Colonized[,1:2], sp = T) #mean variable value at Colonized habitat in the future
    future_ssp126_Colonized["Traits"]<-varnames[n]
    ssp126_Colonized<-melt(future_ssp126_Colonized,id.vars = c("Longitude","Latitude","Traits"))
    ssp126_Colonized_result<-rbind(ssp126_Colonized_result,ssp126_Colonized)
    
    
    test_Colonized_past<-unlist(lapply(ssp126_Colonized_past,resp_fun))
    test_Colonized_past[is.na(test_Colonized_past)]<- dat_resp_n[order(dat_resp_n[,1]),][nrow(dat_resp_n),2]
    test_Colonized_future<-unlist(lapply(ssp126_Colonized_future,resp_fun))
    test_Colonized_future[is.na(test_Colonized_future)]<- dat_resp_n[order(dat_resp_n[,1]),][nrow(dat_resp_n),2]
    
    future_ssp126_Colonized_past<-cbind.data.frame(future_ssp126_Colonized_past,test_Colonized_past)
    future_ssp126_Colonized_future<-cbind.data.frame(future_ssp126_Colonized_future,test_Colonized_future)
  }
  
  colnames(future_ssp126_Lost_past)<-c("Longitude","Latitude",varnames)
  colnames(future_ssp126_Lost_future)<-c("Longitude","Latitude",varnames)
  
  
  future_ssp126_Lost_change<-future_ssp126_Lost_future[,3:8]-future_ssp126_Lost_past[,3:8]
  future_ssp126_Lost_change1<-melt(future_ssp126_Lost_change)
  future_ssp126_Lost_change2<-future_ssp126_Lost_change1 %>%
    group_by(variable) %>%
    #summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))
    summarise(mean=mean(value/1000,na.rm=T), sd=sd(value/1000,na.rm=T),se=sd/sqrt(n))
  future_ssp126_Lost_change2<-as.data.frame(future_ssp126_Lost_change2)
  
  colnames(future_ssp126_Colonized_past)<-c("Longitude","Latitude",varnames)
  colnames(future_ssp126_Colonized_future)<-c("Longitude","Latitude",varnames)
  future_ssp126_Colonized_past["Elevation"]<-extract(elmat,future_ssp126_Colonized_past[,c(1:2)])
  future_ssp126_Colonized_future["Elevation"]<-extract(elmat,future_ssp126_Colonized_future[,c(1:2)])
  
  
  future_ssp126_Colonized_change_Low<-future_ssp126_Colonized_future[future_ssp126_Colonized_future$Elevation<limit,][,3:8]-future_ssp126_Colonized_past[future_ssp126_Colonized_past$Elevation<limit,][,3:8]
  future_ssp126_Colonized_change_High<-future_ssp126_Colonized_future[future_ssp126_Colonized_future$Elevation>=limit,][,3:8]-future_ssp126_Colonized_past[future_ssp126_Colonized_past$Elevation>=limit,][,3:8]
  
   future_ssp126_Colonized_change1_Low<-melt(future_ssp126_Colonized_change_Low)
  future_ssp126_Colonized_change2_Low<-future_ssp126_Colonized_change1_Low %>%
    group_by(variable) %>%
    #summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))
    summarise(mean=mean(value/1000,na.rm=T), sd=sd(value/1000,na.rm=T),se=sd/sqrt(n))
  future_ssp126_Colonized_change2_Low<-as.data.frame(future_ssp126_Colonized_change2_Low)
  
  future_ssp126_Colonized_change1_High<-melt(future_ssp126_Colonized_change_High)
  future_ssp126_Colonized_change2_High<-future_ssp126_Colonized_change1_High %>%
    group_by(variable) %>%
    #summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))
    summarise(mean=mean(value/1000,na.rm=T), sd=sd(value/1000,na.rm=T),se=sd/sqrt(n))
  future_ssp126_Colonized_change2_High<-as.data.frame(future_ssp126_Colonized_change2_High)
  
  
  
  
  future_ssp126_change<-rbind(future_ssp126_Lost_change2,future_ssp126_Colonized_change2_High,future_ssp126_Colonized_change2_Low)
  future_ssp126_change$Population<-rep(pop,18)
  future_ssp126_change$Habitat<-c(rep("Lost",6),rep("Colonized_High",6),rep("Colonized_Low",6))
  future_ssp126_change$varchange<-c(ssp126change_Lost,ssp126change_Colonized_High,ssp126change_Colonized_Low)
  
  ssp126_change<-rbind.data.frame(ssp126_change,future_ssp126_change)  
  
  
  
  
  #####SSP585
  ##Occurrence probability induced by each ssp585 no adaptation assuming all other variables are constant(mean value) at Lost/Colonizeded suitable habitat (induced by only LUC)
  future_ssp585_Lost_past<-future_ssp585_Lost
  future_ssp585_Lost_future<-future_ssp585_Lost
  future_ssp585_Colonized_past<-future_ssp585_Colonized
  future_ssp585_Colonized_future<-future_ssp585_Colonized
  
  
  ssp585change_Lost<-vector()
  ssp585change_Colonized_Low<-vector()
  ssp585change_Colonized_High<-vector()
  
  ssp585_Lost_result<-data.frame()
  ssp585_Colonized_result<-data.frame()
  
  
  for(n in 1:6){
    
   # n=4
    dat_resp_n<-dat_resp[,c((n*5-4),n*5)]
    resp_fun<-approxfun(dat_resp_n[,1],dat_resp_n[,2]) #Linear interpolation
    
    ssp585_Lost_past<-raster::extract(env_1970_2000.Reduced[[n]],future_ssp585_Lost[,1:2], sp = T) #mean variable value at Lost habitat in the past
    ssp585_Lost_future<-raster::extract(env_2081_2100_ssp585_Reduce[[n]],future_ssp585_Lost[,1:2], sp = T) #mean variable value at Lost habitat in the future
    ssp585change_Lost[n]<-ifelse((mean(na.omit(ssp585_Lost_future))-mean(na.omit(ssp585_Lost_past)))>0,1,-1) #increase=1, decrease=-1
    
    future_ssp585_Lost["current"]<-raster::extract(env_1970_2000.Reduced[[n]],future_ssp585_Lost[,1:2], sp = T) #mean variable value at Lost habitat in the past
    future_ssp585_Lost["ssp585"]<-raster::extract(env_2081_2100_ssp585_Reduce[[n]],future_ssp585_Lost[,1:2], sp = T) #mean variable value at Lost habitat in the future
    future_ssp585_Lost["Traits"]<-varnames[n]
    ssp585_Lost<-melt(future_ssp585_Lost,id.vars = c("Longitude","Latitude","Traits"))
    ssp585_Lost_result<-rbind(ssp585_Lost_result,ssp585_Lost)
    
    
    ###replace the na by the last value
    test_Lost_past<-unlist(lapply(ssp585_Lost_past,resp_fun))
    test_Lost_past[is.na(test_Lost_past)]<- dat_resp_n[order(dat_resp_n[,1]),][nrow(dat_resp_n),2]
    test_Lost_future<-unlist(lapply(ssp585_Lost_future,resp_fun))
    test_Lost_future[is.na(test_Lost_future)]<- dat_resp_n[order(dat_resp_n[,1]),][nrow(dat_resp_n),2]
    
    future_ssp585_Lost_past<-cbind.data.frame(future_ssp585_Lost_past,test_Lost_past)
    future_ssp585_Lost_future<-cbind.data.frame(future_ssp585_Lost_future,test_Lost_future)
    
    Ele_colonized<-future_ssp585_Colonized
    Ele_colonized["Ele"]<-extract(elmat,Ele_colonized[,1:2])
    Colonized_low_ssp585<-Ele_colonized[Ele_colonized$Ele<limit,]
    Colonized_high_ssp585<-Ele_colonized[Ele_colonized$Ele>=limit,]
    
    ssp585_Colonized_past<-raster::extract(env_1970_2000.Reduced[[n]],future_ssp585_Colonized[,1:2], sp = T) #mean variable value at Colonizeded habitat in the past
    ssp585_Colonized_future<-raster::extract(env_2081_2100_ssp585_Reduce[[n]],future_ssp585_Colonized[,1:2], sp = T) #mean variable value at Colonizeded habitat in the future
    
    
    ssp585_Colonized_past_Low<-raster::extract(env_1970_2000.Reduced[[n]],Colonized_low_ssp585[,1:2], sp = T) #mean variable value at Colonizeded habitat in the past
    ssp585_Colonized_future_Low<-raster::extract(env_2081_2100_ssp585_Reduce[[n]],Colonized_low_ssp585[,1:2], sp = T) #mean variable value at Colonizeded habitat in the future
    ssp585change_Colonized_Low[n]<-ifelse((mean(na.omit(ssp585_Colonized_future_Low))-mean(na.omit(ssp585_Colonized_past_Low)))>0,1,-1) #increase=1, decrease=-1
    
    ssp585_Colonized_past_High<-raster::extract(env_1970_2000.Reduced[[n]],Colonized_high_ssp585[,1:2], sp = T) #mean variable value at Colonizeded habitat in the past
    ssp585_Colonized_future_High<-raster::extract(env_2081_2100_ssp585_Reduce[[n]],Colonized_high_ssp585[,1:2], sp = T) #mean variable value at Colonizeded habitat in the future
    ssp585change_Colonized_High[n]<-ifelse((mean(na.omit(ssp585_Colonized_future_High))-mean(na.omit(ssp585_Colonized_past_High)))>0,1,-1) #increase=1, decrease=-1
    
    
    future_ssp585_Colonized["current"]<-raster::extract(env_1970_2000.Reduced[[n]],future_ssp585_Colonized[,1:2], sp = T) #mean variable value at Colonized habitat in the past
    future_ssp585_Colonized["ssp585"]<-raster::extract(env_2081_2100_ssp585_Reduce[[n]],future_ssp585_Colonized[,1:2], sp = T) #mean variable value at Colonized habitat in the future
    future_ssp585_Colonized["Traits"]<-varnames[n]
    ssp585_Colonized<-melt(future_ssp585_Colonized,id.vars = c("Longitude","Latitude","Traits"))
    ssp585_Colonized_result<-rbind(ssp585_Colonized_result,ssp585_Colonized)
  
    
    test_Colonized_past<-unlist(lapply(ssp585_Colonized_past,resp_fun))
    test_Colonized_past[is.na(test_Colonized_past)]<- dat_resp_n[order(dat_resp_n[,1]),][nrow(dat_resp_n),2]
    test_Colonized_future<-unlist(lapply(ssp585_Colonized_future,resp_fun))
    test_Colonized_future[is.na(test_Colonized_future)]<- dat_resp_n[order(dat_resp_n[,1]),][nrow(dat_resp_n),2]
    
    future_ssp585_Colonized_past<-cbind.data.frame(future_ssp585_Colonized_past,test_Colonized_past)
    future_ssp585_Colonized_future<-cbind.data.frame(future_ssp585_Colonized_future,test_Colonized_future)
  }
  
  colnames(future_ssp585_Lost_past)<-c("Longitude","Latitude",varnames)
  colnames(future_ssp585_Lost_future)<-c("Longitude","Latitude",varnames)
  
  
  future_ssp585_Lost_change<-future_ssp585_Lost_future[,3:8]-future_ssp585_Lost_past[,3:8]
  future_ssp585_Lost_change1<-melt(future_ssp585_Lost_change)
  future_ssp585_Lost_change2<-future_ssp585_Lost_change1 %>%
    group_by(variable) %>%
    #summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))
    summarise(mean=mean(value/1000,na.rm=T), sd=sd(value/1000,na.rm=T),se=sd/sqrt(n))
  future_ssp585_Lost_change2<-as.data.frame(future_ssp585_Lost_change2)
  
  colnames(future_ssp585_Colonized_past)<-c("Longitude","Latitude",varnames)
  colnames(future_ssp585_Colonized_future)<-c("Longitude","Latitude",varnames)
  future_ssp585_Colonized_past["Elevation"]<-extract(elmat,future_ssp585_Colonized_past[,c(1:2)])
  future_ssp585_Colonized_future["Elevation"]<-extract(elmat,future_ssp585_Colonized_future[,c(1:2)])
  
  
  future_ssp585_Colonized_change_Low<-future_ssp585_Colonized_future[future_ssp585_Colonized_future$Elevation<limit,][,3:8]-future_ssp585_Colonized_past[future_ssp585_Colonized_past$Elevation<limit,][,3:8]
  future_ssp585_Colonized_change_High<-future_ssp585_Colonized_future[future_ssp585_Colonized_future$Elevation>=limit,][,3:8]-future_ssp585_Colonized_past[future_ssp585_Colonized_past$Elevation>=limit,][,3:8]
  
  future_ssp585_Colonized_change1_Low<-melt(future_ssp585_Colonized_change_Low)
  future_ssp585_Colonized_change2_Low<-future_ssp585_Colonized_change1_Low %>%
    group_by(variable) %>%
    #summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))
    summarise(mean=mean(value/1000,na.rm=T), sd=sd(value/1000,na.rm=T),se=sd/sqrt(n))
  future_ssp585_Colonized_change2_Low<-as.data.frame(future_ssp585_Colonized_change2_Low)
  
  future_ssp585_Colonized_change1_High<-melt(future_ssp585_Colonized_change_High)
  future_ssp585_Colonized_change2_High<-future_ssp585_Colonized_change1_High %>%
    group_by(variable) %>%
    #summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))
    summarise(mean=mean(value/1000,na.rm=T), sd=sd(value/1000,na.rm=T),se=sd/sqrt(n))
  future_ssp585_Colonized_change2_High<-as.data.frame(future_ssp585_Colonized_change2_High)
  
  
  
  
  future_ssp585_change<-rbind(future_ssp585_Lost_change2,future_ssp585_Colonized_change2_High,future_ssp585_Colonized_change2_Low)
  future_ssp585_change$Population<-rep(pop,18)
  future_ssp585_change$Habitat<-c(rep("Lost",6),rep("Colonized_High",6),rep("Colonized_Low",6))
  future_ssp585_change$varchange<-c(ssp585change_Lost,ssp585change_Colonized_High,ssp585change_Colonized_Low)
  
  ssp585_change<-rbind.data.frame(ssp585_change,future_ssp585_change)  
  
  
  ssp126_Lost_result["Pop"]<-pop
  ssp126_Colonized_result["Pop"]<-pop
  ssp585_Lost_result["Pop"]<-pop
  ssp585_Colonized_result["Pop"]<-pop
  
  ssp126_Lost_result["ssp"]<-"ssp126"
  ssp126_Colonized_result["ssp"]<-"ssp126"
  ssp585_Lost_result["ssp"]<-"ssp585"
  ssp585_Colonized_result["ssp"]<-"ssp585"
  
  ssp126_Lost_result["class"]<-"Lost"
  ssp126_Colonized_result["class"]<-"Colonized"
  ssp585_Lost_result["class"]<-"Lost"
  ssp585_Colonized_result["class"]<-"Colonized"
  
  
  Change_result<-rbind(ssp126_Lost_result,ssp126_Colonized_result,ssp585_Lost_result,ssp585_Colonized_result)
  
  all_changeresult<-rbind(all_changeresult,Change_result)
  
  print(pop)
}


all_changeresult["Elevation"]<-extract(elmat,all_changeresult[,c(1:2)])


all_low<-all_changeresult[which(all_changeresult$Elevation<limit),]
all_high<-all_changeresult[which(all_changeresult$Elevation>=limit),]
all_low["Habitats"]<-"Low_colonized"
all_high["Habitats"]<-"High_colonized"
all_result<-rbind(all_low,all_high)
all_result[which(all_result$class=="Lost"),]$Habitats<-"Lost"


write.csv(all_result,file = paste0("E:/QT/SDM/",limit,"m/SDM_Group/change_traits.csv"),row.names=F)

ssp126_change["ssp"]<-"ssp126"
ssp585_change["ssp"]<-"ssp585"

Habitat_change<-rbind(ssp126_change,ssp585_change)
write.csv(Habitat_change,paste0(path_rds,"Habchan_varcontribution.csv"),row.names=F)


#####plot Figure S3


Low_habchange<-Habitat_change[Habitat_change$Population=="Low",]
Low_habchange_ssp126<-Low_habchange[Low_habchange$ssp=="ssp126",]
Low_Lost_ssp126<-Low_habchange_ssp126[Low_habchange_ssp126$Habitat =="Lost",]
Low_Colonized_ssp126_Low<-Low_habchange_ssp126[Low_habchange_ssp126$Habitat =="Colonized_Low",]
Low_Colonized_ssp126_High<-Low_habchange_ssp126[Low_habchange_ssp126$Habitat =="Colonized_High",]

Low_habchange_ssp585<-Low_habchange[Low_habchange$ssp=="ssp585",]
Low_Lost_ssp585<-Low_habchange_ssp585[Low_habchange_ssp585$Habitat =="Lost",]
Low_Colonized_ssp585_Low<-Low_habchange_ssp585[Low_habchange_ssp585$Habitat =="Colonized_Low",]
Low_Colonized_ssp585_High<-Low_habchange_ssp585[Low_habchange_ssp585$Habitat =="Colonized_High",]


High_habchange<-Habitat_change[Habitat_change$Population=="High",]
High_habchange_ssp126<-High_habchange[High_habchange$ssp=="ssp126",]
High_Lost_ssp126<-High_habchange_ssp126[High_habchange_ssp126$Habitat =="Lost",]
High_Colonized_ssp126_Low<-High_habchange_ssp126[High_habchange_ssp126$Habitat =="Colonized_Low",]
High_Colonized_ssp126_High<-High_habchange_ssp126[High_habchange_ssp126$Habitat =="Colonized_High",]


High_habchange_ssp585<-High_habchange[High_habchange$ssp=="ssp585",]
High_Lost_ssp585<-High_habchange_ssp585[High_habchange_ssp585$Habitat =="Lost",]
High_Colonized_ssp585_Low<-High_habchange_ssp585[High_habchange_ssp585$Habitat =="Colonized_Low",]
High_Colonized_ssp585_High<-High_habchange_ssp585[High_habchange_ssp585$Habitat =="Colonized_High",]







#####Low elevation#### 
Fig_Low_Lost_ssp126<-ggplot(Low_Lost_ssp126)+ #,aes(Var,Change,fill = Change < 0)  #position=pd
  geom_pointrange(aes(y=mean,ymin=mean-se, ymax=mean+se,x=variable,color=factor(varchange)),size=0.3)+
  scale_colour_manual(values=c("#D55E00",'cornflowerblue'))+
  #scale_y_continuous(trans = 'signed')+
  #coord_trans(y=logx_trans())+
  #annotation_logticks(sides="b")+  
  geom_hline(yintercept=0,linetype="dashed", color = "grey",size=0.9)+
  ylab("Change of probability")+
  #geom_bar(stat="identity",width=0.6)+
  coord_flip()+ 
  theme_light()+
  #ylab("Change of probability")+
  xlab("Low_Elevation
  Variables")+
 # labs(title = 'Lost_ssp126')+
theme(strip.text = element_blank(),panel.grid.major = element_blank(),
      plot.title = element_text(size=18,hjust=0.5,color = 'black'),
      panel.grid.minor = element_blank(),legend.position = "none",
      axis.title=element_text(angle=0,size=16),
      #axis.title.y=element_blank(),
      #axis.title.x=element_blank(),
      axis.text.y=element_text(size=12),
      axis.text.x=element_text(size=12),
      plot.margin = unit(c(1,0.2,1,1), "lines"))



Fig_Low_Colonized_ssp126_Low<-ggplot(Low_Colonized_ssp126_Low)+ #,aes(Var,Change,fill = Change < 0)  #position=pd
  geom_pointrange(aes(y=mean,ymin=mean-se, ymax=mean+se,x=variable,color=factor(varchange)),size=0.3)+
  scale_colour_manual(values=c("#DC143C","#4169E1"))+
  #scale_y_continuous(trans = 'signed')+
  #coord_trans(y=logx_trans())+
  #annotation_logticks(sides="b")+  
  geom_hline(yintercept=0,linetype="dashed", color = "grey",size=0.9)+
  
  #geom_bar(stat="identity",width=0.6)+
  coord_flip()+ 
  theme_light()+
   ylab("Change of probability")+
  # xlab("Variables")+
 # labs(title = 'Colonized_ssp126')+
  theme(strip.text = element_blank(),panel.grid.major = element_blank(),
        plot.title = element_text(size=18,hjust=0.5,color = 'black'),
        panel.grid.minor = element_blank(),legend.position = "none",
        axis.title=element_text(angle=0,size=16),
        axis.title.y=element_blank(),
        #axis.title.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        plot.margin = unit(c(1,0.2,1,1), "lines"))

Fig_Low_Colonized_ssp126_High<-ggplot(Low_Colonized_ssp126_High)+ #,aes(Var,Change,fill = Change < 0)  #position=pd
  geom_pointrange(aes(y=mean,ymin=mean-se, ymax=mean+se,x=variable,color=factor(varchange)),size=0.3)+
  scale_colour_manual(values=c("#DC143C","#4169E1"))+
  #scale_y_continuous(trans = 'signed')+
  #coord_trans(y=logx_trans())+
  #annotation_logticks(sides="b")+  
  geom_hline(yintercept=0,linetype="dashed", color = "grey",size=0.9)+
  
  #geom_bar(stat="identity",width=0.6)+
  coord_flip()+ 
  theme_light()+
  ylab("Change of probability")+
  # xlab("Variables")+
  # labs(title = 'Colonized_ssp126')+
  theme(strip.text = element_blank(),panel.grid.major = element_blank(),
        plot.title = element_text(size=18,hjust=0.5,color = 'black'),
        panel.grid.minor = element_blank(),legend.position = "none",
        axis.title=element_text(angle=0,size=16),
        axis.title.y=element_blank(),
        #axis.title.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        plot.margin = unit(c(1,0.2,1,1), "lines"))



Fig_Low_Lost_ssp585<-ggplot(Low_Lost_ssp585)+ #,aes(Var,Change,fill = Change < 0)  #position=pd
  geom_pointrange(aes(y=mean,ymin=mean-se, ymax=mean+se,x=variable,color=factor(varchange)),size=0.3)+
  scale_colour_manual(values=c("#DC143C","#4169E1"))+
  #scale_y_continuous(trans = 'signed')+
  #coord_trans(y=logx_trans())+
  #annotation_logticks(sides="b")+  
  geom_hline(yintercept=0,linetype="dashed", color = "grey",size=0.9)+
  
  #geom_bar(stat="identity",width=0.6)+
  coord_flip()+ 
  theme_light()+
  ylab("Change of probability")+
  xlab("Variables")+
 # labs(title = 'Lost_ssp585')+
  theme(strip.text = element_blank(),panel.grid.major = element_blank(),
        plot.title = element_text(size=18,hjust=0.5,color = 'black'),
        panel.grid.minor = element_blank(),legend.position = "none",
        axis.title=element_text(angle=0,size=16),
        axis.title.y=element_blank(),
        #axis.title.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        plot.margin = unit(c(1,0.2,1,1), "lines"))

Fig_Low_Colonized_ssp585_Low<-ggplot(Low_Colonized_ssp585_Low)+
  geom_pointrange(aes(y=mean,ymin=mean-se, ymax=mean+se,x=variable,color=factor(varchange)), size=0.3) +#position=pd
  scale_colour_manual(values=c("#DC143C","#4169E1"))+
  geom_hline(yintercept=0,linetype="dashed", color = "grey",size=0.9)+
  coord_flip()+ 
  theme_light()+
  ylab("Change of probability")+
  #xlab("Climate variables")+
#  labs(title = 'Colonized_ssp585')+
  theme(strip.text = element_blank(),panel.grid.major = element_blank(),
        plot.title = element_text(size=18,hjust=0.5,color = 'black'),
        panel.grid.minor = element_blank(),legend.position = "none",
        axis.title=element_text(angle=0,size=16),
        axis.title.y=element_blank(),
        #axis.title.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        plot.margin = unit(c(1,0.2,1,1), "lines"))

Fig_Low_Colonized_ssp585_High<-ggplot(Low_Colonized_ssp585_High)+
  geom_pointrange(aes(y=mean,ymin=mean-se, ymax=mean+se,x=variable,color=factor(varchange)), size=0.3) +#position=pd
  scale_colour_manual(values=c("#DC143C","#4169E1"))+
  geom_hline(yintercept=0,linetype="dashed", color = "grey",size=0.9)+
  coord_flip()+ 
  theme_light()+
  ylab("Change of probability")+
  #xlab("Climate variables")+
  #  labs(title = 'Colonized_ssp585')+
  theme(strip.text = element_blank(),panel.grid.major = element_blank(),
        plot.title = element_text(size=18,hjust=0.5,color = 'black'),
        panel.grid.minor = element_blank(),legend.position = "none",
        axis.title=element_text(angle=0,size=16),
        axis.title.y=element_blank(),
        #axis.title.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        plot.margin = unit(c(1,0.2,1,1), "lines"))


#####High elevation#####


Fig_High_Lost_ssp126<-ggplot(High_Lost_ssp126)+ #,aes(Var,Change,fill = Change < 0)  #position=pd
  geom_pointrange(aes(y=mean,ymin=mean-se, ymax=mean+se,x=variable,color=factor(varchange)),size=0.3)+
  scale_colour_manual(values=c("#DC143C","#4169E1"))+
  #scale_y_continuous(trans = 'signed')+
  #coord_trans(y=logx_trans())+
  #annotation_logticks(sides="b")+  
  geom_hline(yintercept=0,linetype="dashed", color = "grey",size=0.9)+
  
  #geom_bar(stat="identity",width=0.6)+
  coord_flip()+ 
  theme_light()+
 # ylab("Change of probability")+
  xlab("High_Elevation
  Variables")+
  labs(title = 'Lost_ssp126')+
  theme(strip.text = element_blank(),panel.grid.major = element_blank(),
        plot.title = element_text(size=18,hjust=0.5,color = 'black'),
        panel.grid.minor = element_blank(),legend.position = "none",
        axis.title=element_text(angle=0,size=16),
        #axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        plot.margin = unit(c(1,0.2,1,1), "lines"))

Fig_High_Colonized_ssp126_Low<-ggplot(High_Colonized_ssp126_Low)+ #,aes(Var,Change,fill = Change < 0)  #position=pd
  geom_pointrange(aes(y=mean,ymin=mean-se, ymax=mean+se,x=variable,color=factor(varchange)),size=0.3)+
  scale_colour_manual(values=c("#DC143C","#4169E1"))+
  #scale_y_continuous(trans = 'signed')+
  #coord_trans(y=logx_trans())+
  #annotation_logticks(sides="b")+  
  geom_hline(yintercept=0,linetype="dashed", color = "grey",size=0.9)+
  
  #geom_bar(stat="identity",width=0.6)+
  coord_flip()+ 
  theme_light()+
  labs(title = 'Colonized_ssp126_Low')+
  #ylab("Change of probability")+
  xlab("Variables")+
  #labs(title = 'High_Colonized_ssp126')+
  theme(strip.text = element_blank(),panel.grid.major = element_blank(),
        plot.title = element_text(size=18,hjust=0.5,color = 'black'),
        panel.grid.minor = element_blank(),legend.position = "none",
        axis.title=element_text(angle=0,size=16),
        axis.text.x=element_text(size=12),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        plot.margin = unit(c(1,0.2,1,1), "lines"))

Fig_High_Colonized_ssp126_High<-ggplot(High_Colonized_ssp126_High)+ #,aes(Var,Change,fill = Change < 0)  #position=pd
  geom_pointrange(aes(y=mean,ymin=mean-se, ymax=mean+se,x=variable,color=factor(varchange)),size=0.3)+
  scale_colour_manual(values=c("#DC143C","#4169E1"))+
  #scale_y_continuous(trans = 'signed')+
  #coord_trans(y=logx_trans())+
  #annotation_logticks(sides="b")+  
  geom_hline(yintercept=0,linetype="dashed", color = "grey",size=0.9)+
  
  #geom_bar(stat="identity",width=0.6)+
  coord_flip()+ 
  theme_light()+
  labs(title = 'Colonized_ssp126_High')+
  #ylab("Change of probability")+
  xlab("Variables")+
  #labs(title = 'High_Colonized_ssp126')+
  theme(strip.text = element_blank(),panel.grid.major = element_blank(),
        plot.title = element_text(size=18,hjust=0.5,color = 'black'),
        panel.grid.minor = element_blank(),legend.position = "none",
        axis.title=element_text(angle=0,size=16),
        axis.text.x=element_text(size=12),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        plot.margin = unit(c(1,0.2,1,1), "lines"))

Fig_High_Lost_ssp585<-ggplot(High_Lost_ssp585)+ #,aes(Var,Change,fill = Change < 0)  #position=pd
  geom_pointrange(aes(y=mean,ymin=mean-se, ymax=mean+se,x=variable,color=factor(varchange)),size=0.3)+
  scale_colour_manual(values=c("#DC143C","#4169E1"))+
  #scale_y_continuous(trans = 'signed')+
  #coord_trans(y=logx_trans())+
  #annotation_logticks(sides="b")+  
  geom_hline(yintercept=0,linetype="dashed", color = "grey",size=0.9)+
  labs(title = 'Lost_ssp585')+
  #geom_bar(stat="identity",width=0.6)+
  coord_flip()+ 
  theme_light()+
 # ylab("Change of probability")+
  xlab("Variables")+
  #labs(title = 'High_Lost_ssp585')+
  theme(strip.text = element_blank(),panel.grid.major = element_blank(),
        plot.title = element_text(size=18,hjust=0.5,color = 'black'),
        panel.grid.minor = element_blank(),legend.position = "none",
        axis.title=element_text(angle=0,size=16),
        axis.text.x=element_text(size=12),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        plot.margin = unit(c(1,0.2,1,1), "lines"))

Fig_High_Colonized_ssp585_Low<-ggplot(High_Colonized_ssp585_Low)+
  geom_pointrange(aes(y=mean,ymin=mean-se, ymax=mean+se,x=variable,color=factor(varchange)), size=0.3) +#position=pd
  scale_colour_manual(values=c("#DC143C","#4169E1"))+
  geom_hline(yintercept=0,linetype="dashed", color = "grey",size=0.9)+
  coord_flip()+ 
  theme_light()+
 # ylab("Change of probability")+
  labs(title = 'Colonized_ssp585_Low')+
  #xlab("Climate variables")+
  #labs(title = 'High_Colonized_ssp585')+
  theme(strip.text = element_blank(),panel.grid.major = element_blank(),
        plot.title = element_text(size=18,hjust=0.5,color = 'black'),
        panel.grid.minor = element_blank(),legend.position = "none",
        axis.title=element_text(angle=0,size=16),
        axis.text.x=element_text(size=12),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        plot.margin = unit(c(1,0.2,1,1), "lines"))

Fig_High_Colonized_ssp585_High<-ggplot(High_Colonized_ssp585_High)+
  geom_pointrange(aes(y=mean,ymin=mean-se, ymax=mean+se,x=variable,color=factor(varchange)), size=0.3) +#position=pd
  scale_colour_manual(values=c("#DC143C","#4169E1"))+
  geom_hline(yintercept=0,linetype="dashed", color = "grey",size=0.9)+
  coord_flip()+ 
  theme_light()+
  # ylab("Change of probability")+
  labs(title = 'Colonized_ssp585_High')+
  #xlab("Climate variables")+
  #labs(title = 'High_Colonized_ssp585')+
  theme(strip.text = element_blank(),panel.grid.major = element_blank(),
        plot.title = element_text(size=18,hjust=0.5,color = 'black'),
        panel.grid.minor = element_blank(),legend.position = "none",
        axis.title=element_text(angle=0,size=16),
        axis.text.x=element_text(size=12),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        plot.margin = unit(c(1,0.2,1,1), "lines"))


library(ggpubr)

pdf(paste(path_rds,"Fig_S3.pdf",sep=""),height=15,width=20)

ggarrange(Fig_High_Lost_ssp126,Fig_High_Colonized_ssp126_Low,Fig_High_Colonized_ssp126_High,Fig_High_Lost_ssp585,Fig_High_Colonized_ssp585_Low,Fig_High_Colonized_ssp585_High,Fig_Low_Lost_ssp126,Fig_Low_Colonized_ssp126_Low,Fig_Low_Colonized_ssp126_High,Fig_Low_Lost_ssp585,Fig_Low_Colonized_ssp585_Low,Fig_Low_Colonized_ssp585_High,ncol=6,nrow=2,labels=c('(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)',"(I)","(J)","(K)","(M)","(L)","(N)"),
          font.label = list(size = 18, color = "black"))

dev.off()























#####Boxplot of change variable in the change habitats for Response curve####


change_var<-read.csv(paste0("E:/QT/SDM/",limit,"m/SDM_Group/change_traits.csv"))


change_var$Traits<-factor(change_var$Traits,levels = c("Bio3","Bio4","Bio14","Bio15","OC","WL"))
change_var$Habitats<-factor(change_var$Habitats,levels = c("Low_colonized","High_colonized","Lost"))

Change_var_Low<-change_var[change_var$Pop=="Low",]
Change_var_High<-change_var[change_var$Pop=="High",]


mean(Change_var_Low[Change_var_Low$Traits=="WL"&Change_var_Low$variable=="current"&Change_var_Low$ssp=="ssp126"&Change_var_Low$class=="Lost",]$value)
mean(Change_var_Low[Change_var_Low$Traits=="WL"&Change_var_Low$variable=="current"&Change_var_Low$ssp=="ssp585"&Change_var_Low$class=="Lost",]$value)
mean(Change_var_Low[Change_var_Low$Traits=="WL"&Change_var_Low$variable=="ssp126"&Change_var_Low$ssp=="ssp126"&Change_var_Low$class=="Lost",]$value)
mean(Change_var_Low[Change_var_Low$Traits=="WL"&Change_var_Low$variable=="ssp585"&Change_var_Low$ssp=="ssp585"&Change_var_Low$class=="Lost",]$value)


mean(Change_var_Low[Change_var_Low$Traits=="WL"&Change_var_Low$variable=="current"&Change_var_Low$ssp=="ssp126"&Change_var_Low$class=="Colonized",]$value)
mean(Change_var_Low[Change_var_Low$Traits=="WL"&Change_var_Low$variable=="current"&Change_var_Low$ssp=="ssp585"&Change_var_Low$class=="Colonized",]$value)
mean(Change_var_Low[Change_var_Low$Traits=="WL"&Change_var_Low$variable=="ssp126"&Change_var_Low$ssp=="ssp126"&Change_var_Low$class=="Colonized",]$value)
mean(Change_var_Low[Change_var_Low$Traits=="WL"&Change_var_Low$variable=="ssp585"&Change_var_Low$ssp=="ssp585"&Change_var_Low$class=="Colonized",]$value)






mean(Change_var_High[Change_var_High$Traits=="WL"&Change_var_High$variable=="current"&Change_var_High$ssp=="ssp126"&Change_var_High$class=="Lost",]$value)
mean(Change_var_High[Change_var_High$Traits=="WL"&Change_var_High$variable=="current"&Change_var_High$ssp=="ssp585"&Change_var_High$class=="Lost",]$value)
mean(Change_var_High[Change_var_High$Traits=="WL"&Change_var_High$variable=="ssp126"&Change_var_High$ssp=="ssp126"&Change_var_High$class=="Lost",]$value)
mean(Change_var_High[Change_var_High$Traits=="WL"&Change_var_High$variable=="ssp585"&Change_var_High$ssp=="ssp585"&Change_var_High$class=="Lost",]$value)


mean(Change_var_High[Change_var_High$Traits=="WL"&Change_var_High$variable=="current"&Change_var_High$ssp=="ssp126"&Change_var_High$class=="Colonized",]$value)
mean(Change_var_High[Change_var_High$Traits=="WL"&Change_var_High$variable=="current"&Change_var_High$ssp=="ssp585"&Change_var_High$class=="Colonized",]$value)
mean(Change_var_High[Change_var_High$Traits=="WL"&Change_var_High$variable=="ssp126"&Change_var_High$ssp=="ssp126"&Change_var_High$class=="Colonized",]$value)
mean(Change_var_High[Change_var_High$Traits=="WL"&Change_var_High$variable=="ssp585"&Change_var_High$ssp=="ssp585"&Change_var_High$class=="Colonized",]$value)


fig_chavar_Low<-ggplot(data=Change_var_Low,aes(x=variable,y=value,fill=Habitats))+
  scale_fill_manual(values=c("#24BEB4","#4169E1","#DC143C"))+
  geom_boxplot()+
 # facet_wrap(Traits~ssp,scales = "free",nrow=2)+
  facet_grid(Traits~ssp,scales = "free")+
  theme_classic()+
  xlab("Scenario")+
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

  
  ggsave(fig_chavar_Low,file=paste("E:/QT/SDM/",limit,"m/SDM_Group/Fig_S4_Low.pdf",sep=""),device = "pdf",height=9,width=9,dpi = 300)


  

  
  fig_chavar_High<-ggplot(data=Change_var_High,aes(x=variable,y=value,fill=Habitats))+
    scale_fill_manual(values=c("#24BEB4","#4169E1","#DC143C"))+
    geom_boxplot()+
    # facet_wrap(Traits~ssp,scales = "free",nrow=2)+
    facet_grid(Traits~ssp,scales = "free")+
    theme_classic()+
    xlab("Scenario")+
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
  
  
  ggsave(fig_chavar_High,file=paste("E:/QT/SDM/",limit,"m/SDM_Group/Fig_S4_High.pdf",sep=""),device = "pdf",height=9,width=9,dpi = 300)
  
  
  
}
  
  
  
  limit=3600
  #####model evaluation
  library(tidyr)
  Low_eva<-read.csv(paste0("E:/QT/SDM/",limit,"m/SDM_Group/Low/Phrynocephalus vlangalii_Low_EM_Eval_rep10.csv"))
  High_eva<-read.csv(paste0("E:/QT/SDM/",limit,"m/SDM_Group/High/Phrynocephalus vlangalii_High_EM_Eval_rep10.csv"))
  
  
  Low_eva<-separate(Low_eva,model,into = c("RUN","model"),sep = "_")
  High_eva<-separate(High_eva,model,into = c("RUN","model"),sep = "_")
  
  
  
  head(Low_eva)
  Low_eva_ANN_AUC<-mean(Low_eva[Low_eva$model=="ANN",]$AUC)
  Low_eva_ANN_TSS<-mean(Low_eva[Low_eva$model=="ANN",]$TSS)
  Low_eva_ANN_Boyce<-mean(Low_eva[Low_eva$model=="ANN",]$Boyce)
  Low_sd_ANN_AUC<-sd(Low_eva[Low_eva$model=="ANN",]$AUC)
  Low_sd_ANN_TSS<-sd(Low_eva[Low_eva$model=="ANN",]$TSS)
  Low_sd_ANN_Boyce<-sd(Low_eva[Low_eva$model=="ANN",]$Boyce)
  
  Low_eva_GLM_AUC<-mean(Low_eva[Low_eva$model=="GLM",]$AUC)
  Low_eva_GLM_TSS<-mean(Low_eva[Low_eva$model=="GLM",]$TSS)
  Low_eva_GLM_Boyce<-mean(Low_eva[Low_eva$model=="GLM",]$Boyce)
  Low_sd_GLM_AUC<-sd(Low_eva[Low_eva$model=="GLM",]$AUC)
  Low_sd_GLM_TSS<-sd(Low_eva[Low_eva$model=="GLM",]$TSS)
  Low_sd_GLM_Boyce<-sd(Low_eva[Low_eva$model=="GLM",]$Boyce)
    
  
  Low_eva_MAXENT.Phillips_AUC<-mean(Low_eva[Low_eva$model=="MAXENT.Phillips",]$AUC)
  Low_eva_MAXENT.Phillips_TSS<-mean(Low_eva[Low_eva$model=="MAXENT.Phillips",]$TSS)
  Low_eva_MAXENT.Phillips_Boyce<-mean(Low_eva[Low_eva$model=="MAXENT.Phillips",]$Boyce)
  Low_sd_MAXENT.Phillips_AUC<-sd(Low_eva[Low_eva$model=="MAXENT.Phillips",]$AUC)
  Low_sd_MAXENT.Phillips_TSS<-sd(Low_eva[Low_eva$model=="MAXENT.Phillips",]$TSS)
  Low_sd_MAXENT.Phillips_Boyce<-sd(Low_eva[Low_eva$model=="MAXENT.Phillips",]$Boyce)
  
  Low_eva_EF_AUC<-mean(Low_eva[Low_eva$model=="EF",]$AUC)
  Low_eva_EF_TSS<-mean(Low_eva[Low_eva$model=="EF",]$TSS)
  Low_eva_EF_Boyce<-mean(Low_eva[Low_eva$model=="EF",]$Boyce)
  Low_sd_EF_AUC<-sd(Low_eva[Low_eva$model=="EF",]$AUC)
  Low_sd_EF_TSS<-sd(Low_eva[Low_eva$model=="EF",]$TSS)
  Low_sd_EF_Boyce<-sd(Low_eva[Low_eva$model=="EF",]$Boyce)
  
  
  
  
  
  High_eva_ANN_AUC<-mean(High_eva[High_eva$model=="ANN",]$AUC)
  High_eva_ANN_TSS<-mean(High_eva[High_eva$model=="ANN",]$TSS)
  High_eva_ANN_Boyce<-mean(High_eva[High_eva$model=="ANN",]$Boyce)
  High_sd_ANN_AUC<-sd(High_eva[High_eva$model=="ANN",]$AUC)
  High_sd_ANN_TSS<-sd(High_eva[High_eva$model=="ANN",]$TSS)
  High_sd_ANN_Boyce<-sd(High_eva[High_eva$model=="ANN",]$Boyce)
  
  High_eva_GLM_AUC<-mean(High_eva[High_eva$model=="GLM",]$AUC)
  High_eva_GLM_TSS<-mean(High_eva[High_eva$model=="GLM",]$TSS)
  High_eva_GLM_Boyce<-mean(High_eva[High_eva$model=="GLM",]$Boyce)
  High_sd_GLM_AUC<-sd(High_eva[High_eva$model=="GLM",]$AUC)
  High_sd_GLM_TSS<-sd(High_eva[High_eva$model=="GLM",]$TSS)
  High_sd_GLM_Boyce<-sd(High_eva[High_eva$model=="GLM",]$Boyce)
  
  
  High_eva_MAXENT.Phillips_AUC<-mean(High_eva[High_eva$model=="MAXENT.Phillips",]$AUC)
  High_eva_MAXENT.Phillips_TSS<-mean(High_eva[High_eva$model=="MAXENT.Phillips",]$TSS)
  High_eva_MAXENT.Phillips_Boyce<-mean(High_eva[High_eva$model=="MAXENT.Phillips",]$Boyce)
  High_sd_MAXENT.Phillips_AUC<-sd(High_eva[High_eva$model=="MAXENT.Phillips",]$AUC)
  High_sd_MAXENT.Phillips_TSS<-sd(High_eva[High_eva$model=="MAXENT.Phillips",]$TSS)
  High_sd_MAXENT.Phillips_Boyce<-sd(High_eva[High_eva$model=="MAXENT.Phillips",]$Boyce)
  
  High_eva_EF_AUC<-mean(High_eva[High_eva$model=="EF",]$AUC)
  High_eva_EF_TSS<-mean(High_eva[High_eva$model=="EF",]$TSS)
  High_eva_EF_Boyce<-mean(High_eva[High_eva$model=="EF",]$Boyce)
  High_sd_EF_AUC<-sd(High_eva[High_eva$model=="EF",]$AUC)
  High_sd_EF_TSS<-sd(High_eva[High_eva$model=="EF",]$TSS)
  High_sd_EF_Boyce<-sd(High_eva[High_eva$model=="EF",]$Boyce)
  
  
  