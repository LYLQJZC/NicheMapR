##########Analysis for the results of ESM SDM###############
#####----Zhong-Wen Jiang----####
#####----10th March 2022----####

library(raster)
library(rgdal)
library(sp)
library(leaflet)
library(foreign)
library(ggplot2)
library(rgeos)


options(scipen=200)

Vlangalii_IUCN<-readOGR("E:/QT/SDM/P. vlangalii/data_P.shp")
IUCN_buffer<-gBuffer(Vlangalii_IUCN,width=0.5)
spe<-"Phrynocephalus.vlangalii"

#Elevation
elmat = raster("E:/QT/SDM/wc2.1_10m_elev.tif")
elmat=crop(elmat,c(79,110,24,43))
elmat.df<-data.frame(rasterToPoints(elmat))
colnames(elmat.df) <- c("lon", "lat", "ele")
Net_change_mean<-data.frame()




for (limit in c(3500,3600,3700)) {
  
limit=3600

  
  current_Low_pre<-raster(paste0("E:/QT/SDM/",limit,"m/SDM_Group/Low/TIFresult/binary_current_Low.tif"))
  current_High_pre<-raster(paste0("E:/QT/SDM/",limit,"m/SDM_Group/High/TIFresult/binary_current_High.tif"))
  
  current_crop_Low<-crop(current_Low_pre,IUCN_buffer)
  current_mask_Low<-mask(current_crop_Low,IUCN_buffer)
  
  current_crop_High<-crop(current_High_pre,IUCN_buffer)
  current_mask_High<-mask(current_crop_High,IUCN_buffer)
  
  ##########limit the current distributions of different population by elevation#####
  dat_current_mask_Low<-as.data.frame(rasterToPoints(current_mask_Low))
  colnames(dat_current_mask_Low)<-c("Longitude","Latitude","Preabs") #Preabs: presence (1) vs absence (0);
  dat_current_mask_Low["Elevation"]<-extract(elmat,dat_current_mask_Low[,c(1,2)])
  dat_current_mask_Low<-dat_current_mask_Low[dat_current_mask_Low$Elevation<limit,] 
  
  dat_current_mask_High<-as.data.frame(rasterToPoints(current_mask_High))
  colnames(dat_current_mask_High)<-c("Longitude","Latitude","Preabs") #Preabs: presence (1) vs absence (0);
  dat_current_mask_High["Elevation"]<-extract(elmat,dat_current_mask_High[,c(1,2)])
  dat_current_mask_High<-dat_current_mask_High[dat_current_mask_High$Elevation>=limit,]
  #######################
  
  current_raster_Low<-rasterFromXYZ(dat_current_mask_Low[dat_current_mask_Low$Preabs==1,]) 
  current_polygen_Low<-rasterToPolygons(current_raster_Low)
  future_buffer_Low<-gBuffer(current_polygen_Low,width=0.5)######the future dispersed boundary limited by predicted current with a 50km buffer n 
  
  current_raster_High<-rasterFromXYZ(dat_current_mask_High[dat_current_mask_High$Preabs==1,]) 
  current_polygen_High<-rasterToPolygons(current_raster_High)
  future_buffer_High<-gBuffer(current_polygen_High,width=0.5)######the future dispersed boundary limited by predicted current with a 50km buffer n 
  
  
  future_ssp126_Low_pre<-raster(paste0("E:/QT/SDM/",limit,"m/SDM_Group/Low/TIFresult/","binary_ssp126_wc_Low.tif"))
  future_ssp585_Low_pre<-raster(paste0("E:/QT/SDM/",limit,"m/SDM_Group/Low/TIFresult/","binary_ssp585_wc_Low.tif"))
  future_ssp126_High_pre<-raster(paste0("E:/QT/SDM/",limit,"m/SDM_Group/High/TIFresult/","binary_ssp126_wc_High.tif"))
  future_ssp585_High_pre<-raster(paste0("E:/QT/SDM/",limit,"m/SDM_Group/High/TIFresult/","binary_ssp585_wc_High.tif"))
  
  
  current_Low<-mask(current_Low_pre,future_buffer_Low)
  future_ssp126_Low<-mask(future_ssp126_Low_pre,future_buffer_Low)
  future_ssp585_Low<-mask(future_ssp585_Low_pre,future_buffer_Low)
  #####Elevation<limit =Low population,Elevation>=limit =High population#####
  future_ssp126_Low_ad1<-mask(future_ssp126_High_pre,future_buffer_Low)
  future_ssp585_Low_ad1<-mask(future_ssp585_High_pre,future_buffer_Low)
  ######
  
  current_High<-mask(current_High_pre,future_buffer_High)
  future_ssp126_High<-mask(future_ssp126_High_pre,future_buffer_High)
  future_ssp585_High<-mask(future_ssp585_High_pre,future_buffer_High)
  #####Elevation<limit =Low population,Elevation>=limit =High population#####
  future_ssp126_High_ad1<-mask(future_ssp126_Low_pre,future_buffer_High)
  future_ssp585_High_ad1<-mask(future_ssp585_Low_pre,future_buffer_High)
  ######
  
  
  #####Low elevation raster to point#########
  dat_current_Low<-as.data.frame(rasterToPoints(current_Low))
  colnames(dat_current_Low)<-c("Longitude","Latitude","Preabs") #Preabs: presence (1) vs absence (0);
  dat_ssp126_Low<-as.data.frame(rasterToPoints(future_ssp126_Low))
  colnames(dat_ssp126_Low)<-c("Longitude","Latitude","Preabs") #Preabs: presence (1) vs absence (0);
  dat_ssp585_Low<-as.data.frame(rasterToPoints(future_ssp585_Low))
  colnames(dat_ssp585_Low)<-c("Longitude","Latitude","Preabs") #Preabs: presence (1) vs absence (0);
  dat_ssp126_Low_ad1<-as.data.frame(rasterToPoints(future_ssp126_Low_ad1))
  colnames(dat_ssp126_Low_ad1)<-c("Longitude","Latitude","Preabs") #Preabs: presence (1) vs absence (0);
  dat_ssp585_Low_ad1<-as.data.frame(rasterToPoints(future_ssp585_Low_ad1))
  colnames(dat_ssp585_Low_ad1)<-c("Longitude","Latitude","Preabs") #Preabs: presence (1) vs absence (0);
  
  dat_current_Low["Elevation"]<-extract(elmat,dat_current_Low[,c(1,2)])
  dat_ssp126_Low["Elevation"]<-extract(elmat,dat_ssp126_Low[,c(1,2)])
  dat_ssp585_Low["Elevation"]<-extract(elmat,dat_ssp585_Low[,c(1,2)])
  dat_ssp126_Low_ad1["Elevation"]<-extract(elmat,dat_ssp126_Low_ad1[,c(1,2)])
  dat_ssp585_Low_ad1["Elevation"]<-extract(elmat,dat_ssp585_Low_ad1[,c(1,2)])
  
  dat_current_Low_col<-dat_current_Low
  dat_current_Low[dat_current_Low$Preabs==1&dat_current_Low$Elevation>=limit,]$Preabs<-0
  dat_ssp126_Low_ad<-rbind(dat_ssp126_Low[dat_ssp126_Low$Elevation<limit,],dat_ssp126_Low_ad1[dat_ssp126_Low_ad1$Elevation>=limit,])
  dat_ssp585_Low_ad<-rbind(dat_ssp585_Low[dat_ssp585_Low$Elevation<limit,],dat_ssp585_Low_ad1[dat_ssp585_Low_ad1$Elevation>=limit,])
  
  dat_current_Low_col<-dat_current_Low_col[order(dat_current_Low_col$Elevation),]
  dat_current_Low<-dat_current_Low[order(dat_current_Low$Elevation),]
  dat_ssp126_Low<-dat_ssp126_Low[order(dat_ssp126_Low$Elevation),]
  dat_ssp126_Low_ad<-dat_ssp126_Low_ad[order(dat_ssp126_Low_ad$Elevation),]
  dat_ssp585_Low<-dat_ssp585_Low[order(dat_ssp585_Low$Elevation),]
  dat_ssp585_Low_ad<-dat_ssp585_Low_ad[order(dat_ssp585_Low_ad$Elevation),]
  
  
  
  
  
  ######High elevation raster to point####
  dat_current_High<-as.data.frame(rasterToPoints(current_High))
  colnames(dat_current_High)<-c("Longitude","Latitude","Preabs") #Preabs: presence (1) vs absence (0);
  dat_ssp126_High<-as.data.frame(rasterToPoints(future_ssp126_High))
  colnames(dat_ssp126_High)<-c("Longitude","Latitude","Preabs") #Preabs: presence (1) vs absence (0);
  dat_ssp585_High<-as.data.frame(rasterToPoints(future_ssp585_High))
  colnames(dat_ssp585_High)<-c("Longitude","Latitude","Preabs") #Preabs: presence (1) vs absence (0);
  dat_ssp126_High_ad1<-as.data.frame(rasterToPoints(future_ssp126_High_ad1))
  colnames(dat_ssp126_High_ad1)<-c("Longitude","Latitude","Preabs") #Preabs: presence (1) vs absence (0);
  dat_ssp585_High_ad1<-as.data.frame(rasterToPoints(future_ssp585_High_ad1))
  colnames(dat_ssp585_High_ad1)<-c("Longitude","Latitude","Preabs") #Preabs: presence (1) vs absence (0);
  
  dat_current_High["Elevation"]<-extract(elmat,dat_current_High[,c(1,2)])
  dat_ssp126_High["Elevation"]<-extract(elmat,dat_ssp126_High[,c(1,2)])
  dat_ssp585_High["Elevation"]<-extract(elmat,dat_ssp585_High[,c(1,2)])
  dat_ssp126_High_ad1["Elevation"]<-extract(elmat,dat_ssp126_High_ad1[,c(1,2)])
  dat_ssp585_High_ad1["Elevation"]<-extract(elmat,dat_ssp585_High_ad1[,c(1,2)])
  
  dat_current_High_col<-dat_current_High
  dat_current_High[dat_current_High$Preabs==1&dat_current_High$Elevation<limit,]$Preabs<-0
  dat_ssp126_High_ad<-rbind(dat_ssp126_High[dat_ssp126_High$Elevation>=limit,],dat_ssp126_High_ad1[dat_ssp126_High_ad1$Elevation<limit,])
  dat_ssp585_High_ad<-rbind(dat_ssp585_High[dat_ssp585_High$Elevation>=limit,],dat_ssp585_High_ad1[dat_ssp585_High_ad1$Elevation<limit,])
  
  dat_current_High_col<-dat_current_High_col[order(dat_current_High_col$Elevation),]
  dat_current_High<-dat_current_High[order(dat_current_High$Elevation),]
  dat_ssp126_High<-dat_ssp126_High[order(dat_ssp126_High$Elevation),]
  dat_ssp126_High_ad<-dat_ssp126_High_ad[order(dat_ssp126_High_ad$Elevation),]
  dat_ssp585_High<-dat_ssp585_High[order(dat_ssp585_High$Elevation),]
  dat_ssp585_High_ad<-dat_ssp585_High_ad[order(dat_ssp585_High_ad$Elevation),]


  ##### Low population calculate the suitable habitats change#####
  ##1.Lost habitat (=1 for current, =0 for future)
  Lost_ssp126_Low<-dat_current_Low[dat_current_Low$Preabs==1 & dat_ssp126_Low$Preabs==0,1:2]
  Lost_ssp126a_Low<-dat_current_Low[dat_current_Low$Preabs==1 & dat_ssp126_Low_ad$Preabs==0,1:2]
  Lost_ssp585_Low<-dat_current_Low[dat_current_Low$Preabs==1 & dat_ssp585_Low$Preabs==0,1:2]
  Lost_ssp585a_Low<-dat_current_Low[dat_current_Low$Preabs==1 & dat_ssp585_Low_ad$Preabs==0,1:2]
  
  Lost_ssp126_rate_Low<-nrow(Lost_ssp126_Low)/nrow(dat_current_Low[dat_current_Low$Preabs==1,])*100
  Lost_ssp126a_rate_Low<-nrow(Lost_ssp126a_Low)/nrow(dat_current_Low[dat_current_Low$Preabs==1,])*100
  Lost_ssp585_rate_Low<-nrow(Lost_ssp585_Low)/nrow(dat_current_Low[dat_current_Low$Preabs==1,])*100
  Lost_ssp585a_rate_Low<-nrow(Lost_ssp585a_Low)/nrow(dat_current_Low[dat_current_Low$Preabs==1,])*100
  
  ##2. Colonized habitat (=0 for current, =1 for future)
  Colonized_ssp126_Low<-dat_current_Low_col[dat_current_Low_col$Preabs==0 & dat_ssp126_Low$Preabs==1,1:2]
  Colonized_ssp126a_Low<-dat_current_Low_col[dat_current_Low_col$Preabs==0 & dat_ssp126_Low_ad$Preabs==1,1:2]
  Colonized_ssp585_Low<-dat_current_Low_col[dat_current_Low_col$Preabs==0 & dat_ssp585_Low$Preabs==1,1:2]
  Colonized_ssp585a_Low<-dat_current_Low_col[dat_current_Low_col$Preabs==0 & dat_ssp585_Low_ad$Preabs==1,1:2]
  
  Colonized_ssp126_rate_Low<-nrow(Colonized_ssp126_Low)/nrow(dat_current_Low[dat_current_Low$Preabs==1,])*100
  Colonized_ssp126a_rate_Low<-nrow(Colonized_ssp126a_Low)/nrow(dat_current_Low[dat_current_Low$Preabs==1,])*100
  Colonized_ssp585_rate_Low<-nrow(Colonized_ssp585_Low)/nrow(dat_current_Low[dat_current_Low$Preabs==1,])*100
  Colonized_ssp585a_rate_Low<-nrow(Colonized_ssp585a_Low)/nrow(dat_current_Low[dat_current_Low$Preabs==1,])*100
  
  ##3. Preserved habitat
  Preserved_ssp126_Low<-dat_current_Low[dat_current_Low$Preabs==1 & dat_ssp126_Low$Preabs==1,1:2]
  Preserved_ssp126a_Low<-dat_current_Low[dat_current_Low$Preabs==1 & dat_ssp126_Low_ad$Preabs==1,1:2]
  Preserved_ssp585_Low<-dat_current_Low[dat_current_Low$Preabs==1 & dat_ssp585_Low$Preabs==1,1:2]
  Preserved_ssp585a_Low<-dat_current_Low[dat_current_Low$Preabs==1 & dat_ssp585_Low_ad$Preabs==1,1:2]
  
  Preserved_ssp126_rate_Low<-nrow(Preserved_ssp126_Low)/nrow(dat_current_Low[dat_current_Low$Preabs==1,])*100
  Preserved_ssp126a_rate_Low<-nrow(Preserved_ssp126a_Low)/nrow(dat_current_Low[dat_current_Low$Preabs==1,])*100
  Preserved_ssp585_rate_Low<-nrow(Preserved_ssp585_Low)/nrow(dat_current_Low[dat_current_Low$Preabs==1,])*100
  Preserved_ssp585a_rate_Low<-nrow(Preserved_ssp585a_Low)/nrow(dat_current_Low[dat_current_Low$Preabs==1,])*100

  ##4. Net change in amount of habitat
  net_ssp126_Low<-length(dat_ssp126_Low$Preabs[dat_ssp126_Low$Preabs==1])-length(dat_current_Low_col$Preabs[dat_current_Low_col$Preabs==1])
  print(net_ssp126_Low)
  net_ssp126a_Low<-length(dat_ssp126_Low_ad$Preabs[dat_ssp126_Low_ad$Preabs==1])-length(dat_current_Low_col$Preabs[dat_current_Low_col$Preabs==1])
  print(net_ssp126a_Low)
  net_ssp585_Low<-length(dat_ssp585_Low$Preabs[dat_ssp585_Low$Preabs==1])-length(dat_current_Low_col$Preabs[dat_current_Low_col$Preabs==1])
  print(net_ssp585_Low)
  net_ssp585a_Low<-length(dat_ssp585_Low_ad$Preabs[dat_ssp585_Low_ad$Preabs==1])-length(dat_current_Low_col$Preabs[dat_current_Low_col$Preabs==1])
  print(net_ssp585a_Low)
  
  Net_change_Low<-data.frame(pop="Low",Scenario=c("SSP126","SSP126a","SSP585","SSP585a"),Lost_rate=c(Lost_ssp126_rate_Low,Lost_ssp126a_rate_Low,Lost_ssp585_rate_Low,Lost_ssp585a_rate_Low),Colonized_rate=c(Colonized_ssp126_rate_Low,Colonized_ssp126a_rate_Low,Colonized_ssp585_rate_Low,Colonized_ssp585a_rate_Low),Preserved_rate=c(Preserved_ssp126_rate_Low,Preserved_ssp126a_rate_Low,Preserved_ssp585_rate_Low,Preserved_ssp585a_rate_Low),Netchange=c(net_ssp126_Low,net_ssp126a_Low,net_ssp585_Low,net_ssp585a_Low),Change_rate=c(net_ssp126_Low/nrow(dat_current_Low[dat_current_Low$Preabs==1,]),net_ssp126a_Low/nrow(dat_current_Low[dat_current_Low$Preabs==1,]),net_ssp585_Low/nrow(dat_current_Low[dat_current_Low$Preabs==1,]),net_ssp585a_Low/nrow(dat_current_Low[dat_current_Low$Preabs==1,])))
  write.csv(Net_change_Low,file = paste0("E:/QT/SDM/",limit,"m/SDM_Group/Low_Habitat_Change.csv"))
 
  
  Lost_ssp126_Low["Change"]<-"Lost"
  Lost_ssp126a_Low["Change"]<-"Lost"
  Lost_ssp585_Low["Change"]<-"Lost"
  Lost_ssp585a_Low["Change"]<-"Lost"
  Colonized_ssp126_Low["Change"]<-"Colonized"
  Colonized_ssp126a_Low["Change"]<-"Colonized"
  Colonized_ssp585_Low["Change"]<-"Colonized"
  Colonized_ssp585a_Low["Change"]<-"Colonized"
  
  
  Lost_ssp126_Low["Scenario"]<-"ssp126"
  Lost_ssp126a_Low["Scenario"]<-"ssp126"
  Lost_ssp585_Low["Scenario"]<-"ssp585"
  Lost_ssp585a_Low["Scenario"]<-"ssp585"
  Colonized_ssp126_Low["Scenario"]<-"ssp126"
  Colonized_ssp126a_Low["Scenario"]<-"ssp126"
  Colonized_ssp585_Low["Scenario"]<-"ssp585"
  Colonized_ssp585a_Low["Scenario"]<-"ssp585"
  
  Lost_ssp126_Low["Adaptation"]<-"No"
  Lost_ssp126a_Low["Adaptation"]<-"Yes"
  Lost_ssp585_Low["Adaptation"]<-"No"
  Lost_ssp585a_Low["Adaptation"]<-"Yes"
  Colonized_ssp126_Low["Adaptation"]<-"No"
  Colonized_ssp126a_Low["Adaptation"]<-"Yes"
  Colonized_ssp585_Low["Adaptation"]<-"No"
  Colonized_ssp585a_Low["Adaptation"]<-"Yes"
  
  Change_Low<-rbind(Lost_ssp126_Low,Lost_ssp126a_Low,Lost_ssp585_Low,Lost_ssp585a_Low,Colonized_ssp126_Low,Colonized_ssp126a_Low,Colonized_ssp585_Low,Colonized_ssp585a_Low)
  Change_Low["pop"]<-"Low"
  
  
  
  ##### High population calculate the suitable habitats change#####
  ##1.Lost habitat (=1 for current, =0 for future)
  

  
  Lost_ssp126_High<-dat_current_High[dat_current_High$Preabs==1 & dat_ssp126_High$Preabs==0,1:2]
  Lost_ssp126a_High<-dat_current_High[dat_current_High$Preabs==1 & dat_ssp126_High_ad$Preabs==0,1:2]
  Lost_ssp585_High<-dat_current_High[dat_current_High$Preabs==1 & dat_ssp585_High$Preabs==0,1:2]
  Lost_ssp585a_High<-dat_current_High[dat_current_High$Preabs==1 & dat_ssp585_High_ad$Preabs==0,1:2]

  
  Lost_ssp126_rate_High<-nrow(Lost_ssp126_High)/nrow(dat_current_High[dat_current_High$Preabs==1,])*100
  Lost_ssp126a_rate_High<-nrow(Lost_ssp126a_High)/nrow(dat_current_High[dat_current_High$Preabs==1,])*100
  Lost_ssp585_rate_High<-nrow(Lost_ssp585_High)/nrow(dat_current_High[dat_current_High$Preabs==1,])*100
  Lost_ssp585a_rate_High<-nrow(Lost_ssp585a_High)/nrow(dat_current_High[dat_current_High$Preabs==1,])*100
  
  ##2. Colonized habitat (=0 for current, =1 for future)
  Colonized_ssp126_High<-dat_current_High_col[dat_current_High_col$Preabs==0 & dat_ssp126_High$Preabs==1,1:2]
  Colonized_ssp126a_High<-dat_current_High_col[dat_current_High_col$Preabs==0 & dat_ssp126_High_ad$Preabs==1,1:2]
  Colonized_ssp585_High<-dat_current_High_col[dat_current_High_col$Preabs==0 & dat_ssp585_High$Preabs==1,1:2]
  Colonized_ssp585a_High<-dat_current_High_col[dat_current_High_col$Preabs==0 & dat_ssp585_High_ad$Preabs==1,1:2]
  
  Colonized_ssp126_rate_High<-nrow(Colonized_ssp126_High)/nrow(dat_current_High[dat_current_High$Preabs==1,])*100
  Colonized_ssp126a_rate_High<-nrow(Colonized_ssp126a_High)/nrow(dat_current_High[dat_current_High$Preabs==1,])*100
  Colonized_ssp585_rate_High<-nrow(Colonized_ssp585_High)/nrow(dat_current_High[dat_current_High$Preabs==1,])*100
  Colonized_ssp585a_rate_High<-nrow(Colonized_ssp585a_High)/nrow(dat_current_High[dat_current_High$Preabs==1,])*100
  
  ##3. Preserved habitat
  Preserved_ssp126_High<-dat_current_High[dat_current_High$Preabs==1 & dat_ssp126_High$Preabs==1,1:2]
  Preserved_ssp126a_High<-dat_current_High[dat_current_High$Preabs==1 & dat_ssp126_High_ad$Preabs==1,1:2]
  Preserved_ssp585_High<-dat_current_High[dat_current_High$Preabs==1 & dat_ssp585_High$Preabs==1,1:2]
  Preserved_ssp585a_High<-dat_current_High[dat_current_High$Preabs==1 & dat_ssp585_High_ad$Preabs==1,1:2]
  
  Preserved_ssp126_rate_High<-nrow(Preserved_ssp126_High)/nrow(dat_current_High[dat_current_High$Preabs==1,])*100
  Preserved_ssp126a_rate_High<-nrow(Preserved_ssp126a_High)/nrow(dat_current_High[dat_current_High$Preabs==1,])*100
  Preserved_ssp585_rate_High<-nrow(Preserved_ssp585_High)/nrow(dat_current_High[dat_current_High$Preabs==1,])*100
  Preserved_ssp585a_rate_High<-nrow(Preserved_ssp585a_High)/nrow(dat_current_High[dat_current_High$Preabs==1,])*100
  
  ##4. Net change in amount of habitat
  net_ssp126_High<-length(dat_ssp126_High$Preabs[dat_ssp126_High$Preabs==1])-length(dat_current_High_col$Preabs[dat_current_High_col$Preabs==1])
  print(net_ssp126_High)
  net_ssp126a_High<-length(dat_ssp126_High_ad$Preabs[dat_ssp126_High_ad$Preabs==1])-length(dat_current_High_col$Preabs[dat_current_High_col$Preabs==1])
  print(net_ssp126a_High)
  net_ssp585_High<-length(dat_ssp585_High$Preabs[dat_ssp585_High$Preabs==1])-length(dat_current_High_col$Preabs[dat_current_High_col$Preabs==1])
  print(net_ssp585_High)
  net_ssp585a_High<-length(dat_ssp585_High_ad$Preabs[dat_ssp585_High_ad$Preabs==1])-length(dat_current_High_col$Preabs[dat_current_High_col$Preabs==1])
  print(net_ssp585a_High)
  
  Net_change_High<-data.frame(pop="High",Scenario=c("SSP126","SSP126a","SSP585","SSP585a"),Lost_rate=c(Lost_ssp126_rate_High,Lost_ssp126a_rate_High,Lost_ssp585_rate_High,Lost_ssp585a_rate_High),Colonized_rate=c(Colonized_ssp126_rate_High,Colonized_ssp126a_rate_High,Colonized_ssp585_rate_High,Colonized_ssp585a_rate_High),Preserved_rate=c(Preserved_ssp126_rate_High,Preserved_ssp126a_rate_High,Preserved_ssp585_rate_High,Preserved_ssp585a_rate_High),Netchange=c(net_ssp126_High,net_ssp126a_High,net_ssp585_High,net_ssp585a_High),Change_rate=c(net_ssp126_High/nrow(dat_current_High[dat_current_High$Preabs==1,]),net_ssp126a_High/nrow(dat_current_High[dat_current_High$Preabs==1,]),net_ssp585_High/nrow(dat_current_High[dat_current_High$Preabs==1,]),net_ssp585a_High/nrow(dat_current_High[dat_current_High$Preabs==1,])))
  write.csv(Net_change_High,file =paste0("E:/QT/SDM/",limit,"m/SDM_Group/High_Habitat_Change.csv") )
  
  
  
  Lost_ssp126_High["Change"]<-"Lost"
  Lost_ssp126a_High["Change"]<-"Lost"
  Lost_ssp585_High["Change"]<-"Lost"
  Lost_ssp585a_High["Change"]<-"Lost"
  Colonized_ssp126_High["Change"]<-"Colonized"
  Colonized_ssp126a_High["Change"]<-"Colonized"
  Colonized_ssp585_High["Change"]<-"Colonized"
  Colonized_ssp585a_High["Change"]<-"Colonized"
  
  
  Lost_ssp126_High["Scenario"]<-"ssp126"
  Lost_ssp126a_High["Scenario"]<-"ssp126"
  Lost_ssp585_High["Scenario"]<-"ssp585"
  Lost_ssp585a_High["Scenario"]<-"ssp585"
  Colonized_ssp126_High["Scenario"]<-"ssp126"
  Colonized_ssp126a_High["Scenario"]<-"ssp126"
  Colonized_ssp585_High["Scenario"]<-"ssp585"
  Colonized_ssp585a_High["Scenario"]<-"ssp585"
  
  Lost_ssp126_High["Adaptation"]<-"No"
  Lost_ssp126a_High["Adaptation"]<-"Yes"
  Lost_ssp585_High["Adaptation"]<-"No"
  Lost_ssp585a_High["Adaptation"]<-"Yes"
  Colonized_ssp126_High["Adaptation"]<-"No"
  Colonized_ssp126a_High["Adaptation"]<-"Yes"
  Colonized_ssp585_High["Adaptation"]<-"No"
  Colonized_ssp585a_High["Adaptation"]<-"Yes"
  
  Change_High<-rbind(Lost_ssp126_High,Lost_ssp126a_High,Lost_ssp585_High,Lost_ssp585a_High,Colonized_ssp126_High,Colonized_ssp126a_High,Colonized_ssp585_High,Colonized_ssp585a_High)
  Change_High["pop"]<-"High"
  
  
  Change<-rbind(Change_Low,Change_High)
  Change["Elevation"]<-extract(elmat,Change[,c(1,2)])
  write.csv(Change,file =paste0("E:/QT/SDM/",limit,"m/SDM_Group/Habitat_Change.csv"))
  

  #######Making figures############
  library(RColorBrewer)
  library(ggspatial)
  library(ggThemeAssist)
  
  colormap<-colorRampPalette(brewer.pal(9,"Greys"))(5)
  
  
  
  #########Low population####
  ##Figure 1: Range shift under ssp126
  
  Figure1_Low<-ggplot()+
    geom_raster(data=elmat.df, aes(lon, lat, fill = ele),alpha=0.8)+
    scale_fill_gradientn(values = scales::rescale(c(0,2500,4000,5500,8000)),colors = colormap)+
    geom_point(data=Colonized_ssp126_Low,aes(y=Latitude, x=Longitude,colour="Colonized"),size=0.4,shape = 15)+
    geom_point(data=Lost_ssp126_Low,aes(y=Latitude, x=Longitude,colour="Lost"),size=0.4, shape = 15)+
    geom_point(data=Preserved_ssp126_Low,aes(y=Latitude, x=Longitude,colour="Preserved"),size=0.4, shape = 15)+
    scale_color_manual(values=c("#4169E1","#DC143C","#3CB371"))+
    #guides(colour=guide_legend(title = NULL))+
    labs(col = "Habitat category")+
    guides(col = guide_legend(override.aes = list(size=4)))+
    labs(fill = "Elevation (m)")+
    theme(legend.text = element_text(size = 18,family = "serif",face = "bold"), 
          legend.background = element_rect(fill = NA))+ 
    theme_void()+
    coord_equal()
  #annotation_scale(location = "bl", width_hint = 0.4, 
  #            pad_x = unit(0.8, "in"), pad_y = unit(0.5, "in")) +
  
  #  annotation_north_arrow(location = "tr", which_north = "true", height = unit(1.5,"cm"),width=unit(1.5,"cm"),
  # pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"),
  #style = north_arrow_nautical) 
  
  ggsave(Figure1_Low,file=paste("E:/QT/SDM/",limit,"m/SDM_Group/Fig_ssp126_Low.pdf",sep=""),device = "pdf",height=3,width=4,dpi = 300)
  
  
  ##Figure 2: Range shift under ssp126a
  
  Figure2_Low<-ggplot()+
    geom_raster(data=elmat.df, aes(lon, lat, fill = ele),alpha=0.8)+
    scale_fill_gradientn(values = scales::rescale(c(0,2500,4000,5500,8000)),colors = colormap)+
    geom_point(data=Colonized_ssp126a_Low,aes(y=Latitude, x=Longitude,colour="Colonized"),size=0.4,shape = 15)+
    geom_point(data=Lost_ssp126a_Low,aes(y=Latitude, x=Longitude,colour="Lost"),size=0.4, shape = 15)+
    geom_point(data=Preserved_ssp126a_Low,aes(y=Latitude, x=Longitude,colour="Preserved"),size=0.4, shape = 15)+
    scale_color_manual(values=c("#4169E1","#DC143C","#3CB371"))+
    #guides(colour=guide_legend(title = NULL))+
    labs(col = "Habitat category")+
    guides(col = guide_legend(override.aes = list(size=4)))+
    labs(fill = "Elevation (m)")+
    theme(legend.text = element_text(size = 18,family = "serif",face = "bold"), 
          legend.background = element_rect(fill = NA))+ 
    theme_void()+
    coord_equal()
  #annotation_scale(location = "bl", width_hint = 0.4, 
  #            pad_x = unit(0.8, "in"), pad_y = unit(0.5, "in")) +
  
  #  annotation_north_arrow(location = "tr", which_north = "true", height = unit(1.5,"cm"),width=unit(1.5,"cm"),
  # pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"),
  #style = north_arrow_nautical) 
  
  ggsave(Figure2_Low,file=paste("E:/QT/SDM/",limit,"m/SDM_Group/Fig_ssp126a_Low.pdf",sep=""),device = "pdf",height=3,width=4,dpi = 300)
  
  
  ##Figure 3: Range shift under ssp585
  
  Figure3_Low<-ggplot()+
    geom_raster(data=elmat.df, aes(lon, lat, fill = ele),alpha=0.8)+
    scale_fill_gradientn(values = scales::rescale(c(0,2500,4000,5500,8000)),colors = colormap)+
    geom_point(data=Colonized_ssp585_Low,aes(y=Latitude, x=Longitude,colour="Colonized"),size=0.4,shape = 15)+
    geom_point(data=Lost_ssp585_Low,aes(y=Latitude, x=Longitude,colour="Lost"),size=0.4, shape = 15)+
    geom_point(data=Preserved_ssp585_Low,aes(y=Latitude, x=Longitude,colour="Preserved"),size=0.4, shape = 15)+
    scale_color_manual(values=c("#4169E1","#DC143C","#3CB371"))+
    #guides(colour=guide_legend(title = NULL))+
    labs(col = "Habitat category")+
    guides(col = guide_legend(override.aes = list(size=4)))+
    labs(fill = "Elevation (m)")+
    theme(legend.text = element_text(size = 18,family = "serif",face = "bold"), 
          legend.background = element_rect(fill = NA))+ 
    theme_void()+
    coord_equal()
  #annotation_scale(location = "bl", width_hint = 0.4, 
  #            pad_x = unit(0.8, "in"), pad_y = unit(0.5, "in")) +
  
  #  annotation_north_arrow(location = "tr", which_north = "true", height = unit(1.5,"cm"),width=unit(1.5,"cm"),
  # pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"),
  #style = north_arrow_nautical) 
  
  ggsave(Figure3_Low,file=paste("E:/QT/SDM/",limit,"m/SDM_Group/Fig_ssp585_Low.pdf",sep=""),device = "pdf",height=3,width=4,dpi = 300)
  
  
  
  ##Figure 4: Range shift under ssp585a
  
  Figure4_Low<-ggplot()+
    geom_raster(data=elmat.df, aes(lon, lat, fill = ele),alpha=0.8)+
    scale_fill_gradientn(values = scales::rescale(c(0,2500,4000,5500,8000)),colors = colormap)+
    geom_point(data=Colonized_ssp585a_Low,aes(y=Latitude, x=Longitude,colour="Colonized"),size=0.4,shape = 15)+
    geom_point(data=Lost_ssp585a_Low,aes(y=Latitude, x=Longitude,colour="Lost"),size=0.4, shape = 15)+
    geom_point(data=Preserved_ssp585a_Low,aes(y=Latitude, x=Longitude,colour="Preserved"),size=0.4, shape = 15)+
    scale_color_manual(values=c("#4169E1","#DC143C","#3CB371"))+
    #guides(colour=guide_legend(title = NULL))+
    labs(col = "Habitat category")+
    guides(col = guide_legend(override.aes = list(size=4)))+
    labs(fill = "Elevation (m)")+
    theme(legend.text = element_text(size = 18,family = "serif",face = "bold"), 
          legend.background = element_rect(fill = NA))+ 
    theme_void()+
    coord_equal()
  #annotation_scale(location = "bl", width_hint = 0.4, 
  #            pad_x = unit(0.8, "in"), pad_y = unit(0.5, "in")) +
  
  #  annotation_north_arrow(location = "tr", which_north = "true", height = unit(1.5,"cm"),width=unit(1.5,"cm"),
  # pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"),
  #style = north_arrow_nautical) 
  
  ggsave(Figure4_Low,file=paste("E:/QT/SDM/",limit,"m/SDM_Group/Fig_ssp585a_Low.pdf",sep=""),device = "pdf",height=3,width=4,dpi = 300)
  
  
  ##Figure : Amount of habitat change
  
  fig<-ggplot(Net_change_Low,aes(x=Scenario,y=Change_rate*100,fill=Scenario))+
    geom_bar(stat="identity",width=0.6)+
    scale_fill_manual(values = c("#FFD700","#4169E1","#DC143C","#9370DB"))+
    # ylim(-2,2)+
    theme_light()+
    #scale_x_discrete(labels=c("lu"="LUC","wc"="GCC", "wclu"="BOTH"))+
    ylab("Net habitat increase (%)")+
    geom_hline(yintercept=0,linetype="dashed",size=0.8)+
    theme(strip.text = element_text(face = "italic"),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+ 
    theme(axis.title = element_text(size = 18), 
          axis.text = element_text(size = 18), 
          legend.position = "none")+
    theme(axis.ticks = element_line(colour = "black"), 
          axis.text = element_text(colour = "black"))
  
  ggsave(fig,file=paste("E:/QT/SDM/",limit,"m/SDM_Group/Fig_Netchange_Low.pdf",sep=""),device = "pdf",height=9,width=10,dpi = 300)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #########High population####
  ##Figure 1: Range shift under ssp126
  
  Figure1_High<-ggplot()+
    geom_raster(data=elmat.df, aes(lon, lat, fill = ele),alpha=0.8)+
    scale_fill_gradientn(values = scales::rescale(c(0,2500,4000,5500,8000)),colors = colormap)+
    geom_point(data=Colonized_ssp126_High,aes(y=Latitude, x=Longitude,colour="Colonized"),size=0.4,shape = 15)+
    geom_point(data=Lost_ssp126_High,aes(y=Latitude, x=Longitude,colour="Lost"),size=0.4, shape = 15)+
    geom_point(data=Preserved_ssp126_High,aes(y=Latitude, x=Longitude,colour="Preserved"),size=0.4, shape = 15)+
    scale_color_manual(values=c("#4169E1","#DC143C","#3CB371"))+
    #guides(colour=guide_legend(title = NULL))+
    labs(col = "Habitat category")+
    guides(col = guide_legend(override.aes = list(size=4)))+
    labs(fill = "Elevation (m)")+
    theme(legend.text = element_text(size = 18,family = "serif",face = "bold"), 
          legend.background = element_rect(fill = NA))+ 
    theme_void()+
    coord_equal()
  #annotation_scale(location = "bl", width_hint = 0.4, 
  #            pad_x = unit(0.8, "in"), pad_y = unit(0.5, "in")) +
  
  #  annotation_north_arrow(location = "tr", which_north = "true", height = unit(1.5,"cm"),width=unit(1.5,"cm"),
  # pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"),
  #style = north_arrow_nautical) 
  
  ggsave(Figure1_High,file=paste("E:/QT/SDM/",limit,"m/SDM_Group/Fig_ssp126_High.pdf",sep=""),device = "pdf",height=3,width=4,dpi = 300)
  
  
  ##Figure 2: Range shift under ssp126a
  
  Figure2_High<-ggplot()+
    geom_raster(data=elmat.df, aes(lon, lat, fill = ele),alpha=0.8)+
    scale_fill_gradientn(values = scales::rescale(c(0,2500,4000,5500,8000)),colors = colormap)+
    geom_point(data=Colonized_ssp126a_High,aes(y=Latitude, x=Longitude,colour="Colonized"),size=0.4,shape = 15)+
    geom_point(data=Lost_ssp126a_High,aes(y=Latitude, x=Longitude,colour="Lost"),size=0.4, shape = 15)+
    geom_point(data=Preserved_ssp126a_High,aes(y=Latitude, x=Longitude,colour="Preserved"),size=0.4, shape = 15)+
    scale_color_manual(values=c("#4169E1","#DC143C","#3CB371"))+
    #guides(colour=guide_legend(title = NULL))+
    labs(col = "Habitat category")+
    guides(col = guide_legend(override.aes = list(size=4)))+
    labs(fill = "Elevation (m)")+
    theme(legend.text = element_text(size = 18,family = "serif",face = "bold"), 
          legend.background = element_rect(fill = NA))+ 
    theme_void()+
    coord_equal()
  #annotation_scale(location = "bl", width_hint = 0.4, 
  #            pad_x = unit(0.8, "in"), pad_y = unit(0.5, "in")) +
  
  #  annotation_north_arrow(location = "tr", which_north = "true", height = unit(1.5,"cm"),width=unit(1.5,"cm"),
  # pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"),
  #style = north_arrow_nautical) 
  
  ggsave(Figure2_High,file=paste("E:/QT/SDM/",limit,"m/SDM_Group/Fig_ssp126a_High.pdf",sep=""),device = "pdf",height=3,width=4,dpi = 300)
  
  
  ##Figure 3: Range shift under ssp585
  
  Figure3_High<-ggplot()+
    geom_raster(data=elmat.df, aes(lon, lat, fill = ele),alpha=0.8)+
    scale_fill_gradientn(values = scales::rescale(c(0,2500,4000,5500,8000)),colors = colormap)+
    geom_point(data=Colonized_ssp585_High,aes(y=Latitude, x=Longitude,colour="Colonized"),size=0.4,shape = 15)+
    geom_point(data=Lost_ssp585_High,aes(y=Latitude, x=Longitude,colour="Lost"),size=0.4, shape = 15)+
    geom_point(data=Preserved_ssp585_High,aes(y=Latitude, x=Longitude,colour="Preserved"),size=0.4, shape = 15)+
    scale_color_manual(values=c("#4169E1","#DC143C","#3CB371"))+
    #guides(colour=guide_legend(title = NULL))+
    labs(col = "Habitat category")+
    guides(col = guide_legend(override.aes = list(size=4)))+
    labs(fill = "Elevation (m)")+
    theme(legend.text = element_text(size = 18,family = "serif",face = "bold"), 
          legend.background = element_rect(fill = NA))+ 
    theme_void()+
    coord_equal()
  #annotation_scale(location = "bl", width_hint = 0.4, 
  #            pad_x = unit(0.8, "in"), pad_y = unit(0.5, "in")) +
  
  #  annotation_north_arrow(location = "tr", which_north = "true", height = unit(1.5,"cm"),width=unit(1.5,"cm"),
  # pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"),
  #style = north_arrow_nautical) 
  
  ggsave(Figure3_High,file=paste("E:/QT/SDM/",limit,"m/SDM_Group/Fig_ssp585_High.pdf",sep=""),device = "pdf",height=3,width=4,dpi = 300)
  
  
  
  ##Figure 4: Range shift under ssp585a
  
  Figure4_High<-ggplot()+
    geom_raster(data=elmat.df, aes(lon, lat, fill = ele),alpha=0.8)+
    scale_fill_gradientn(values = scales::rescale(c(0,2500,4000,5500,8000)),colors = colormap)+
    geom_point(data=Colonized_ssp585a_High,aes(y=Latitude, x=Longitude,colour="Colonized"),size=0.4,shape = 15)+
    geom_point(data=Lost_ssp585a_High,aes(y=Latitude, x=Longitude,colour="Lost"),size=0.4, shape = 15)+
    geom_point(data=Preserved_ssp585a_High,aes(y=Latitude, x=Longitude,colour="Preserved"),size=0.4, shape = 15)+
    scale_color_manual(values=c("#4169E1","#DC143C","#3CB371"))+
    #guides(colour=guide_legend(title = NULL))+
    labs(col = "Habitat category")+
    guides(col = guide_legend(override.aes = list(size=4)))+
    labs(fill = "Elevation (m)")+
    theme(legend.text = element_text(size = 18,family = "serif",face = "bold"), 
          legend.background = element_rect(fill = NA))+ 
    theme_void()+
    coord_equal()
  #annotation_scale(location = "bl", width_hint = 0.4, 
  #            pad_x = unit(0.8, "in"), pad_y = unit(0.5, "in")) +
  
  #  annotation_north_arrow(location = "tr", which_north = "true", height = unit(1.5,"cm"),width=unit(1.5,"cm"),
  # pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"),
  #style = north_arrow_nautical) 
  
  ggsave(Figure4_High,file=paste("E:/QT/SDM/",limit,"m/SDM_Group/Fig_ssp585a_High.pdf",sep=""),device = "pdf",height=3,width=4,dpi = 300)
  
  
  ##Figure : Amount of habitat change
  
  fig<-ggplot(Net_change_High,aes(x=Scenario,y=Change_rate*100,fill=Scenario))+
    geom_bar(stat="identity",width=0.6)+
    scale_fill_manual(values = c("#FFD700","#4169E1","#DC143C","#9370DB"))+
    # ylim(-2,2)+
    theme_light()+
    #scale_x_discrete(labels=c("lu"="LUC","wc"="GCC", "wclu"="BOTH"))+
    ylab("Net habitat increase (%)")+
    geom_hline(yintercept=0,linetype="dashed",size=0.8)+
    theme(strip.text = element_text(face = "italic"),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+ 
    theme(axis.title = element_text(size = 18), 
          axis.text = element_text(size = 18), 
          legend.position = "none")+
    theme(axis.ticks = element_line(colour = "black"), 
          axis.text = element_text(colour = "black"))
  
  ggsave(fig,file=paste("E:/QT/SDM/",limit,"m/SDM_Group/Fig_Netchange_High.pdf",sep=""),device = "pdf",height=9,width=10,dpi = 300)
  
  
  



  #limit=3600
  
  ################Habitat change#############
  hab_change<-read.csv(paste0("E:/QT/SDM/",limit,"m/SDM_Group/Habitat_Change.csv"))

  hab_change_Clim<-read.csv(paste0("E:/QT/SDM/",limit,"m/SDM_Clim/Habitat_Change.csv"))
  hab_change_Clim<-hab_change_Clim[hab_change_Clim$Adaptation=="No",]
  hab_change_Clim$Adaptation<-"Only Climate"
  
  change_data<-data.frame()
  
  for (pop in c("Low","High")) {
   # pop="High"
    
    pop_change<-hab_change[hab_change$pop==pop,]
    
    for (ssp in c("ssp126","ssp585")) {
      #ssp="ssp585"
      
      pop_change_ssp<-pop_change[pop_change$Scenario==ssp,
                                 ]
      number_lost<- nrow(pop_change_ssp[pop_change_ssp$Change=="Lost",])
      
      for (adap in c("No","Yes")) {
       # adap="Yes"
        
      number_colonized_high<-nrow(pop_change_ssp[pop_change_ssp$Change=="Colonized"&pop_change_ssp$Adaptation==adap&pop_change_ssp$Elevation>=limit,])
      number_colonized_low<-nrow(pop_change_ssp[pop_change_ssp$Change=="Colonized"&pop_change_ssp$Adaptation==adap&pop_change_ssp$Elevation<limit,])

      change_number<-data.frame(pop=pop,ssp=ssp,adaptation=adap,change=c("Lost","Colonized_high","Colonized_low"),change_number=c(number_lost/2,number_colonized_high,number_colonized_low))
      change_data<-rbind(change_data,change_number)
      }
    }
  }
  
  
  low_change<-hab_change[hab_change$pop=="Low",]
  low_change_ssp126<-low_change[low_change$Scenario=="ssp126",]
  low_change_ssp585<-low_change[low_change$Scenario=="ssp585",]
  
  low_change_Clim<-hab_change_Clim[hab_change_Clim$pop=="Low",]
  low_change_ssp126_Clim<-low_change_Clim[low_change_Clim$Scenario=="ssp126",]
  low_change_ssp585_Clim<-low_change_Clim[low_change_Clim$Scenario=="ssp585",]
  
  nrow(low_change_ssp126_Clim[low_change_ssp126_Clim$Change=="Lost",])
  nrow(low_change_ssp585_Clim[low_change_ssp585_Clim$Change=="Lost",])
  
  nrow(low_change_ssp126_Clim[low_change_ssp126_Clim$Change=="Colonized",])
  nrow(low_change_ssp126_Clim[low_change_ssp126_Clim$Change=="Colonized"&low_change_ssp126_Clim$Elevation>=limit,])
  
  nrow(low_change_ssp585_Clim[low_change_ssp585_Clim$Change=="Colonized",])
  nrow(low_change_ssp585_Clim[low_change_ssp585_Clim$Change=="Colonized"&low_change_ssp585_Clim$Elevation>=limit,])
  
  nrow(low_change_ssp126[low_change_ssp126$Change=="Lost",])
  nrow(low_change_ssp585[low_change_ssp585$Change=="Lost",])
  
  nrow(low_change_ssp126[low_change_ssp126$Change=="Colonized"&low_change_ssp126$Adaptation=="No",])
  nrow(low_change_ssp126[low_change_ssp126$Change=="Colonized"&low_change_ssp126$Adaptation=="No"&low_change_ssp126$Elevation>=limit,])
  
  nrow(low_change_ssp585[low_change_ssp585$Change=="Colonized"&low_change_ssp585$Adaptation=="No",])
  nrow(low_change_ssp585[low_change_ssp585$Change=="Colonized"&low_change_ssp585$Adaptation=="No"&low_change_ssp585$Elevation>=limit,])

  
  high_change<-hab_change[hab_change$pop=="High",]
  high_change_ssp126<-high_change[high_change$Scenario=="ssp126",]
  high_change_ssp585<-high_change[high_change$Scenario=="ssp585",]
  
  
  high_change_Clim<-hab_change_Clim[hab_change_Clim$pop=="High",]
  high_change_ssp126_Clim<-high_change_Clim[high_change_Clim$Scenario=="ssp126",]
  high_change_ssp585_Clim<-high_change_Clim[high_change_Clim$Scenario=="ssp585",]
  
  nrow(high_change_ssp126_Clim[high_change_ssp126_Clim$Change=="Lost",])
  nrow(high_change_ssp585_Clim[high_change_ssp585_Clim$Change=="Lost",])
  
  nrow(high_change_ssp126_Clim[high_change_ssp126_Clim$Change=="Colonized",])
  nrow(high_change_ssp126_Clim[high_change_ssp126_Clim$Change=="Colonized"&high_change_ssp126_Clim$Elevation>=limit,])
  
  nrow(high_change_ssp585_Clim[high_change_ssp585_Clim$Change=="Colonized",])
  nrow(high_change_ssp585_Clim[high_change_ssp585_Clim$Change=="Colonized"&high_change_ssp585_Clim$Elevation>=limit,])
  
  nrow(high_change_ssp126[high_change_ssp126$Change=="Lost",])
  nrow(high_change_ssp585[high_change_ssp585$Change=="Lost",])
  
  nrow(high_change_ssp126[high_change_ssp126$Change=="Colonized"&high_change_ssp126$Adaptation=="No",])
  nrow(high_change_ssp126[high_change_ssp126$Change=="Colonized"&high_change_ssp126$Adaptation=="No"&high_change_ssp126$Elevation>=limit,])
  
  nrow(high_change_ssp585[high_change_ssp585$Change=="Colonized"&high_change_ssp585$Adaptation=="No",])
  nrow(high_change_ssp585[high_change_ssp585$Change=="Colonized"&high_change_ssp585$Adaptation=="No"&high_change_ssp585$Elevation>=limit,])
  
 
 #hab_change_all<-rbind(hab_change,hab_change_Clim)
  
# hab_change_all$Adaptation<-factor(hab_change_all$Adaptation,levels = c("Only Climate","No","Yes"))
  
   for (pop in c("Low","High")) {
    
    pop_change<-hab_change[hab_change$pop==pop,]
    
      fig_change<- ggplot() + 
        geom_vline(xintercept = limit,colour="#990000",size=1)+
        
        geom_histogram(data = subset(pop_change, Change== "Colonized"), 
                       bins = 50, 
                       aes(x = Elevation, y = ..count..), 
                       fill = "#4169E1",
                       color="white")+
        
        
        geom_histogram(data = subset(pop_change, Change== "Lost"), 
                       bins = 50, aes(x = Elevation, 
                                      y = -..count..),
                       fill = "#DC143C",
                       color="white")+
      
      facet_grid(Adaptation~Scenario )+
        scale_y_continuous(breaks = seq(-50,100,25))+
        xlim(1000,5000)+
        scale_x_continuous(limits=c(1000,5000),breaks = seq(1000,5000,1000))+
        ylab("Habitat Changes (n)")+
        theme_minimal()+
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
      
      ggsave(fig_change,file=paste("E:/QT/SDM/",limit,"m/SDM_Group/Fig_Habitatchanges_",pop,".pdf",sep=""),device = "pdf",height=4,width=12,dpi = 300)
      
    
  }
  
  
  


  
  
  }
  
  
