#######Point to Raster/ACT/Warming Tolerance/Waterloss/MR/ENB########
####################Zhong-wen Jiang 2022/01/09#######################
library(raster)

Traits<-c("ACT","WT","CT","ENB","OC","QMET","WL","RWL")

lizard_traits<-read.csv("E:/QT/NicheMapR/Data/Lizard_traits.csv")

for (time in c("1970_2000","ssp126_2081_2100","ssp585_2081_2100")) {
  
  for (Ele in c("3400m","4200m")) {
    
  #  time="1970_2000"
  #  Ele="4200m"
    
    CTmax<-lizard_traits[lizard_traits$﻿Ele==Ele,]$CT_max
    CTmin<-lizard_traits[lizard_traits$﻿Ele==Ele,]$CT_min
    
    for (i in 4:11) {
      
   
    if(i==5) {
      
      Trait<-Traits[i-3]
      Ras_point<-read.csv(paste0("E:/QT/NicheMapR/Result/Eco_result_",time,"_",Ele,".csv"))[,c(1,2,i)]
      Ras_point["WT"]<-CTmax-Ras_point[,3]
      Ras_point<-Ras_point[,c(1,2,4)]
      Raster_traits<-rasterFromXYZ(Ras_point)   
      writeRaster(Raster_traits,paste0("E:/QT/NicheMapR/Result/Eco_result_",time,"_",Ele,"_",Trait,".tif"),format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
      
        
    } else if(i==6){
      
        Trait<-Traits[i-3]
        Ras_point<-read.csv(paste0("E:/QT/NicheMapR/Result/Eco_result_",time,"_",Ele,".csv"))[,c(1,2,i)]
        Ras_point["CT"]<-Ras_point[,3]-CTmin
        Ras_point<-Ras_point[,c(1,2,4)]
        Raster_traits<-rasterFromXYZ(Ras_point)   
        writeRaster(Raster_traits,paste0("E:/QT/NicheMapR/Result/Eco_result_",time,"_",Ele,"_",Trait,".tif"),format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
        
    }else
      {
      
       Trait<-Traits[i-3]
     Ras_point<-read.csv(paste0("E:/QT/NicheMapR/Result/Eco_result_",time,"_",Ele,".csv"))[,c(1,2,i)]
     Raster_traits<-rasterFromXYZ(Ras_point)   
     writeRaster(Raster_traits,paste0("E:/QT/NicheMapR/Result/Eco_result_",time,"_",Ele,"_",Trait,".tif"),format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
     
    }
      
  print(paste0(time,"_",Ele,"_",Trait))
    }
    
  }
}






#####When in Romae, do as the Romans do.#######
##### Elevation: 2000-4500m IUCN/ 2232-4576 our collected point
#### Here, the lizard group divided into 2600m&3400m and 4200m.
#### Maybe  Ele< 3x00m is 2600m and 3400m, Ele>3x00m is 4200m?



Coordinate_ras<-read.csv("E:/QT/NicheMapR/Result/Eco_result_ssp585_2081_2100_4200m.csv")[,1:3]


ele_map<-rasterFromXYZ(Coordinate_ras) 


for (time in c("1970_2000","ssp126_2081_2100","ssp585_2081_2100")) {
  
 # time="1970_2000"
  
  for (i in 4:11) {
    
   # i=4
    
    Trait<-Traits[i-3]
    
    ras_map<-data.frame()
    
    if(i==5) {
    
    for (p in 1:nrow(Coordinate_ras)) {
      lon<-Coordinate_ras[p,]$longitude
      lat<-Coordinate_ras[p,]$latitude
      
      if(Coordinate_ras[p,]$ele<3500){
        
        CTmax<-lizard_traits[lizard_traits$﻿Ele=="3400m",]$CT_max
        CTmin<-lizard_traits[lizard_traits$﻿Ele=="3400m",]$CT_min
        
        c<-read.csv(paste0("E:/QT/NicheMapR/Result/Eco_result_",time,"_3400m.csv"))
        a<-subset(c,longitude==lon&latitude==lat)[,c(1,2,i)]
        a["WT"]<-CTmax-a[,3]
        a<-a[,c(1,2,4)]
        
        ras_map<-rbind(ras_map,a)
      }else{
        CTmax<-lizard_traits[lizard_traits$﻿Ele=="4200m",]$CT_max
        CTmin<-lizard_traits[lizard_traits$﻿Ele=="4200m",]$CT_min
        
        c<-read.csv(paste0("E:/QT/NicheMapR/Result/Eco_result_",time,"_4200m.csv"))
        a<-subset(c,longitude==lon&latitude==lat)[,c(1,2,i)]
        a["WT"]<-CTmax-a[,3]
        a<-a[,c(1,2,4)]
        ras_map<-rbind(ras_map,a)
      }
    }
    
    } else if(i==6){
      for (p in 1:nrow(Coordinate_ras)) {
        lon<-Coordinate_ras[p,]$longitude
        lat<-Coordinate_ras[p,]$latitude
        if(Coordinate_ras[p,]$ele<3500){
          
          CTmax<-lizard_traits[lizard_traits$﻿Ele=="3400m",]$CT_max
          CTmin<-lizard_traits[lizard_traits$﻿Ele=="3400m",]$CT_min
          
          c<-read.csv(paste0("E:/QT/NicheMapR/Result/Eco_result_",time,"_3400m.csv"))
          a<-subset(c,longitude==lon&latitude==lat)[,c(1,2,i)]
          a["CT"]<-a[,3]-CTmin
          a<-a[,c(1,2,4)]
          
          ras_map<-rbind(ras_map,a)
        }else{
          CTmax<-lizard_traits[lizard_traits$﻿Ele=="4200m",]$CT_max
          CTmin<-lizard_traits[lizard_traits$﻿Ele=="4200m",]$CT_min
          
          c<-read.csv(paste0("E:/QT/NicheMapR/Result/Eco_result_",time,"_4200m.csv"))
          a<-subset(c,longitude==lon&latitude==lat)[,c(1,2,i)]
          
          a["CT"]<-a[,3]-CTmin
          a<-a[,c(1,2,4)]
          ras_map<-rbind(ras_map,a)
        }
      }
    } else {
      for (p in 1:nrow(Coordinate_ras)) {
        
        lon<-Coordinate_ras[p,]$longitude
        lat<-Coordinate_ras[p,]$latitude
        
        if(Coordinate_ras[p,]$ele<3500){
          
          c<-read.csv(paste0("E:/QT/NicheMapR/Result/Eco_result_",time,"_3400m.csv"))
          a<-subset(c,longitude==lon&latitude==lat)[,c(1,2,i)]
         
          ras_map<-rbind(ras_map,a)
        }else{
          c<-read.csv(paste0("E:/QT/NicheMapR/Result/Eco_result_",time,"_4200m.csv"))
          a<-subset(c,longitude==lon&latitude==lat)[,c(1,2,i)]
          
           ras_map<-rbind(ras_map,a)
        }
      }
      
    }
    
    Raster_map<-rasterFromXYZ(ras_map)   
    
    writeRaster(Raster_map,paste0("E:/QT/NicheMapR/Result/Eco_result_",time,"_",Trait,".tif"),format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
    print(paste0(time,"_",Trait))
  }
}
