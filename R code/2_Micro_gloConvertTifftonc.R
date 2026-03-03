######################Convert Stack tiff to nc file##############
######################Zhong wen Jiang 2021 helped by MCR
#############
rm(list = ls())
library(raster)
library(ncdf4)

#ref<-raster("/Volumes/My Passport/QT/SDM/CMIP6/Current_10m/wc2.1_10m_bio_1.tif")
global_climate_tif<-stack("/Volumes/My Passport/QT/NicheMapR/Data/global_climate_JZW/global_climate.nc")
#Current<-mask(global_climate_tif,ref)
#Current<-stack(Current)
Current<-global_climate_tif

Altitude <- Current[[1]]      #####Rain
#RAINFALL <- Current[[2:13]]      #####Rain
RAINYDAYS <- Current[[14:25]]
WNMAXX <-  Current[[26:37]]    ###Wind
WNMINN <- WNMAXX * 0.1
#TMINN <-  Current[[38:49]]     ####Tempera
#TMAXX <-  Current[[50:61]]    
#ALLMINTEMPS <- TMINN
#ALLMAXTEMPS <- TMAXX
#ALLTEMPS <- cbind(ALLMAXTEMPS, ALLMINTEMPS) 
RHMINN <-  Current[[62:73]]    ####Relative hum
RHMAXX <-  Current[[74:85]]
CCMINN <-  Current[[86:97]]    ###cloud cover?




###############Current1970_2000
path_current_pre<-"/Volumes/My Passport/QT/NicheMapR/Data/Mon_Precipitation/1970-2000/"
current_prefile<-list.files(path_current_pre,pattern="*.tif")
current_pre<-lapply(current_prefile,function(x) paste(path_current_pre,x,sep=""))
current_pre<-lapply(current_pre,raster)
current_pre<-do.call(stack,current_pre)

path_current_Tmin<-"/Volumes/My Passport/QT/NicheMapR/Data/Mon_Tmin/1970-2000/"
current_Tminfile<-list.files(path_current_Tmin,pattern="*.tif")
current_Tmin<-lapply(current_Tminfile,function(x) paste(path_current_Tmin,x,sep=""))
current_Tmin<-lapply(current_Tmin,raster)
current_Tmin<-do.call(stack,current_Tmin)
current_Tmin@layers<-lapply(current_Tmin@layers,function(x)x*10)


path_current_Tmax<-"/Volumes/My Passport/QT/NicheMapR/Data/Mon_Tmax/1970-2000/"
current_Tmaxfile<-list.files(path_current_Tmax,pattern="*.tif")
current_Tmax<-lapply(current_Tmaxfile,function(x) paste(path_current_Tmax,x,sep=""))
current_Tmax<-lapply(current_Tmax,raster)
current_Tmax<-do.call(stack,current_Tmax)
current_Tmax@layers<-lapply(current_Tmax@layers,function(x)x*10)

RAINFALL <- current_pre
TMINN <-  current_Tmin     ####Tempera
TMAXX <-  current_Tmax   
ALLMINTEMPS <- TMINN
ALLMAXTEMPS <- TMAXX
ALLTEMPS <- cbind(ALLMAXTEMPS, ALLMINTEMPS) 

Current1970_2000<-stack(Altitude,RAINFALL,RAINYDAYS,WNMAXX,TMINN,TMAXX,RHMINN,RHMAXX,CCMINN)

setwd("/Volumes/My Passport/QT/NicheMapR/Data/global_climate_JZW")
writeRaster(Current1970_2000, "current1970_2000.nc", overwrite=TRUE, format="CDF",     varname="Temperature", varunit="degC", 
            longname="Temperature -- raster stack to netCDF", xname="Longitude",   yname="Latitude", zname="Time (Month)")








#############################2081-2100 ssp126
path_tn126<-"/Volumes/My Passport/QT/NicheMapR/Data/Mon_Tmin/2081-2100/ssp126/"
tn126<-stack()
for (i in 1:12) {
  stack_mon<-stack()
  for (GCMs in c("BCC-CSM2-MR","CNRM-CM6-1","CNRM-ESM2-1","CanESM5","IPSL-CM6A-LR","MIROC-ES2L","MIROC6","MRI-ESM2-0")) {
    tn126file<-raster(paste0(path_tn126,"wc2.1_10m_tmin_",GCMs,"_ssp126_2081-2100.tif"),band=i)
    stack_mon<-stack(stack_mon,tn126file)
    
  }
  stack_mean<-mean(stack_mon)
  tn126<-stack(tn126,stack_mean)

}


path_tx126<-"/Volumes/My Passport/QT/NicheMapR/Data/Mon_Tmax/2081-2100/ssp126/"
tx126<-stack()
for (i in 1:12) {
  stack_mon<-stack()
  for (GCMs in c("BCC-CSM2-MR","CNRM-CM6-1","CNRM-ESM2-1","CanESM5","IPSL-CM6A-LR","MIROC-ES2L","MIROC6","MRI-ESM2-0")) {
    tx126file<-raster(paste0(path_tx126,"wc2.1_10m_tmax_",GCMs,"_ssp126_2081-2100.tif"),band=i)
    stack_mon<-stack(stack_mon,tx126file)
    
  }
  stack_mean<-mean(stack_mon)
  tx126<-stack(tx126,stack_mean)
}

path_pr126<-"/Volumes/My Passport/QT/NicheMapR/Data/Mon_Precipitation/2081-2100/ssp126/"
pr126<-stack()
for (i in 1:12) {
  stack_mon<-stack()
  for (GCMs in c("BCC-CSM2-MR","CNRM-CM6-1","CNRM-ESM2-1","CanESM5","IPSL-CM6A-LR","MIROC-ES2L","MIROC6","MRI-ESM2-0")) {
    pr126file<-raster(paste0(path_pr126,"wc2.1_10m_prec_",GCMs,"_ssp126_2081-2100.tif"),band=i)
    stack_mon<-stack(stack_mon,pr126file)
    
  }
  stack_mean<-mean(stack_mon)
  pr126<-stack(pr126,stack_mean)
}

tn126@layers<-lapply(tn126@layers,function(x)x*10)
tx126@layers<-lapply(tx126@layers,function(x)x*10)



RAINFALL <- pr126 
TMINN <-  tn126     ####Tempera
TMAXX <-  tx126
ALLMINTEMPS <- TMINN
ALLMAXTEMPS <- TMAXX
ALLTEMPS <- cbind(ALLMAXTEMPS, ALLMINTEMPS) 

ssp126_2081_2100<-stack(Altitude,RAINFALL,RAINYDAYS,WNMAXX,TMINN,TMAXX,RHMINN,RHMAXX,CCMINN)

setwd("/Volumes/My Passport/QT/NicheMapR/Data/global_climate_JZW")
writeRaster(ssp126_2081_2100, "ssp126_2081_2100.nc", overwrite=TRUE, format="CDF",     varname="Temperature", varunit="degC", 
            longname="Temperature -- raster stack to netCDF", xname="Longitude",   yname="Latitude", zname="Time (Month)")




#############################2081-2100 ssp585
path_tn585<-"/Volumes/My Passport/QT/NicheMapR/Data/Mon_Tmin/2081-2100/ssp585/"
tn585<-stack()
for (i in 1:12) {
  stack_mon<-stack()
  for (GCMs in c("BCC-CSM2-MR","CNRM-CM6-1","CNRM-ESM2-1","CanESM5","IPSL-CM6A-LR","MIROC-ES2L","MIROC6","MRI-ESM2-0")) {
    tn585file<-raster(paste0(path_tn585,"wc2.1_10m_tmin_",GCMs,"_ssp585_2081-2100.tif"),band=i)
    stack_mon<-stack(stack_mon,tn585file)
    
  }
  stack_mean<-mean(stack_mon)
  tn585<-stack(tn585,stack_mean)
}


path_tx585<-"/Volumes/My Passport/QT/NicheMapR/Data/Mon_Tmax/2081-2100/ssp585/"
tx585<-stack()
for (i in 1:12) {
  stack_mon<-stack()
  for (GCMs in c("BCC-CSM2-MR","CNRM-CM6-1","CNRM-ESM2-1","CanESM5","IPSL-CM6A-LR","MIROC-ES2L","MIROC6","MRI-ESM2-0")) {
    tx585file<-raster(paste0(path_tx585,"wc2.1_10m_tmax_",GCMs,"_ssp585_2081-2100.tif"),band=i)
    stack_mon<-stack(stack_mon,tx585file)
    
  }
  stack_mean<-mean(stack_mon)
  tx585<-stack(tx585,stack_mean)
}

path_pr585<-"/Volumes/My Passport/QT/NicheMapR/Data/Mon_Precipitation/2081-2100/ssp585/"
pr585<-stack()
for (i in 1:12) {
  stack_mon<-stack()
  for (GCMs in c("BCC-CSM2-MR","CNRM-CM6-1","CNRM-ESM2-1","CanESM5","IPSL-CM6A-LR","MIROC-ES2L","MIROC6","MRI-ESM2-0")) {
    pr585file<-raster(paste0(path_pr585,"wc2.1_10m_prec_",GCMs,"_ssp585_2081-2100.tif"),band=i)
    stack_mon<-stack(stack_mon,pr585file)
    
  }
  stack_mean<-mean(stack_mon)
  pr585<-stack(pr585,stack_mean)
}


tn585@layers<-lapply(tn585@layers,function(x)x*10)
tx585@layers<-lapply(tx585@layers,function(x)x*10)

RAINFALL <- pr585 
TMINN <-  tn585     ####Tempera
TMAXX <-  tx585
ALLMINTEMPS <- TMINN
ALLMAXTEMPS <- TMAXX
ALLTEMPS <- cbind(ALLMAXTEMPS, ALLMINTEMPS) 

ssp585_2081_2100<-stack(Altitude,RAINFALL,RAINYDAYS,WNMAXX,TMINN,TMAXX,RHMINN,RHMAXX,CCMINN)

setwd("/Volumes/My Passport/QT/NicheMapR/Data/global_climate_JZW")
writeRaster(ssp585_2081_2100, "ssp585_2081_2100.nc", overwrite=TRUE, format="CDF",     varname="Temperature", varunit="degC", 
            longname="Temperature -- raster stack to netCDF", xname="Longitude",   yname="Latitude", zname="Time (Month)")


