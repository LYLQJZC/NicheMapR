#####-----Ensemble environmental data for SDMs-----#####
#####-----Zhong-wen Jiang -----#####
#####-----2021-12-25-----#####


####Import packages and functions####
library(rgdal)
library(raster) # Raster data management
library(stringr)
library(units)
library(rlist)
library(rgeos)
library(sf)

rasterOptions(maxmemory=2.5e+10,tmpdir="E:/TempR/")


####Import IUCN polygons###
#download from IUCN(detail in supporting information 1)
vlangalii<-readOGR("E:/QT/SDM/P. vlangalii/data_P.shp") #Phrynocephalus vlangalii
vlangalii_buffer<-gBuffer(vlangalii,width=1) ####add a 100km buffer


####Species to model###
#if here have many species，need to build a loop such as for(spe_list in c("spe1","spe2"...)){}
spe_list<-"Phrynocephalus vlangalii"

####1. Prepare environmental data####
###1.1 Current environmental data
###1.1.1 Current climate from worldclim (https://worldclim.org) 1970-2000 10min BIO1-19 (detail in supporting information2)
### it's same for other macroclimate databases, such as CHELSA, which is  closer to current time 

path_wc_1970_2000<-"E:/QT/SDM/CMIP6/Current_10m/" ##the path of your downloaded files
wc_1970_2000<-list.files(path_wc_1970_2000,pattern="*.tif")
wc_1970_2000<-str_sort(sapply(wc_1970_2000,function(x)x[1]), numeric = TRUE)###order the vairable from bio1-bio19, here you need to make sure the default order is right 
wc_1970_2000<-lapply(wc_1970_2000,function(x) paste(path_wc_1970_2000,x,sep=""))
wc_1970_2000<-lapply(wc_1970_2000,raster)
env_1970_2000<-do.call(stack,wc_1970_2000)
  
path="E:/QT/SDM/United/"###set the output path by yourself

env_1970_2000_crop<-crop(env_1970_2000,vlangalii_buffer)
env_1970_2000_crop<-mask(env_1970_2000_crop,vlangalii_buffer)

writeRaster(env_1970_2000_crop, filename=paste(path,spe_list,"env_1970_2000_10m.tif",sep=""), format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))

##(Omit this step，if you just consider climate variables)
###1.1.2 Ecological consequence data (ACT/WT/CT/WL/QMET/END/OC/RWL) from NicheMapR 10min
en_1970_2000<-stack(paste(path,spe_list,"env_1970_2000_10m.tif",sep=""))
ref<-raster(paste(path,spe_list,"env_1970_2000_10m.tif",sep=""),band=1)

for (pop in c("Low","High")) {
  
 # pop="High"
 
  path_ec<-paste0("E:/QT/NicheMapR/Result/",pop,"_pop/1970_2000/")
  
  ec_1970_2000<-list.files(path_ec,pattern="*.tif")
  ec_1970_2000<-lapply(ec_1970_2000,function(x) paste(path_ec,x,sep=""))
  ec_1970_2000<-lapply(ec_1970_2000,raster)
  ec_1970_2000<-lapply(ec_1970_2000, function(x)resample(x,ref))
  ecc_1970_2000<-do.call(stack,ec_1970_2000)
  
  con_var_1970_2000<-stack(en_1970_2000,ecc_1970_2000)
  
  writeRaster(con_var_1970_2000, filename=paste("E:/QT/SDM/United/Niche_var/Eco_",pop,"_convar_1970_2000_10m.tif",sep=""), format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
  
}

##(Omit this step，if you just consider climate variables)
###1.1.3 Adaptation Ecological consequence data (ACT/WT/CT/WL/QMET/END/OC/RWL) from NicheMapR 10min

  path_ec<-paste0("E:/QT/NicheMapR/Result/Adaptation_1970_2000/")
  
  ec_1970_2000<-list.files(path_ec,pattern="*.tif")
  ec_1970_2000<-lapply(ec_1970_2000,function(x) paste(path_ec,x,sep=""))
  ec_1970_2000<-lapply(ec_1970_2000,raster)
  ec_1970_2000<-lapply(ec_1970_2000, function(x)resample(x,ref))
  ecc_1970_2000<-do.call(stack,ec_1970_2000)
  
  con_var_1970_2000<-stack(en_1970_2000,ecc_1970_2000)
  
  writeRaster(con_var_1970_2000, filename=paste("E:/QT/SDM/United/Niche_var/Eco_convar_1970_2000_10m.tif",sep=""), format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
  


##1.2 Future environmental data
#1.2.1 Ensemble environmental variables for future condition by average GCMs
###Future climate from worldclim (https://worldclim.org) 10min BIO1-19 (detail in supporting information 3)
time="2081-2100" ##depend on your study periods
for(ssp in c("ssp126","ssp585")){#depend on your study climate scenarios
  Clim_bio<-raster()
  for (bio in 1:19) {
    Bio_stack<-raster()
   # bio=1
    for(mod in c("BCC-CSM2-MR","CNRM-CM6-1","CNRM-ESM2-1","CanESM5","IPSL-CM6A-LR","MIROC-ES2L","MIROC6","MRI-ESM2-0")){ #select some GCMs, who perform better on your study regions
      # ssp<-"ssp245"
      #mod<-"BCC-CSM2-MR"
      ##1.3.2.1 Future climate
      path_wc<-paste("E:/QT/SDM/CMIP6/",time,"/",ssp,"/",mod,"/",sep="")
      wc_2081_2100<-list.files(path_wc,pattern="*.tif")
      path_wc_2081_2100<-lapply(wc_2081_2100, function(x) paste(path_wc,x,sep=""))
      # wc_2061_2080<-str_sort(sapply(wc_2061_2080,function(x)x[1]), numeric = TRUE)###order the vairable from bio1-bio19
      wc_2081_2100<-raster(path_wc_2081_2100[[1]],band=bio)
      Bio_stack<-stack(Bio_stack,wc_2081_2100)
    }
    Bio_mean<-mean(Bio_stack)
    names(Bio_mean)<-paste0("Bio",bio)
    Clim_bio<-stack(Clim_bio,Bio_mean)
    }
       Clim_bio_crop<-crop(Clim_bio,vlangalii_buffer)
       Clim_bio_crop<-mask(Clim_bio_crop,vlangalii_buffer)
    
  
    writeRaster(Clim_bio_crop, filename=paste(path,spe_list,"env_",ssp,"_2081_2100_10m.tif",sep=""), format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
    
    print(ssp)  
}

##(Omit this step，if you just consider climate variables)
#1.2.2 Ensemble environmental variables and ecological consequences for future condition 
for(ssp in c("ssp126","ssp585")){
  Clim_data<-stack(paste(path,spe_list,"env_",ssp,"_2081_2100_10m.tif",sep=""))
  
  for (pop in c("Low","High")) {
    path_ec<-paste0("E:/QT/NicheMapR/Result/",pop,"_pop/",ssp,"_2081_2100/")
    
    ec_2081_2100<-list.files(path_ec,pattern="*.tif")
    ec_2081_2100<-lapply(ec_2081_2100,function(x) paste(path_ec,x,sep=""))
    ec_2081_2100<-lapply(ec_2081_2100,raster)
    ec_2081_2100<-lapply(ec_2081_2100, function(x)resample(x,ref))
    ecc_2081_2100<-do.call(stack,ec_2081_2100)
    
    con_var_2081_2100<-stack(Clim_data,ecc_2081_2100)
    
    writeRaster(con_var_2081_2100, filename=paste("E:/QT/SDM/United/Niche_var/Eco_",pop,"_",ssp,"_convar_2081_2100_10m.tif",sep=""), format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
  }
  
}

##(Omit this step，if you just consider climate variables)
#1.2.3 Ensemble environmental variables and adaptation consequences for future condition 
for(ssp in c("ssp126","ssp585")){
  Clim_data<-stack(paste(path,spe_list,"env_",ssp,"_2081_2100_10m.tif",sep=""))
  
 
    path_ac<-paste0("E:/QT/NicheMapR/Result/Adaptation_",ssp,"/")
    
    ac_2081_2100<-list.files(path_ac,pattern="*.tif")
    ac_2081_2100<-lapply(ac_2081_2100,function(x) paste(path_ac,x,sep=""))
    ac_2081_2100<-lapply(ac_2081_2100,raster)
    ac_2081_2100<-lapply(ac_2081_2100, function(x)resample(x,ref))
    acc_2081_2100<-do.call(stack,ac_2081_2100)
    
    con_vaar_2081_2100<-stack(Clim_data,acc_2081_2100)
    
    writeRaster(con_vaar_2081_2100, filename=paste("E:/QT/SDM/United/Niche_var/Eco_",ssp,"_convar_2081_2100_10m.tif",sep=""), format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
 
  
}



