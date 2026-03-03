#####-----Collecting species occurrence data from sources-----#####
#####-----Zhong-wen Jiang and Liang Ma-----#####
#####-----2020-12-07-----#####
rm(list = ls())
###import libraries 这些包有的需要用，有的不需要，没有用到的就删掉
library(rgbif)
library(dplyr)
library(countrycode)
library(rnaturalearth)
#library(ggplot2)
library(rgdal)
library(raster)
library(sp)
library(dismo)


###load following data for later use, if many species build a loop, same to the 1.clim_vars.R

species<-"Phrynocephalus vlangalii"

#IUCN polygon (download IUCN polygons from IUCN redlist) #记得保存每个物种polygon的引用方???
vlangalii<-readOGR("E:/QT/SDM/P. vlangalii/data_P.shp") #Phrynocephalus vlangalii 用IUCN的polygon来限制分布点，只保留落在polygon中的分布点

#IUCN elevation
##use the elevation upper and lower data to limit your occurrences，this elevation data could from IUCN，or other database

Species_ele<-read.csv("E:/QT/SDM/Species_ele.csv") #Here based IUCNredlist，built a list，the colnames is Species, Lower, Upper，IUCN也可以提供物种分布海拔信息，可以根据海拔信息进一步筛选分布点
colnames(Species_ele)<-c("Species","Lower","Upper")
#Global DEM data
#Fick, S.E. and R.J. Hijmans, 2017. WorldClim 2: new 1km spatial resolution climate surfaces for global land areas. International Journal of Climatology 37 (12): 4302-4315.

ele<-raster("E:/QT/SDM/wc2.1_30s_elev.tif") ##global DEM data，could download from Worldclim，use the best resolution。分辨率你自己确定吧，分辨率越高，算力越高
###import occurrence data

  
  ##get occurrence records from all sources
  #a matrix for holding lonlats from sources
  lonlat_all<-matrix(NA,nrow=0,ncol=4)
  lonlat_all<-as.data.frame(lonlat_all)
  colnames(lonlat_all)<-c("Longitude","Latitude","Species","Source")
  
  #Records from survey
  lonlat_survey<-read.csv(paste("E:/QT/SDM/Survey/","Phrynocephalus vlangalii",".csv",sep=""))[,2:3]
  lonlat_survey<-cbind(lonlat_survey,rep("Phrynocephalus vlangalii",nrow(lonlat_survey)))
  lonlat_survey<-cbind(lonlat_survey,rep("Survey",nrow(lonlat_survey)))
  colnames(lonlat_survey)<-c("Longitude","Latitude","Species","Source")
  lonlat_all<-na.omit(rbind(lonlat_all,lonlat_survey))
  
  
  #Records from GBIF
    dat_GBIF<-occ_search(scientificName="Phrynocephalus vlangalii",return="data",hasCoordinate=T)
    if(!is.null(nrow(dat_GBIF$data))){
      
      dat_GBIF<-data.frame(dat_GBIF$data)
      
      #remove records without year data
      dat_GBIF<-dat_GBIF[!is.na(dat_GBIF$year),]
      
      #Keep GBIF records collected after 1979
      dat_GBIF<-dat_GBIF[dat_GBIF$year>1979,]
      
      lonlat_GBIF<-cbind(dat_GBIF$decimalLongitude,dat_GBIF$decimalLatitude)
      lonlat_GBIF<-cbind(lonlat_GBIF,rep("Phrynocephalus vlangalii",nrow(lonlat_GBIF)))
      lonlat_GBIF<-cbind(lonlat_GBIF,rep("GBIF",nrow(lonlat_GBIF)))
      colnames(lonlat_GBIF)<-c("Longitude","Latitude","Species","Source")
      lonlat_all<-rbind(lonlat_all,lonlat_GBIF)
    }
    
    

  
  #Records from literature
  lonlat_liter<-read.csv(paste("E:/QT/SDM/Literature/","Phrynocephalus vlangalii",".csv",sep=""))[,2:3]
  lonlat_liter<-cbind(lonlat_liter,rep("Phrynocephalus vlangalii",nrow(lonlat_liter)))
  lonlat_liter<-cbind(lonlat_liter,rep("Literature",nrow(lonlat_liter)))
  colnames(lonlat_liter)<-c("Longitude","Latitude","Species","Source")
  lonlat_all<-rbind(lonlat_all,lonlat_liter)
  
  lonlat_all$Longitude<-as.numeric(lonlat_all$Longitude)
  lonlat_all$Latitude<-as.numeric(lonlat_all$Latitude)
  
  
  ##Filtering data
  
  #remove sites out of IUCN polygon
  if(nrow(lonlat_all)>1){
    
    lonlat_all$Longitude<-as.numeric(lonlat_all$Longitude)
    lonlat_all$Latitude<-as.numeric(lonlat_all$Latitude)
    lonlat_all<-na.omit(lonlat_all)
    lonlat_sub<-lonlat_all
    coordinates(lonlat_sub) <- ~ Longitude + Latitude
    proj4string(lonlat_sub) <- vlangalii@proj4string
    in_out<-!is.na(over(lonlat_sub,vlangalii)[,1]) #a True/False list indicating in/out of the IUCN range
    lonlat_all<-lonlat_all[in_out,]
  }
  
  #remove sites out of IUCN elevation range
  Species_ele$Lower<-as.numeric(Species_ele$Lower)
  Species_ele$Upper<-as.numeric(Species_ele$Upper)
  
  if(nrow(lonlat_all)>1){
    ele_all<-raster::extract(ele,lonlat_all[,1:2],sp=TRUE)
    
    ele_lower<-Species_ele$Lower[Species_ele$Species=="Phrynocephalus vlangalii"]-100 #we set a buffer of 100 m here
    ele_upper<-Species_ele$Upper[Species_ele$Species=="Phrynocephalus vlangalii"]+100
    
    if(!is.na(ele_lower)&!is.na(ele_upper))lonlat_all<-lonlat_all[ele_all>ele_lower&ele_all<ele_upper,]
    if(is.na(ele_lower)&!is.na(ele_upper))lonlat_all<-lonlat_all[ele_all<ele_upper,]
    if(!is.na(ele_lower)&is.na(ele_upper))lonlat_all<-lonlat_all[ele_all>ele_lower,]
  }
  
  #remove duplicated records and reduce sampling bias
  
  if(nrow(lonlat_all)>1){
    #remove duplicated records
    lonlat_all<-na.omit(distinct(lonlat_all))
    
    #reduce sampling bias
    occ<-SpatialPoints(coords=lonlat_all[,1:2],proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    #create a RasterLayer with the extent of acgeo
    r<-raster(occ)
    #Set the resolution of the cells to 30s
    res(r)<-0.0083*5
    #Expand (extend) the extent of the RasterLayer a little
    r<-extend(r, extent(r)+0.0083*5)
    #Subsample:
    occ<-gridSample(occ, r, n=1)
    lonlat_all<-merge(lonlat_all,occ,by=c("Longitude","Latitude"))
  }
  

 write.csv(lonlat_all,file=paste("E:/QT/SDM/SpeciesRecord/","Phrynocephalus vlangalii",".csv",sep=""))



vlangalii<-read.csv("E:/QT/SDM/SpeciesRecord/Phrynocephalus vlangalii.csv") #Phrynocephalus vlangalii


plot(vlangalii$Longitude,vlangalii$Latitude)



#Omit this step if you just consider one population
#####It is diveded into three populations by elevation
#####2600m:1900-2800m;3400m:2800-3700m;4200m:3700-4600m

ele<-raster("E:/QT/SDM/wc2.1_10m_elev.tif") 
Species_ele<-read.csv("E:/QT/SDM/Species_ele.csv")[4,]
colnames(Species_ele)<-c("Species","Lower","Upper")
Lizard_coordinate<-read.csv("E:/QT/SDM/SpeciesRecord/Phrynocephalus vlangalii.csv") #Phrynocephalus vlangalii
Lizard_coordinate["Elevation"]<-raster::extract(ele,Lizard_coordinate[,2:3],sp=TRUE)


for (limit in c("3600m","3700m")) {
  Lizard_low<-subset(Lizard_coordinate,Elevation<limit)
Lizard_High<-subset(Lizard_coordinate,Elevation>=limit)

write.csv(Lizard_low,file=paste("E:/QT/SDM/",limit,"/SpeciesRecord/","Lizard_low",".csv",sep=""))
write.csv(Lizard_High,file=paste("E:/QT/SDM/",limit,"/SpeciesRecord/","Lizard_High",".csv",sep=""))



}

