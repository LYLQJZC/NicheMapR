#############Lizard_mod (ActiveTime/WarmingTolerance/WaterLoss/OxygenConsumption)###############
##################################Zhong-wen Jiang 2021/12/28##############################

rm(list = ls())
library(NicheMapR)
library(raster)
clmfolder<-paste("E:/QT/NicheMapR/Data/global_climate_JZW")
global_climate<-raster::brick(paste(clmfolder,"/global_climate.nc",sep=""))
soilmoisture<-raster::brick(paste(clmfolder,"/soilw.mon.ltm.v2.nc",sep=""))
source("E:/QT/NicheMapR/Code/micro_worldclim_version/3_micro_worldclim.R")

Lizard_traits<-read.csv("E:/QT/NicheMapR/Data/Lizard_traits.csv")
Lizard_extent<-raster("E:/QT/SDM/United/Phrynocephalus vlangaliienv_1970_2000_10m.tif")
ele<-raster("E:/QT/SDM/wc2.1_10m_elev.tif")
slope<-terrain(ele,'slope',unit='degrees',neighbors=8)
aspect<-terrain(ele,'aspect',unit='degrees',neighbors=8)

Lizard_point<-as.data.frame(rasterToPoints(Lizard_extent))[1:2]
colnames(Lizard_point)<-c("lon","lat")
Lizard_point["slope"]<-extract(slope,Lizard_point[,(1:2)])
Lizard_point["aspect"]<-extract(aspect,Lizard_point[,(1:2)])




#####lizard_mod function####

cal.liz<-function(micro,Ele){
  
  trait_list<-Lizard_traits[which(Lizard_traits$﻿Ele==Ele),]
  
  lizard_loc<-ectotherm(
    live=1,
    Ww_g	= trait_list$Ww_g, #default=40; Wet weight of animal (g), note this model is 'steady state' so no lags in heating/cooling due to mass
    alpha_max	= 0.95, #default=0.85; Maximum solar absorptivity, 0-1; keep this variable for sensitivity analysis
    alpha_min	= 0.95, #default=0.85; Minimum solar absorptivity, 0-1; keep this variable for sensitivity analysis
    T_F_min	= trait_list$T_F_min, #default=24; Minimum foraging temperature, °C (also affects burrow depth selection)
    T_F_max	= trait_list$T_F_max, #default=34; Maximum foraging temperature, °C
    T_B_min = trait_list$T_B_min, #default=17.5; Minimum basking temperature, °C
    T_RB_min	= trait_list$T_B_min, #default=17.5; Minimum temperature at which animal will move from retreat to basking site, °C
    T_pref	= trait_list$T_pref, #default=30; Preferred body temperature, °C
    CT_max	= trait_list$CT_max, #default=40; Critical thermal maximum, °C (affects burrow depth selection and may be used to impose death from heat stress)
    CT_min	= trait_list$CT_min, #default=6; Critical thermal minimum, °C (affects burrow depth selection and may be used to impose death from cold stress)
    
    M_1	= trait_list$M_1, #Metabolic rate parameter 1 V_O2=M_1*M^M_2*10^(M_3*Tb), in ml O2 / h, default parameters for lizards based on Eq. 2 from Andrews & Pough 1985. Physiol. Zool. 58:214-231
    M_2	= trait_list$M_2, #Metabolic rate parameter 2
    M_3	= trait_list$M_3, #Metabolic rate parameter 3
    
    
    pct_wet=1,#%of surface area acting as a free-water exchanger, for computing cutaneous water loss,calculate by get_p_wet/ data from (Marco 2018 Zoology)
    #pct_eyes=0.03,#% of surface area taken up by open eyes, for computing ocular water loss (only when active)
    #pct_mouth=5,#% of surface area taken up by open mouth, for computing panting water loss
    #pantmax=1,# maximum multiplier on breathing rate, for respiratory water loss via panting
    #F_O2=20,#% oxygen extraction efficiency, for respiratory water loss
    #delta_air=0.1,#Temperature difference (??C) between expired and inspired air, for computing respiratory water loss
    mindepth=1,
    maxdepth=9,
    aestdepth	= 9, #Depth (soil node #) to which animal retreats if burrowing and aestivating due to desiccation
    #transient = 1,#default=0,Run a transient (i.e. include heat storage) simulation (1=yes, 0=no)? No behaviour yet - assums full sun
    
    #Details
    #Parameters controling how the model runs:
    
    nyears = 1, #micro$nyears, #Number of years the simulation runs for - must be consistent with dimensions of environmental input data
    
    #Environmental inputs:
    minshades =micro[["minshade"]], #Vector of daily minimum shade values - can be different to value used in microclimate model (e.g. to simulate sunspot tracking on a forest floor) (%)
    maxshades =micro[["maxshade"]], #Vector of daily maximum shade values - can be different to value used in microclimate model (e.g. to simulate use of fine-scale shade in a generally unshaded habitat) (%)
    
    alpha_sub = 1-micro[["REFL"]], #Vector of daily substrate reflectances (0-1)?ù̶?ɳ?ӵ?
    DEP = micro$DEP, #Depths available from the microclimate model simulation
    metout = micro$metout, #Microclimate model output for above ground, minimum shade conditions
    shadmet = micro$shadmet, #Microclimate model output for above ground, maximum shade conditions
    soil = micro$soil, #Microclimate model output for soil temperature, minimum shade conditions
    shadsoil = micro$shadsoil, #Microclimate model output for soil temperature, maximum shade conditions
    soilmoist = micro$soilmoist, #Microclimate model output for soil moisture, minimum shade conditions
    shadmoist = micro$shadmoist, #Microclimate model output for soil moisture, maximum shade conditions
    humid = micro$humid, #Microclimate model output for soil humidity, minimum shade conditions
    shadhumid = micro$shadhumid, #Microclimate model output for soil humidity, maximum shade conditions
    soilpot = micro$soilpot, #Microclimate model output for soil water potential, minimum shade conditions
    shadpot = micro$shadpot, #Microclimate model output for soil water potential, maximum shade conditions
    rainfall = as.integer(micro$RAINFALL), #Vector of daily rainfall (mm)
    rainhr = rep(-1,length(micro$metout)), #Vector of hourly rainfall (mm), overwrites rainfall if not negative
    elev = as.numeric(micro$elev), #Elevation of simulation (m), obtained from microclimate model output by default
    longitude = micro$longlat[1], #Longitude (decimal degrees), obtained from microclimate model output by default
    latitude = micro$longlat[2], #Latitude (decimal degrees), obtained from microclimate model output by default
    shade_seek =0,
    burrow = 1
    #Morphological parameters:
    
    #custom_shape = c(10.4713,.688,0.425,0.85,3.798,.683,0.694,.743), Custom shape coefficients. Operates if shape=5, and consists of 4 pairs of values representing the parameters a and b of a relationship AREA=a*mass^b, where AREA is in cm2 and mass is in g. The first pair are a and b for total surface area, then a and b for ventral area, then for sillhouette area normal to the sun, then sillhouette area perpendicular to the sun
    #shape_a = 1, Proportionality factor (-) for going from volume to area, keep this 1 (redundant parameter that should be removed)
    #shape_b = 3, Proportionality factor (-) for going from volume to area, represents ratio of width:height for a plate, length:diameter for cylinder, b axis:a axis for ellipsoid 
    #shape_c = 0.6666666667, Proportionality factor (-) for going from volume to area, represents ratio of length:height for a plate, c axis:a axis for ellipsoid
    #fatosk = 0.4, Configuration factor to sky (-) for infrared calculations
    #fatosb = 0.4, Configuration factor to subsrate for infrared calculations
    #rinsul = 0, Insulative fat layer thickness (m)
    #pct_cond = 10, Percentage of surface contacting the substrate (%)
    #c_body = 3073, #Specific heat of flesh J/(kg-K)
    #k_flesh = 0.5, #Thermal conductivity of flesh (W/mC, range: 0.412-2.8)
    #rho_body = 1000, Density of flesh (kg/m3)
    #epsilon = 0.95, Emissivity of animal (0-1)
    
    #Behavioural parameters:
  
    #postur = 1, postural orientation to sun, 1 = perpendicular, 2 = parallel, 0 = half way between, relevant if live = 0
    #warmsig = 0, Warming signal for emergence? °C/h (if in burrow deeper than node 2, change in burrow temp must be exceed warmsig)
    #fossorial = 0, Fossorial activity? 1=yes, 0=no (this option hasn't been properly implemented)
    #rainact =1 , #Activity is limited by rainfall? 1=yes, 0=no, threshold rainfall for activity set by actrainthresh
    #actrainthresh = 0.1, #Threshold (mm) of rain causing activity if rainact=1
    #soilnode = 4, #Soil node (1-10, corresponding to values in DEP) at which eggs are laid (overridden if amphibreed=1)
    #eggshade = 0, #are eggs laid in shade? 0=no, 1=yes
    #aquabask = 0, #If aquatic, does it bask? 0=no, stay at water temp, 1=yes, when not hungry, 2=all the time

    #Thermal physiological parameters:
    
    #CT_minthresh = 12, Number of consecutive hours below CT_min that leads to death - simulation will terminate beyond this threshold if CT_kill=1
    #CT_kill = 0, Animal dies when it hits critical thermal limits? 1=yes, 0=no
  )

  lizard_environ_out<- data.frame(lizard_loc[["environ"]])
  lizard_enbal_out<- data.frame(lizard_loc[["enbal"]])
  lizard_masbal_out<- data.frame(lizard_loc[["masbal"]])

  # Environ<-data.frame("DOY"=lizard_loc_out$environ.DOY,"DAY"=lizard_loc_out$environ.DAY,"TIME"=lizard_loc_out$environ.TIME,"TC"=lizard_loc_out$environ.TC,"ACT"=lizard_loc_out$environ.ACT,"Ele"= as.numeric(micro$elev))
  
  ACT<-nrow(subset(lizard_environ_out,ACT!=0))
  
  Ele<-as.numeric(micro$elev)
  Coordinate<-data.frame("longitude" = micro$longlat[1],"latitude" = micro$longlat[2])
  
  #365*24
  
  Result<-data.frame("longitude" = micro$longlat[1],"latitude" = micro$longlat[2],"ele"=as.numeric(micro$elev),"ACT"=nrow(subset(lizard_environ_out,ACT!=0)),"Tbmax"=max(lizard_environ_out$TC),"Tbmin"=min(lizard_environ_out$TC),"ENB"=sum(lizard_enbal_out$ENB),"OC"=sum(lizard_masbal_out$O2_ml),"QMET"=sum(lizard_enbal_out$QMET),"WL"=sum(lizard_masbal_out$H2OResp_g+lizard_masbal_out$H2OCut_g+lizard_masbal_out$H2OEye_g),"RWL"=sum(lizard_masbal_out$H2OResp_g))
  
  Result
} 

#time="1970_2000"

for (time in c("1970_2000","ssp126_2081_2100","ssp585_2081_2100")) {
  
  
  for (Ele in c("3400m","4200m")) {
    
    result<-data.frame()
    
    for(p in c(1:nrow(Lizard_point))){
      
      tryCatch(
        
        expr = {
          
          
          micro<-micro_world(loc = c(Lizard_point[p,]$lon,Lizard_point[p,]$lat),slope =Lizard_point[p,]$slope,aspect = Lizard_point[p,]$aspect,
                             run.gads = 2,timeinterval = 365,soiltype = 1,REFL = 0.35,CMH2O = 0.5,maxshade = 10)
          
          a<-cal.liz(micro,Ele)
          
          
          result<-rbind(result,a)
          
          print(p)
          
          write.csv(result, file=paste("E:/QT/NicheMapR/Result/Eco_Result_",time,"_",Ele,".csv",sep=""),row.names = FALSE)
          
        },
        
        error=function(e){print(paste0(Ele,"_",p,"error"))})
      
    }
    print(paste0(Ele,"_",time))
  }
  
}

