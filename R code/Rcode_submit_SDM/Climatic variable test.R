#####Climatic variable####


####Import packages and functions####
library(rgdal)
library(raster)
library(biomod2) # Ensemble SDM package
library(usdm)
library(rJava)
library(ecospat)

memory.limit(size=31500000)
rasterOptions(maxmemory=2.5e+10,tmpdir="E:/TempR/")
####Import global variables and functions####
####Import IUCN polygons####
vlangalii<-readOGR("E:/QT/SDM/P. vlangalii/data_P.shp") #Phrynocephalus vlangalii
species="Phrynocephalus vlangalii"

#pop="Low"
for (pop in c("Low","High")) {
  
  setwd(paste0("E:/QT/SDM/SDM_Group/",pop))
  spe_IUCN<-subset(vlangalii, BINOMIAL==species) #IUCN range of i species
  path_occ="E:/QT/SDM/SpeciesRecord/"
  occ<-as.data.frame(read.csv(paste(path_occ,"Lizard_",pop,".csv",sep=""),header=T))[,3:4]
  colnames(occ)<-c("Longitude","Latitude")
  coordinates(occ) <- ~ Longitude + Latitude
  proj4string(occ) <- proj4string(spe_IUCN)
  
  
  ####Import currernt climate data####
 
  path="E:/QT/SDM/United/Niche_var/"
  env_1970_2000<-stack(paste(path,"Eco_",pop,"_convar_1970_2000_10m.tif",sep=""))
  #env_1970_2000<-stack(paste0(path,"Phrynocephalus vlangaliiHigh_convar_1970_2000_10m.tif"))
  varnames<-read.csv("E:/QT/SDM/varname1.csv")$ 锘縑1
  names(env_1970_2000)<-varnames
  env_1970_2000_sub<-env_1970_2000[[c(1:19)]]
  env_1970_2000_sub.NoCor <- vifcor(env_1970_2000_sub, th = 0.7)
  env_1970_2000.Reduced <- exclude(env_1970_2000_sub, env_1970_2000_sub.NoCor)
  ###bio2/3/4/9/14
  
  data<-values(env_1970_2000.Reduced)
  data<-na.omit(data)
  data<-as.data.frame(data)
  cor_value<-as.data.frame(cor(data,method = "pearson"))
  vif_value<-as.data.frame(vif(data))
  #write.csv(cor_value,file =paste0("E:/QT/SDM/SDM_Group/",pop,"/",pop,"_corvalue.csv"))
  #write.csv(vif_value,file =paste0("E:/QT/SDM/SDM_Group/",pop,"/",pop,"_vifvalue.csv"),row.names = FALSE )
  
  
  ####4.Initialize the datasets for usage in biomod2 (Thuiller et al. 2009; Thuiller et al. 2016)####
  mod.dat<-BIOMOD_FormatingData(resp.var=occ,
                                expl.var=env_1970_2000.Reduced,
                                resp.name=species,
                                PA.nb.rep=1,
                                PA.nb.absences=1000,
                                PA.strategy="random")
  
  myBiomodOptions <- BIOMOD_ModelingOptions(
    MAXENT.Phillips=list(
      path_to_maxent.jar='E:/QT/SDM/SDM_Group/maxent/maxent/maxent.jar'))
  #Download maxent.jar from http://www.cs.princeton.edu/~schapire/maxent/; and store it in the working directory
  #Also need to download and install java (x64)
  
  Biomod.tuning<-BIOMOD_tuning(mod.dat,
                               models=c('ANN','MAXENT.Phillips'),
                               env.ME = env_1970_2000.Reduced)###GLM为什么不能tun，邮件一下。
  
  path_rds<-"E:/QT/SDM/SDM_Group/"
  saveRDS(Biomod.tuning,paste(path_rds,"/",pop,"/",species,"_",pop,"_Biomod_tuning_clim19.rds",sep=""))
  Biomod.tuning<-readRDS(paste(path_rds,"/",pop,"/",species,"_",pop,"_Biomod_tuning_clim19.rds",sep=""))
  
  
  #'GBM','CTA','ANN','FDA','MARS','RF','MAXENT.Phillips'
  
  mod<-ecospat.ESM.Modeling(mod.dat, 
                            models=c('GLM','ANN','MAXENT.Phillips'),
                            models.options = Biomod.tuning$models.options,
                            NbRunEval=10,
                            DataSplit=80, # common practice to validate! 
                            weighting.score = c('SomersD'),
                            parallel = F)
  saveRDS(mod, paste(path_rds,"/",pop,"/",species,"_",pop,"_mod_clim19.rds",sep=""))
  
  #mod<-readRDS(paste(path_rds,"/",pop,"/",species,"_",pop,"_mod_clim.rds",sep=""))
  
  
  ####6.Evaluation and average of simple bivariate models to ESMs####
  mod.EM<-ecospat.ESM.EnsembleModeling(mod,weighting.score="SomersD")
  saveRDS(mod.EM, paste0(path_rds,pop,"/",species,"_",pop,"_mod.EM_clim19.rds"))
  
  EM_Eval<-mod.EM$ESM.evaluations
  write.csv(EM_Eval,paste0(path_rds,pop,"/",species,"_",pop,"_EM_Eval_rep10_clim19.csv"),row.names=F)
  
  EM_Contrib<-ecospat.ESM.VarContrib(mod,mod.EM)
  write.csv(EM_Contrib,paste0(path_rds,pop,"/",species,"_",pop,"_EM_Contrib_clim.csv"))
  
  #pop="Low"
  #mod.EM<-readRDS(paste0(path_rds,pop,"/",species,"_",pop,"_mod.EM_clim.rds"))
  
  threshold<-ecospat.ESM.threshold(mod.EM)
  threshold<-threshold$TSS.th[4]*1000
  print(paste0(pop,"_threshold=",threshold))
  
  ####7.Make projections under current condition####
  mod.proj.current<-ecospat.ESM.Projection(ESM.modeling.output=mod,
                                           new.env=env_1970_2000.Reduced,
                                           #parallel=T,
                                           cleanup=48)
  saveRDS(mod.proj.current, paste0(path_rds,pop,"/",species,"_",pop,"_mod.proj.current_clim19.rds"))
  
  mod.EFproj.current<-ecospat.ESM.EnsembleProjection(ESM.prediction.output=mod.proj.current,
                                                     ESM.EnsembleModeling.output=mod.EM)
  
  saveRDS(mod.EFproj.current, paste0(path_rds,pop,"/",species,"_",pop,"_mod.EFproj.current_clim19.rds"))
  
  binary_current<-ecospat.binary.model(mod.EFproj.current$EF,threshold)##threshold from optimal.thresholds
  
  writeRaster(mod.EFproj.current$EF,paste0("E:/QT/SDM/SDM_Group/",pop,"/TIFresult/","EF_current_",pop,"_clim19.tif"),format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
  writeRaster(binary_current,paste0("E:/QT/SDM/SDM_Group/",pop,"/TIFresult/","binary_current_",pop,"_clim19.tif"),format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
  
  #####.8 Make projections under future condition####
  
  
  env_2081_2100_ssp126<-stack(paste(path,"Eco_",pop,"_ssp126_convar_2081_2100_10m.tif",sep=""))
  env_2081_2100_ssp585<-stack(paste(path,"Eco_",pop,"_ssp585_convar_2081_2100_10m.tif",sep=""))
  env_2081_2100_ssp126_ad<-stack(paste(path,"Eco_ssp126_convar_2081_2100_10m.tif",sep=""))
  env_2081_2100_ssp585_ad<-stack(paste(path,"Eco_ssp585_convar_2081_2100_10m.tif",sep=""))
  
  names(env_2081_2100_ssp126)<-varnames
  names(env_2081_2100_ssp585)<-varnames
  names(env_2081_2100_ssp126_ad)<-varnames
  names(env_2081_2100_ssp585_ad)<-varnames
  
  env_2081_2100_ssp126_Reduce<-exclude(env_2081_2100_ssp126,env_1970_2000_sub.NoCor)
  env_2081_2100_ssp585_Reduce<-exclude(env_2081_2100_ssp585,env_1970_2000_sub.NoCor)
  env_2081_2100_ssp126ad_Reduce<-exclude(env_2081_2100_ssp126_ad,env_1970_2000_sub.NoCor)
  env_2081_2100_ssp585ad_Reduce<-exclude(env_2081_2100_ssp585_ad,env_1970_2000_sub.NoCor)
  

  #####Just climate ssp126
  
  mod.proj_wc<-ecospat.ESM.Projection(ESM.modeling.output=mod,
                                      new.env=env_2081_2100_ssp126_Reduce)
  #eval(parse(text = paste("env_2081_2100_",ssp,"_Reduce",sep="")))
  saveRDS(mod.proj_wc, paste(path_rds,"/",pop,"/",species,"_",pop,"_mod.proj_ssp126_clim19.rds",sep=""))
  
  mod.proj_wc$pred.biva<-mod.proj_wc$pred.biva[grep("wclu",mod.proj_wc$pred.biva,invert=T)]
  
  mod.EFproj_wc<-ecospat.ESM.EnsembleProjection(ESM.prediction.output=mod.proj_wc,
                                                ESM.EnsembleModeling.output=mod.EM)
  saveRDS(mod.EFproj_wc, paste(path_rds,"/",pop,"/",species,"_",pop,"_mod.EFproj_ssp126_clim19.rds",sep=""))
  
  binary_ssp126_wc<-ecospat.binary.model(mod.EFproj_wc$EF,threshold)##threshold from optimal.thresholds
  writeRaster(mod.EFproj_wc$EF,paste0("E:/QT/SDM/SDM_Group/",pop,"/TIFresult/","EF_ssp126_clim19_",pop,".tif"),format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
  writeRaster(binary_ssp126_wc,paste0("E:/QT/SDM/SDM_Group/",pop,"/TIFresult/","binary_ssp126_clim19_",pop,".tif"),format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
  
  
 
  
  #####Just climate ssp585
  mod.proj_wc<-ecospat.ESM.Projection(ESM.modeling.output=mod,
                                      new.env=env_2081_2100_ssp585_Reduce)
  #eval(parse(text = paste("env_2081_2100_",ssp,"_Reduce",sep="")))
  saveRDS(mod.proj_wc, paste(path_rds,"/",pop,"/",species,"_",pop,"_mod.proj_ssp585_clim19.rds",sep=""))
  
  mod.proj_wc$pred.biva<-mod.proj_wc$pred.biva[grep("wclu",mod.proj_wc$pred.biva,invert=T)]
  
  mod.EFproj_wc<-ecospat.ESM.EnsembleProjection(ESM.prediction.output=mod.proj_wc,
                                                ESM.EnsembleModeling.output=mod.EM)
  saveRDS(mod.EFproj_wc, paste(path_rds,"/",pop,"/",species,"_",pop,"_mod.EFproj_ssp585_clim19.rds",sep=""))
  
  
  binary_ssp585_wc<-ecospat.binary.model(mod.EFproj_wc$EF,threshold)##threshold from optimal.thresholds
  writeRaster(mod.EFproj_wc$EF,paste0("E:/QT/SDM/SDM_Group/",pop,"/TIFresult/","EF_ssp585_clim19_",pop,".tif"),format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
  writeRaster(binary_ssp585_wc,paste0("E:/QT/SDM/SDM_Group/",pop,"/TIFresult/","binary_ssp585_clim19_",pop,".tif"),format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
  
  print(pop)
}



