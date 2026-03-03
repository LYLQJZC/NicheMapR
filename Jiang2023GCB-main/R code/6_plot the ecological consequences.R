######plot the ecological consequences#####
library(raster)


limit="3600m"


for (limit in c("3600m","3700m")) {
  
  

    path="E:/QT/SDM/United/Niche_var/"
    path_occ=paste("E:/QT/SDM/",limit,"/SpeciesRecord/",sep="")
    
    varnames<-read.csv("E:/QT/SDM/varname1.csv")$ ﻿V1
#  pop="Low"

  
  all_traits<-data.frame()
  all_change<-data.frame()
  

for (pop in c("Low","High")) {
  
  occ<-as.data.frame(read.csv(paste(path_occ,"Lizard_",pop,".csv",sep=""),header=T))[,3:4]
  colnames(occ)<-c("Longitude","Latitude")
  occ["Pop"]<-pop
  
  for (i in c(20,23,26,27)) {
  #i=20

    occ1<-occ
    occ2<-occ
    occ3<-occ
    occ4<-occ
    occ5<-occ
    
    env_1970_2000<-stack(paste(path,"Eco_",pop,"_convar_1970_2000_10m.tif",sep=""))
    env_2081_2100_ssp126<-stack(paste(path,"Eco_",pop,"_ssp126_convar_2081_2100_10m.tif",sep=""))
    env_2081_2100_ssp585<-stack(paste(path,"Eco_",pop,"_ssp585_convar_2081_2100_10m.tif",sep=""))
    names(env_1970_2000)<-varnames
    names(env_2081_2100_ssp126)<-varnames
    names(env_2081_2100_ssp585)<-varnames
  
    occ1["Time"]<-"Current"
    occ1["Traits"]<-names(env_1970_2000[[i]])
    occ1["Value"]<-extract(env_1970_2000[[i]],occ1[,c(1:2)])
    occ2["Time"]<-"ssp126"
    occ2["Traits"]<-names(env_2081_2100_ssp126[[i]])
    occ2["Value"]<-extract(env_2081_2100_ssp126[[i]],occ2[,c(1:2)])
    occ3["Time"]<-"ssp585"
    occ3["Traits"]<-names(env_2081_2100_ssp585[[i]])
    occ3["Value"]<-extract(env_2081_2100_ssp585[[i]],occ3[,c(1:2)])
    
    occ4["Time"]<-"ssp126_change"
    occ4["Traits"]<-names(env_2081_2100_ssp126[[i]])
    occ4["Value"]<-extract(env_2081_2100_ssp126[[i]],occ4[,c(1:2)])-extract(env_1970_2000[[i]],occ4[,c(1:2)])
    occ5["Time"]<-"ssp585_change"
    occ5["Traits"]<-names(env_2081_2100_ssp585[[i]])
    occ5["Value"]<-extract(env_2081_2100_ssp585[[i]],occ5[,c(1:2)])-extract(env_1970_2000[[i]],occ5[,c(1:2)])
    
    
    result1<-rbind(occ1,occ2,occ3)
    result2<-rbind(occ4,occ5)
      
    all_traits<-rbind(all_traits,result1)
    all_change<-rbind(all_change,result2)
    
  }
}
  
  
  all_traits
  
  library(ggplot2)
  
 Three_traits<-all_traits[all_traits$Traits!="WT",]
 Three_change<-all_change[all_change$Traits!="WT",]
 
 Three_traits_Low<-Three_traits[Three_traits$Pop=="Low",]
 Three_traits_High<-Three_traits[Three_traits$Pop=="High",]
 
 
 
  fig_eco<-ggplot(data=Three_traits,aes(x=Time,y=Value,fill=Pop))+
    scale_fill_manual(values=c("#4169E1","#DC143C"))+
    geom_boxplot()+
    facet_wrap(Traits~.,scales = "free",nrow=1)+
    #facet_grid(.~Traits,scales = "free")+
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
  
    
  ggsave(fig_eco,file=paste("E:/QT/SDM/",limit,"/SDM_Group/Fig_3.pdf",sep=""),device = "pdf",height=4,width=12,dpi = 300)
  
  
  
  
  
  
  
ACT<-Three_traits[Three_traits$Traits=="ACT",]
OC<-Three_traits[Three_traits$Traits=="OC",]
WL<-Three_traits[Three_traits$Traits=="WL",]

library(car)
leveneTest(ACT$Value~ACT$Pop)
mod_ACT<-aov(Value~Time+Pop,data = ACT)
summary(mod_ACT)
shapiro.test(resid(mod_ACT))


kruskal.test(Value~Pop,data = ACT)
kruskal.test(Value~Time,data = ACT)




ACT_current<-ACT[ACT$Time=="Current",]
ACT_ssp126<-ACT[ACT$Time=="ssp126",]
ACT_ssp585<-ACT[ACT$Time=="ssp585",]

kruskal.test(Value~Pop,data = ACT_current)
kruskal.test(Value~Pop,data = ACT_ssp126)
kruskal.test(Value~Pop,data = ACT_ssp585)


leveneTest(log(OC$Value)~OC$Pop)
mod_OC<-lm(log(Value)~Time*Pop,data = OC)
summary(mod_OC)
shapiro.test(resid(mod_OC))  

kruskal.test(Value~Pop,data = OC)
kruskal.test(Value~Time,data = OC)

OC_current<-OC[OC$Time=="Current",]
OC_ssp126<-OC[OC$Time=="ssp126",]
OC_ssp585<-OC[OC$Time=="ssp585",]

kruskal.test(Value~Pop,data = OC_current)
kruskal.test(Value~Pop,data = OC_ssp126)
kruskal.test(Value~Pop,data = OC_ssp585)


leveneTest(WL$Value~WL$Pop)
mod_WL<-lm(Value~Time*Pop,data = WL)
summary(mod_WL)
shapiro.test(resid(mod_WL))  


kruskal.test(Value~Pop,data = WL)
kruskal.test(Value~Time,data = WL)


WL_current<-WL[WL$Time=="Current",]
WL_ssp126<-WL[WL$Time=="ssp126",]
WL_ssp585<-WL[WL$Time=="ssp585",]

kruskal.test(Value~Pop,data = WL_current)
kruskal.test(Value~Pop,data = WL_ssp126)
kruskal.test(Value~Pop,data = WL_ssp585)







fig_ecochange<-ggplot(data=Three_change,aes(x=Time,y=Value,fill=Pop))+
  scale_fill_manual(values=c("#4169E1","#DC143C"))+
  geom_boxplot()+
  facet_wrap(Traits~.,scales = "free",nrow=1)+
  #facet_grid(.~Traits,scales = "free")+
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


ggsave(fig_ecochange,file=paste("E:/QT/SDM/",limit,"/SDM_Group/Fig_3_2.pdf",sep=""),device = "pdf",height=4,width=12,dpi = 300)





}




ACT_chan<-Three_change[Three_change$Traits=="ACT",]
OC_chan<-Three_change[Three_change$Traits=="OC",]
WL_chan<-Three_change[Three_change$Traits=="WL",]

library(car)
leveneTest(ACT_chan$Value~ACT_chan$Pop)
mod_ACT<-lm(Value~Time+Pop,data = ACT_chan)
summary(mod_ACT)
shapiro.test(resid(mod_ACT))  
kruskal.test(Value~Pop,data = ACT_chan)


ACT_cssp126<-ACT_chan[ACT_chan$Time=="ssp126_change",]
ACT_cssp585<-ACT_chan[ACT_chan$Time=="ssp585_change",]
kruskal.test(Value~Pop,data = ACT_cssp126)
kruskal.test(Value~Pop,data = ACT_cssp585)




leveneTest(OC_chan$Value~OC_chan$Pop)
mod_OC<-lm(Value~Time+Pop,data = OC_chan)
summary(mod_OC)
shapiro.test(resid(mod_OC))  
kruskal.test(Value~Pop,data = OC_chan)


OC_cssp126<-OC_chan[OC_chan$Time=="ssp126_change",]
OC_cssp585<-OC_chan[OC_chan$Time=="ssp585_change",]
kruskal.test(Value~Pop,data = OC_cssp126)
kruskal.test(Value~Pop,data = OC_cssp585)



leveneTest(WL_chan$Value~WL_chan$Pop)
mod_WL<-lm(Value~Time+Pop,data = WL_chan)
summary(mod_WL)
shapiro.test(resid(mod_WL))  

WL_cssp126<-WL_chan[WL_chan$Time=="ssp126_change",]
WL_cssp585<-WL_chan[WL_chan$Time=="ssp585_change",]
kruskal.test(Value~Pop,data = WL_cssp126)
kruskal.test(Value~Pop,data = WL_cssp585)




#####plot the distribution of sample

limit="3600m"

for (limit in c("3500m","3600m","3700m")) {
  

path_occ=paste("E:/QT/SDM/",limit,"/SpeciesRecord/",sep="")
sample<-data.frame()

for (pop in c("Low","High")) {
  
  occ<-as.data.frame(read.csv(paste(path_occ,"Lizard_",pop,".csv",sep=""),header=T))[,3:4]
  colnames(occ)<-c("Longitude","Latitude")
  occ["Pop"]<-pop
  
  sample<-rbind(sample,occ)
}

elmat = raster("E:/QT/SDM/wc2.1_10m_elev.tif")
elmat=crop(elmat,c(79,105,24,43))
elmat.df<-data.frame(rasterToPoints(elmat))
colnames(elmat.df) <- c("lon", "lat", "ele")
library(RColorBrewer)
library(ggspatial)
library(ggThemeAssist)

colormap<-colorRampPalette(brewer.pal(9,"Greys"))(5)

Figure_sample<-ggplot()+
  geom_raster(data=elmat.df, aes(lon, lat, fill = ele),alpha=0.8)+
  scale_fill_gradientn(values = scales::rescale(c(0,2500,4000,5500,8000)),colors = colormap)+
  geom_point(data=sample,aes(y=Latitude, x=Longitude,colour=Pop),size=1,shape = 16)+

  scale_color_manual(values=c("#4169E1","#DC143C"))+
  #guides(colour=guide_legend(title = NULL))+
  labs(col = "Population")+
  guides(col = guide_legend(override.aes = list(size=4)))+
  labs(fill = "Elevation (m)")+
  theme(legend.text = element_text(size = 18,family = "serif",face = "bold"), 
        legend.background = element_rect(fill = NA))+ 
  theme_void()+
  coord_equal()

ggsave(Figure_sample,file=paste("E:/QT/SDM/",limit,"/SDM_Group/Popdistribution_",limit,".pdf",sep=""),device = "pdf",height=3,width=4,dpi = 300)
}
