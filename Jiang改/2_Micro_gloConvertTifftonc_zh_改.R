######################将TIFF堆栈转换为nc文件######################
# 适配中国境内研究，更新路径和11个GCM模型
# 核心功能：整合多源气候数据，计算多模式集合平均，输出为NetCDF格式

# -------------------------- 环境初始化 --------------------------
rm(list = ls())
library(raster)   
library(ncdf4)    

# -------------------------- 基础气候数据加载 --------------------------
# 读取基础全球气候NetCDF文件
global_climate_tif <- stack("D:/适生区/气候数据/NicheMapR/Data/global_climate.nc")
Current <- global_climate_tif              

# 从基础气候数据中提取各要素（按图层索引）
Altitude <- Current[[1]]                   # 第1层：海拔数据
RAINYDAYS <- Current[[14:25]]              # 第14-25层：降雨天数（12个月）
WNMAXX <- Current[[26:37]]                 # 第26-37层：最大风速（12个月）
WNMINN <- WNMAXX * 0.1                     # 计算最小风速（假设为最大风速的10%）
RHMINN <- Current[[62:73]]                 # 第62-73层：最小相对湿度（12个月）
RHMAXX <- Current[[74:85]]                 # 第74-85层：最大相对湿度（12个月）
CCMINN <- Current[[86:97]]                 # 第86-97层：最小云量（12个月）

# -------------------------- 1970-2000 当前气候数据处理 --------------------------
# 1. 读取月降雨量数据（12个波段在一个文件中）
path_current_pre <- "D:/适生区/气候数据/worldclim/1970-2000/wc2.1_10m_prec/"
current_pre_file <- list.files(path_current_pre, pattern = "*.tif", full.names = TRUE)
current_pre <- stack(current_pre_file)  # 直接读取为12层的堆栈

# 2. 读取月最低温度数据（12个波段在一个文件中）
path_current_Tmin <- "D:/适生区/气候数据/worldclim/1970-2000/wc2.1_10m_tmin/"
current_Tmin_file <- list.files(path_current_Tmin, pattern = "*.tif", full.names = TRUE)
current_Tmin <- stack(current_Tmin_file)
current_Tmin <- current_Tmin * 10  # 温度值缩放（×10）

# 3. 读取月最高温度数据（12个波段在一个文件中）
path_current_Tmax <- "D:/适生区/气候数据/worldclim/1970-2000/wc2.1_10m_tmax/"
current_Tmax_file <- list.files(path_current_Tmax, pattern = "*.tif", full.names = TRUE)
current_Tmax <- stack(current_Tmax_file)
current_Tmax <- current_Tmax * 10  # 温度值缩放（×10）

# 4. 赋值降雨量、温度数据到全局变量
RAINFALL <- current_pre    
TMINN <- current_Tmin      
TMAXX <- current_Tmax      
ALLMINTEMPS <- TMINN       
ALLMAXTEMPS <- TMAXX       
ALLTEMPS <- cbind(ALLMAXTEMPS, ALLMINTEMPS)  

# 5. 整合所有气候要素并输出NC文件
# 注意：原始global_climate.nc有97层，不包含WNMINN
# WNMINN在micro_global函数中通过 WNMINN <- WNMAXX * 0.1 计算
Current1970_2000 <- stack(Altitude, RAINFALL, RAINYDAYS, WNMAXX, TMINN, TMAXX, RHMINN, RHMAXX, CCMINN)
setwd("D:/worldclimout/")
writeRaster(
  Current1970_2000,                       
  "worldclim_current1970_2000.nc",                   
  overwrite = TRUE,                        
  format = "CDF",                          
  varname = "Temperature",                 
  varunit = "degC",                        
  longname = "Temperature -- raster stack to netCDF", 
  xname = "Longitude",                     
  yname = "Latitude",                      
  zname = "Time (Month)"                   
)

# -------------------------- 定义11个GCM模型列表 --------------------------
GCM_list <- c(
  "ACCESS-CM2", "BCC-CSM2-MR", "CMCC-ESM2", 
  "EC-Earth3-Veg", "GISS-E2-1-G", "INM-CM5-0", 
  "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", 
  "MRI-ESM2-0", "UKESM1-0-LL"
)

# -------------------------- 2081-2100 SSP126 未来气候数据处理 --------------------------
# 1. 处理最低温度（SSP126）
path_tn126 <- "D:/适生区/气候数据/worldclim/2081-2100/ssp126/"
tn126 <- stack()  
for (i in 1:12) {  
  stack_mon <- stack()  
  for (GCMs in GCM_list) {
    # 构建文件名，匹配你的文件格式
    tn126file <- raster(paste0(path_tn126, "wc2.1_10m_tmin_", GCMs, "_ssp126_2081-2100.tif"), band = i)
    stack_mon <- stack(stack_mon, tn126file)  
  }
  stack_mean <- mean(stack_mon)  # 计算当月11个GCM的平均值
  tn126 <- stack(tn126, stack_mean)  
}

# 2. 处理最高温度（SSP126）
path_tx126 <- "D:/适生区/气候数据/worldclim/2081-2100/ssp126/"
tx126 <- stack()
for (i in 1:12) {
  stack_mon <- stack()
  for (GCMs in GCM_list) {
    tx126file <- raster(paste0(path_tx126, "wc2.1_10m_tmax_", GCMs, "_ssp126_2081-2100.tif"), band = i)
    stack_mon <- stack(stack_mon, tx126file)
  }
  stack_mean <- mean(stack_mon)
  tx126 <- stack(tx126, stack_mean)
}

# 3. 处理降雨量（SSP126）
path_pr126 <- "D:/适生区/气候数据/worldclim/2081-2100/ssp126/"
pr126 <- stack()
for (i in 1:12) {
  stack_mon <- stack()
  for (GCMs in GCM_list) {
    pr126file <- raster(paste0(path_pr126, "wc2.1_10m_prec_", GCMs, "_ssp126_2081-2100.tif"), band = i)
    stack_mon <- stack(stack_mon, pr126file)
  }
  stack_mean <- mean(stack_mon)
  pr126 <- stack(pr126, stack_mean)
}

# 4. 温度值缩放 + 赋值
tn126 <- tn126 * 10
tx126 <- tx126 * 10
RAINFALL <- pr126       
TMINN <- tn126          
TMAXX <- tx126          
ALLMINTEMPS <- TMINN
ALLMAXTEMPS <- TMAXX
ALLTEMPS <- cbind(ALLMAXTEMPS, ALLMINTEMPS)

# 5. 整合输出NC文件
ssp126_2081_2100 <- stack(Altitude, RAINFALL, RAINYDAYS, WNMAXX, TMINN, TMAXX, RHMINN, RHMAXX, CCMINN)
setwd("D:/worldclimout/")
writeRaster(
  ssp126_2081_2100,
  "worldclim_ssp126_2081_2100.nc",
  overwrite = TRUE,
  format = "CDF",
  varname = "Temperature",
  varunit = "degC",
  longname = "Temperature -- raster stack to netCDF",
  xname = "Longitude",
  yname = "Latitude",
  zname = "Time (Month)"
)

# -------------------------- 2081-2100 SSP585 未来气候数据处理 --------------------------
# 1. 处理最低温度（SSP585）
path_tn585 <- "D:/适生区/气候数据/worldclim/2081-2100/ssp585/"
tn585 <- stack()
for (i in 1:12) {
  stack_mon <- stack()
  for (GCMs in GCM_list) {
    tn585file <- raster(paste0(path_tn585, "wc2.1_10m_tmin_", GCMs, "_ssp585_2081-2100.tif"), band = i)
    stack_mon <- stack(stack_mon, tn585file)
  }
  stack_mean <- mean(stack_mon)
  tn585 <- stack(tn585, stack_mean)
}

# 2. 处理最高温度（SSP585）
path_tx585 <- "D:/适生区/气候数据/worldclim/2081-2100/ssp585/"
tx585 <- stack()
for (i in 1:12) {
  stack_mon <- stack()
  for (GCMs in GCM_list) {
    tx585file <- raster(paste0(path_tx585, "wc2.1_10m_tmax_", GCMs, "_ssp585_2081-2100.tif"), band = i)
    stack_mon <- stack(stack_mon, tx585file)
  }
  stack_mean <- mean(stack_mon)
  tx585 <- stack(tx585, stack_mean)
}

# 3. 处理降雨量（SSP585）
path_pr585 <- "D:/适生区/气候数据/worldclim/2081-2100/ssp585/"
pr585 <- stack()
for (i in 1:12) {
  stack_mon <- stack()
  for (GCMs in GCM_list) {
    pr585file <- raster(paste0(path_pr585, "wc2.1_10m_prec_", GCMs, "_ssp585_2081-2100.tif"), band = i)
    stack_mon <- stack(stack_mon, pr585file)
  }
  stack_mean <- mean(stack_mon)
  pr585 <- stack(pr585, stack_mean)
}

# 4. 温度值缩放 + 赋值
tn585 <- tn585 * 10
tx585 <- tx585 * 10
RAINFALL <- pr585       
TMINN <- tn585          
TMAXX <- tx585          
ALLMINTEMPS <- TMINN
ALLMAXTEMPS <- TMAXX
ALLTEMPS <- cbind(ALLMAXTEMPS, ALLMINTEMPS)

# 5. 整合输出NC文件
ssp585_2081_2100 <- stack(Altitude, RAINFALL, RAINYDAYS, WNMAXX, TMINN, TMAXX, RHMINN, RHMAXX, CCMINN)
setwd("D:/worldclimout/")
writeRaster(
  ssp585_2081_2100,
  "worldclim_ssp585_2081_2100.nc",
  overwrite = TRUE,
  format = "CDF",
  varname = "Temperature",
  varunit = "degC",
  longname = "Temperature -- raster stack to netCDF",
  xname = "Longitude",
  yname = "Latitude",
  zname = "Time (Month)"
)