#' Global implementation of the microclimate model with custom climate data support
#'
#' An implementation of the NicheMapR microclimate model that allows users to specify
#' custom NetCDF climate data files (e.g., WorldClim v2.1 data). This function is derived
#' from the original micro_global function but adds the ability to use custom climate data
#' paths instead of the default global climate database.
#'
#' The original micro_global function uses the global climate database derived from
#' "New, M., Lister, D., Hulme, M. and Makin, I., 2002: A high-resolution data
#' set of surface climate over global land areas. Climate Research 21:1-25"
#'
#' This modified version allows using WorldClim v2.1 or other custom climate data
#' converted to NetCDF format.
#'
#' @encoding UTF-8
#' @param loc Longitude and latitude (decimal degrees)
#' @param timeinterval The number of time intervals to generate predictions for over a year (must be 12 <= x <=365)
#' @param nyears The number of years to run
#' @param dem A digital elevation model used produced by microclima function 'get_dem' via R package 'elevatr'
#' @param REFL Soil solar reflectance, decimal \%
#' @param elev Elevation, if to be user specified (m)
#' @param slope Slope in degrees
#' @param aspect Aspect in degrees (0 = north)
#' @param DEP Soil depths at which calculations are to be made (cm), must be 10 values starting from 0
#' @param soiltype Soil type: Rock = 0, sand = 1, loamy sand = 2, sandy loam = 3, loam = 4, silt loam = 5, sandy clay loam = 6, clay loam = 7, silt clay loam = 8, sandy clay = 9, silty clay = 10, clay = 11, user-defined = 12
#' @param minshade Minimum shade level to use (\%)
#' @param maxshade Maximum shade level to use (\%)
#' @param Usrhyt Local height (m) at which air temperature, wind speed and humidity are to be computed
#' @param nc_folder Folder path containing custom NetCDF climate data files. If NA, uses default global climate data.
#' @param nc_file Name of the custom NetCDF climate data file (e.g., "worldclim_current1970_2000.nc"). If NA, uses default "global_climate.nc".
#' @param soilmoist_file Name of the custom soil moisture NetCDF file. If NA, uses default "soilw.mon.ltm.v2.nc".
#' @param ... Additional arguments, see Details
#' @usage micro_world(loc = c(-89.40123, 43.07305), timeinterval = 12, nyears = 1, soiltype = 4,
#' REFL = 0.15, slope = 0, aspect = 0,
#' DEP = c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200), minshade = 0, maxshade = 90,
#' Usrhyt = 0.01, nc_folder = NA, nc_file = NA, soilmoist_file = NA, ...)
#' @return metout The above ground micrometeorological conditions under the minimum specified shade
#' @return shadmet The above ground micrometeorological conditions under the maximum specified shade
#' @return soil Hourly predictions of the soil temperatures under the minimum specified shade
#' @return shadsoil Hourly predictions of the soil temperatures under the maximum specified shade
#' @return soilmoist Hourly predictions of the soil moisture under the minimum specified shade
#' @return shadmoist Hourly predictions of the soil moisture under the maximum specified shade
#' @return soilpot Hourly predictions of the soil water potential under the minimum specified shade
#' @return shadpot Hourly predictions of the soil water potential under the maximum specified shade
#' @return humid Hourly predictions of the soil humidity under the minimum specified shade
#' @return shadhumid Hourly predictions of the soil humidity under the maximum specified shade
#' @return plant Hourly predictions of plant transpiration, leaf water potential and root water potential
#' @return shadplant Hourly predictions of plant transpiration, leaf water potential and root water potential
#' @return sunsnow Hourly predictions of snow temperature under the minimum specified shade
#' @return shadsnow Hourly predictions snow temperature under the maximum specified shade
#' @return tcond Hourly predictions of the soil thermal conductivity under the minimum specified shade
#' @return shadtcond Hourly predictions of the soil thermal conductivity under the maximum specified shade
#' @return specheat Hourly predictions of the soil specific heat capacity under the minimum specified shade
#' @return shadspecheat Hourly predictions of soil specific heat capacity under the maximum specified shade
#' @return densit Hourly predictions of the soil density under the minimum specified shade
#' @return shaddensit Hourly predictions of the soil density under the maximum specified shade
#' @details
#' \itemize{
#' \strong{Custom Climate Data Parameters:}\cr\cr
#'
#' \code{nc_folder}{ = NA, Folder path containing custom NetCDF climate data files. If NA, uses default global climate data from NicheMapR.}\cr\cr
#' \code{nc_file}{ = NA, Name of the custom NetCDF climate data file. Examples: "worldclim_current1970_2000.nc", "worldclim_ssp126_2081_2100.nc", "worldclim_ssp585_2081_2100.nc". If NA, uses default "global_climate.nc".}\cr\cr
#' \code{soilmoist_file}{ = NA, Name of the custom soil moisture NetCDF file. If NA, uses default "soilw.mon.ltm.v2.nc".}\cr\cr
#'
#' \strong{Parameters controlling how the model runs:}\cr\cr
#'
#' \code{runshade}{ = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)?}\cr\cr
#' \code{clearsky}{ = 0, Run for clear skies (1) or with observed cloud cover (0)}\cr\cr
#' \code{run.gads}{ = 1, Use the Global Aerosol Database? 1=yes (Fortran version), 2=yes (R version), 0=no}\cr\cr
#' \code{IR}{ = 0, Clear-sky longwave radiation computed using Campbell and Norman (1998) eq. 10.10 (includes humidity) (0) or Swinbank formula (1)}\cr\cr
#' \code{solonly}{ = 0, Only run SOLRAD to get solar radiation? 1=yes, 0=no}\cr\cr
#' \code{lamb}{ = 0, Return wavelength-specific solar radiation output?}\cr\cr
#' \code{IUV}{ = 0, Use gamma function for scattered solar radiation? (computationally intensive)}\cr\cr
#' \code{Soil_Init}{ = NA, initial soil temperature at each soil node, °C (if NA, will use the mean air temperature to initialise)}\cr\cr
#' \code{write_input}{ = 0, Write csv files of final input to folder 'csv input' in working directory? 1=yes, 0=no}\cr\cr
#' \code{writecsv}{ = 0, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{elevatr}{ = 0, Use elevatr package to get high resolution elevation for location? 1 = yes, 0 = no}\cr\cr
#' \code{terrain}{ = 0, Use elevatr package to adjust horizon angles, slope and aspect? 1 = yes, 0 = no}\cr\cr
#' \code{microclima}{ = 0, Use microclima and elevatr package to compute diffuse fraction of solar radiation (1) and adjust solar radiation for terrain (2)? 0 = no}\cr\cr
#' \code{soilgrids}{ = 0, query soilgrids.org database for soil hydraulic properties?}\cr\cr
#' \code{message}{ = 0, allow the Fortran integrator to output warnings? (1) or not (0)}\cr\cr
#' \code{fail}{ = nyears x 24 x 365, how many restarts of the integrator before the Fortran program quits}\cr\cr
#'
#' \strong{ General additional parameters:}\cr\cr
#' \code{ERR}{ = 1.5, Integrator error tolerance for soil temperature calculations}\cr\cr
#' \code{RUF}{ = 0.004, Roughness height (m)}\cr\cr
#' \code{ZH}{ = 0, heat transfer roughness height (m) for Campbell and Norman air temperature/wind speed profile}\cr\cr
#' \code{D0}{ = 0, zero plane displacement correction factor (m)}\cr\cr
#' \code{Z01}{ = 0, Top (1st) segment roughness height(m)}\cr\cr
#' \code{Z02}{ = 0, 2nd segment roughness height(m)}\cr\cr
#' \code{ZH1}{ = 0, Top of (1st) segment, height above surface(m)}\cr\cr
#' \code{ZH2}{ = 0, 2nd segment, height above surface(m)}\cr\cr
#' \code{EC}{ = 0.0167238, Eccentricity of the earth's orbit}\cr\cr
#' \code{SLE}{ = 0.95, Substrate longwave IR emissivity (decimal \%)}\cr\cr
#' \code{Thcond}{ = 2.5, Soil minerals thermal conductivity (W/mK)}\cr\cr
#' \code{Density}{ = 2.56, Soil minerals density (Mg/m3)}\cr\cr
#' \code{SpecHeat}{ = 870, Soil minerals specific heat (J/kg-K)}\cr\cr
#' \code{BulkDensity}{ = 1.3, Soil bulk density (Mg/m3)}\cr\cr
#' \code{PCTWET}{ = 0, \% of ground surface area acting as a free water surface}\cr\cr
#' \code{cap}{ = 1, organic cap present on soil surface?}\cr\cr
#' \code{CMH2O}{ = 1, Precipitable cm H2O in air column}\cr\cr
#' \code{hori}{ = rep(0,24), Horizon angles (degrees)}\cr\cr
#' \code{lapse_min}{ = 0.0039 Lapse rate for minimum air temperature (degrees C/m)}\cr\cr
#' \code{lapse_max}{ = 0.0077 Lapse rate for maximum air temperature (degrees C/m)}\cr\cr
#' \code{TIMAXS}{ = c(1, 1, 0, 0), Time of Maximums for Air Wind RelHum Cloud (h)}\cr\cr
#' \code{TIMINS}{ = c(0, 0, 1, 1), Time of Minimums for Air Wind RelHum Cloud (h)}\cr\cr
#' \code{timezone}{ = 0, Use GNtimezone function to correct to local time zone? 1=yes, 0=no}\cr\cr
#' \code{TAI}{ = 0, Vector of 111 values for solar attenuation}\cr\cr
#' \code{windfac}{ = 1, factor to multiply wind speed by}\cr\cr
#' \code{warm}{ = 0, warming offset vector, °C}\cr\cr
#'
#' \strong{ Soil moisture mode parameters:}\cr\cr
#' \code{runmoist}{ = 0, Run soil moisture model? 1=yes, 0=no}\cr\cr
#' \code{PE}{ = rep(1.1,19), Air entry potential (J/kg)}\cr\cr
#' \code{KS}{ = rep(0.0037,19), Saturated conductivity, (kg s/m^3)}\cr\cr
#' \code{BB}{ = rep(4.5,19), Campbell's soil 'b' parameter (-)}\cr\cr
#' \code{BD}{ = rep(1.3,19), Soil bulk density (Mg/m3)}\cr\cr
#' \code{DD}{ = rep(2.56,19), Soil density (Mg/m3)}\cr\cr
#' \code{maxpool}{ = 10000, Max depth for water pooling on the surface (mm)}\cr\cr
#' \code{rainmult}{ = 1, Rain multiplier for surface soil moisture}\cr\cr
#' \code{evenrain}{ = 0, Spread daily rainfall evenly across 24hrs (1) or one event at midnight (0)}\cr\cr
#' \code{SoilMoist_Init}{ = c(0.1,0.12,0.15,0.2,0.25,0.3,0.3,0.3,0.3,0.3), initial soil water content}\cr\cr
#' \code{L}{ = c(0,0,8.2,8.0,7.8,7.4,7.1,6.4,5.8,4.8,4.0,1.8,0.9,0.6,0.8,0.4,0.4,0,0)*10000, root density}\cr\cr
#' \code{R1}{ = 0.001, root radius, m}\cr\cr
#' \code{RW}{ = 2.5e+10, resistance per unit length of root}\cr\cr
#' \code{RL}{ = 2e+6, resistance per unit length of leaf}\cr\cr
#' \code{PC}{ = -1500, critical leaf water potential for stomatal closure}\cr\cr
#' \code{SP}{ = 10, stability parameter for stomatal closure equation}\cr\cr
#' \code{IM}{ = 1e-06, maximum allowable mass balance error}\cr\cr
#' \code{MAXCOUNT}{ = 500, maximum iterations for mass balance}\cr\cr
#' \code{LAI}{ = 0.1, leaf area index}\cr\cr
#' \code{microclima.LAI}{ = 0, leaf area index for package microclima}\cr\cr
#' \code{microclima.LOR}{ = 1, leaf orientation for package microclima}\cr\cr
#'
#' \strong{ Snow mode parameters:}\cr\cr
#' \code{snowmodel}{ = 0, run the snow model 1=yes, 0=no}\cr\cr
#' \code{snowtemp}{ = 1.5, Temperature (°C) at which precipitation falls as snow}\cr\cr
#' \code{snowdens}{ = 0.375, snow density (Mg/m3)}\cr\cr
#' \code{densfun}{ = c(0.5979, 0.2178, 0.001, 0.0038), snow density function parameters}\cr\cr
#' \code{snowmelt}{ = 1, proportion of calculated snowmelt that doesn't refreeze}\cr\cr
#' \code{undercatch}{ = 1, undercatch multipier for converting rainfall to snow}\cr\cr
#' \code{rainmelt}{ = 0.0125, parameter in equation that melts snow with rainfall}\cr\cr
#' \code{rainfrac}{ = 0.5, fraction of rain that falls on the first day of the month}\cr\cr
#' \code{snowcond}{ = 0, effective snow thermal conductivity W/mC}\cr\cr
#' \code{intercept}{ = max(maxshade) / 100 * 0.3, snow interception fraction}\cr\cr
#' \code{grasshade}{ = 0, if 1, means shade is removed when snow is present}\cr\cr
#' }
#'
#' @examples
#' # Run with default global climate data (original behavior)
#' micro <- micro_world()
#'
#' # Run with custom WorldClim v2.1 climate data
#' micro <- micro_world(
#'   loc = c(100, 35),
#'   nc_folder = "D:/适生区/trae/nc_data",
#'   nc_file = "worldclim_current1970_2000.nc"
#' )
#'
#' # Run with future climate scenario SSP126
#' micro <- micro_world(
#'   loc = c(100, 35),
#'   nc_folder = "D:/适生区/trae/nc_data",
#'   nc_file = "worldclim_ssp126_2081_2100.nc"
#' )
#'
#' # Run with future climate scenario SSP585
#' micro <- micro_world(
#'   loc = c(100, 35),
#'   nc_folder = "D:/适生区/trae/nc_data",
#'   nc_file = "worldclim_ssp585_2081_2100.nc"
#' )
#'
#' @export
micro_world <- function (loc = c(-89.4557, 43.1379), timeinterval = 12, nyears = 1, 
                         soiltype = 4, REFL = 0.15, elev = NA, slope = 0, aspect = 0, 
                         lapse_max = 0.0077, lapse_min = 0.0039, DEP = c(0, 2.5, 
                                                                         5, 10, 15, 20, 30, 50, 100, 200), minshade = 0, maxshade = 90, 
                         dem = NA, Refhyt = 1.2, Usrhyt = 0.01, Z01 = 0, Z02 = 0, 
                         ZH1 = 0, ZH2 = 0, runshade = 1, solonly = 0, clearsky = 0, 
                         run.gads = 1, Soil_Init = NA, write_input = 0, writecsv = 0, 
                         elevatr = 0, terrain = 0, microclima = 0, ERR = 1.5, RUF = 0.004, 
                         ZH = 0, D0 = 0, EC = 0.0167238, SLE = 0.95, Thcond = 2.5, 
                         Density = 2.56, SpecHeat = 870, BulkDensity = 1.3, PCTWET = 0, 
                         cap = 1, CMH2O = 1, hori = rep(0, 24), TIMAXS = c(1, 1, 
                                                                           0, 0), TIMINS = c(0, 0, 1, 1), timezone = 0, runmoist = 0, 
                         PE = rep(1.1, 19), KS = rep(0.0037, 19), BB = rep(4.5, 19), 
                         BD = rep(BulkDensity, 19), DD = rep(Density, 19), maxpool = 10000, 
                         rainmult = 1, evenrain = 0, SoilMoist_Init = c(0.1, 0.12, 
                                                                        0.15, 0.2, 0.25, 0.3, 0.3, 0.3, 0.3, 0.3), L = c(0, 
                                                                                                                         0, 8.2, 8, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4, 1.8, 0.9, 
                                                                                                                         0.6, 0.8, 0.4, 0.4, 0, 0) * 10000, R1 = 0.001, RW = 2.5e+10, 
                         RL = 2e+06, PC = -1500, SP = 10, IM = 1e-06, MAXCOUNT = 500, 
                         LAI = 0.1, microclima.LAI = 0, microclima.LOR = 1, snowmodel = 0, 
                         snowtemp = 1.5, snowdens = 0.375, densfun = c(0.5979, 0.2178, 
                                                                       0.001, 0.0038), snowmelt = 1, undercatch = 1, rainmelt = 0.0125, 
                         rainfrac = 0.5, shore = 0, tides = 0, lamb = 0, IUV = 0, 
                         soilgrids = 0, IR = 0, message = 0, fail = nyears * 24 * 
                           365, TAI = 0, warm = 0, windfac = 1, snowcond = 0, intercept = max(maxshade)/100 * 
                           0.3, grasshade = 0, maxsurf = 95,
                         nc_folder = NA, nc_file = NA, soilmoist_file = NA) 
{
  SoilMoist <- SoilMoist_Init
  errors <- 0
  if (DEP[2] - DEP[1] > 3 | DEP[3] - DEP[2] > 3) {
    message("warning, nodes might be too far apart near the surface, try a different spacing if the program is crashing \n")
  }
  if (DEP[2] - DEP[1] < 2) {
    cat("warning, nodes might be too close near the surface, try a different spacing if the program is crashing \n")
  }
  if (DEP[10] != 200) {
    cat("warning, last depth in soil should not be changed from 200 without good reason \n")
  }
  if (timeinterval < 12 | timeinterval > 365) {
    message("ERROR: the variable 'timeinterval' is out of bounds.\n        Please enter a correct value (12 - 365).", 
            "\n")
    errors <- 1
  }
  if (is.numeric(loc[1])) {
    if (loc[1] > 180 | loc[2] > 90) {
      message("ERROR: Latitude or longitude (longlat) is out of bounds.\n        Please enter a correct value.", 
              "\n")
      errors <- 1
    }
  }
  if (timezone %in% c(0, 1) == FALSE) {
    message("ERROR: the variable 'timezone' be either 0 or 1.\n      Please correct.", 
            "\n")
    errors <- 1
  }
  if (run.gads == 1) {
    message("If program is crashing, try run.gads = 2.", 
            "\n")
  }
  if (run.gads %in% c(0, 1, 2) == FALSE) {
    message("ERROR: the variable 'run.gads' be either 0, 1 or 2.\n      Please correct.", 
            "\n")
    errors <- 1
  }
  if (write_input %in% c(0, 1) == FALSE) {
    message("ERROR: the variable 'write_input' be either 0 or 1.\n      Please correct.", 
            "\n")
    errors <- 1
  }
  if (EC < 0.0034 | EC > 0.058) {
    message("ERROR: the eccentricity variable (EC) is out of bounds.\n        Please enter a correct value (0.0034 - 0.058).", 
            "\n")
    errors <- 1
  }
  if (RUF < 1e-04) {
    message("ERROR: The roughness height (RUF) is too small ( < 0.0001).\n        Please enter a larger value.", 
            "\n")
    errors <- 1
  }
  if (RUF > 2) {
    message("ERROR: The roughness height (RUF) is too large ( > 2).\n        Please enter a smaller value.", 
            "\n")
    errors <- 1
  }
  if (D0 > 0 & D0 < Usrhyt) {
    cat("ERROR: The zero plane displacement height (D0) must be lower than the local height (Usrhyt).\n        Please enter a smaller value.", 
        "\n")
    errors <- 1
  }
  if (DEP[1] != 0) {
    message("ERROR: First soil node (DEP[1]) must = 0 cm.\n        Please correct", 
            "\n")
    errors <- 1
  }
  if (length(DEP) != 10) {
    message("ERROR: You must enter 10 different soil depths.", 
            "\n")
    errors <- 1
  }
  for (i in 1:9) {
    if (DEP[i + 1] <= DEP[i]) {
      message("ERROR: Soil depth (DEP array) is not in ascending size", 
              "\n")
      errors <- 1
    }
  }
  if (DEP[10] > 500) {
    message("ERROR: Deepest soil depth (DEP array) is too large (<=500 cm)", 
            "\n")
    errors <- 1
  }
  if (min(Thcond) < 0) {
    cat("ERROR: Thermal variable conductivity (THCOND) is negative.\n        Please input a positive value.", 
        "\n")
    errors <- 1
  }
  if (min(Density) < 0) {
    cat("ERROR: Density variable (Density) is negative.\n        Please input a positive value.", 
        "\n")
    errors <- 1
  }
  if (min(SpecHeat) < 0) {
    cat("ERROR: Specific heat variable (SpecHeat) is negative.\n        Please input a positive value.", 
        "\n")
    errors <- 1
  }
  if (min(BulkDensity) < 0) {
    message("ERROR: Bulk density value (BulkDensity) is negative.\n        Please input a positive value.", 
            "\n")
    errors <- 1
  }
  if (REFL < 0 | REFL > 1) {
    message("ERROR: Soil reflectivity value (REFL) is out of bounds.\n        Please input a value between 0 and 1.", 
            "\n")
    errors <- 1
  }
  if (slope < 0 | slope > 90) {
    message("ERROR: Slope value (slope) is out of bounds.\n        Please input a value between 0 and 90.", 
            "\n")
    errors <- 1
  }
  if (aspect < 0 | aspect > 365) {
    message("ERROR: Aspect value (aspect) is out of bounds.\n        Please input a value between 0 and 365.", 
            "\n")
    errors <- 1
  }
  if (max(hori) > 90 | min(hori) < 0) {
    message("ERROR: At least one of your horizon angles (hori) is out of bounds.\n        Please input a value between 0 and 90", 
            "\n")
    errors <- 1
  }
  if (length(hori) != 24) {
    message("ERROR: You must enter 24 horizon angle values.", 
            "\n")
    errors <- 1
  }
  if (SLE < 0.05 | SLE > 1) {
    message("ERROR: Emissivity (SLE) is out of bounds.\n        Please enter a correct value (0.05 - 1.00).", 
            "\n")
    errors <- 1
  }
  if (ERR < 0) {
    message("ERROR: Error bound (ERR) is too small.\n        Please enter a correct value (> 0.00).", 
            "\n")
    errors <- 1
  }
  if (Usrhyt < RUF) {
    message("ERROR: Reference height (Usrhyt) smaller than roughness height (RUF).\n        Please use a larger height above the surface.", 
            "\n")
    errors <- 1
  }
  if (Usrhyt > Refhyt) {
    message("ERROR: Reference height is less than local height (Usrhyt) \n")
    errors <- 1
  }
  if (CMH2O < 0.5 | CMH2O > 2) {
    message("ERROR: Preciptable water in air column (CMH2O) is out of bounds.\n        Please enter a correct value (0.1 - 2).", 
            "\n")
    errors <- 1
  }
  if (max(TIMAXS) > 24 | min(TIMAXS) < 0) {
    message("ERROR: At least one of your times of weather maxima (TIMAXS) is out of bounds.\n        Please input a value between 0 and 24", 
            "\n")
    errors <- 1
  }
  if (max(TIMINS) > 24 | min(TIMINS) < 0) {
    message("ERROR: At least one of your times of weather minima (TIMINS) is out of bounds.\n        Please input a value between 0 and 24", 
            "\n")
    errors <- 1
  }
  if (max(minshade - maxshade) >= 0) {
    cat("ERROR: Your value(s) for minimum shade (minshade) is greater than or equal to the maximum shade (maxshade).\n        Please correct this.", 
        "\n")
    errors <- 1
  }
  if (max(minshade) > 100 | min(minshade) < 0) {
    cat("ERROR: Your value(s) for minimum shade (minshade) is out of bounds.\n        Please input a value between 0 and 100.", 
        "\n")
    errors <- 1
  }
  if (max(maxshade) > 100 | min(maxshade) < 0) {
    cat("ERROR: Your value(s) for maximum shade (maxshade) is out of bounds.\n        Please input a value between 0 and 100.", 
        "\n")
    errors <- 1
  }
  if (soiltype < 0 | soiltype > 11) {
    message("ERROR: the soil type must range between 1 and 11.\n      Please correct.", 
            "\n")
    errors <- 1
  }
  if (errors == 0) {
    doys12 <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 
                288, 319, 349)
    doysn <- doys12
    if (nyears > 1 & timeinterval == 365) {
      for (i in 1:(nyears - 1)) {
        doysn <- c(doysn, (doys12 + 365 * i))
      }
    }
    if (timeinterval < 365) {
      microdaily <- 0
    }
    else {
      microdaily <- 1
    }
    daystart <- as.integer(ceiling(365/timeinterval/2))
    if (timeinterval != 12) {
      doys <- seq(daystart, 365, as.integer(floor(365/timeinterval)))
    }
    else {
      doys <- doysn
    }
    doynum <- timeinterval * nyears
    doy <- subset(doys, doys != 0)
    doy <- rep(doy, nyears)
    idayst <- 1
    ida <- timeinterval * nyears
    if (is.numeric(loc) == FALSE) {
      loc = cbind(as.numeric(loc[1]), as.numeric(loc[2]))
    }
    longlat <- loc
    x <- t(as.matrix(as.numeric(c(loc[1], loc[2]))))
    if (timezone == 1) {
      if (!require(geonames, quietly = TRUE)) {
        stop("package \"geonames\" is required to do a specific time zone (timezone=1). Please install it.")
      }
      ALREF <- (geonames::GNtimezone(longlat[2], longlat[1])[4]) * 
        -15
    }
    else {
      ALREF <- abs(trunc(x[1]))
    }
    HEMIS <- ifelse(x[2] < 0, 2, 1)
    ALAT <- abs(trunc(x[2]))
    AMINUT <- (abs(x[2]) - ALAT) * 60
    ALONG <- abs(trunc(x[1]))
    ALMINT <- (abs(x[1]) - ALONG) * 60
    azmuth <- aspect
    if (terrain == 1) {
      elevatr <- 1
    }
    if (is.na(elev) & elevatr == 1) {
      require(microclima)
      require(terra)
      cat("downloading DEM via package elevatr \n")
      dem <- microclima::get_dem(lat = loc[2], long = loc[1])
      dem_terra <- terra::rast(dem)
      xy = sf::st_as_sf(data.frame(lon = loc[1], lat = loc[2]), 
                        coords = c("lon", "lat"))
      xy <- sf::st_set_crs(xy, "EPSG:4326")
      xy <- sf::st_transform(xy, sf::st_crs(dem_terra))
      elev <- as.numeric(terra::extract(dem_terra, xy)[, 
                                                       2])
      if (terrain == 1) {
        cat("computing slope, aspect and horizon angles \n")
        slope <- terra::terrain(dem_terra, v = "slope", 
                                unit = "degrees")
        slope <- as.numeric(terra::extract(slope, xy)[, 
                                                      2])
        aspect <- terra::terrain(dem_terra, v = "aspect", 
                                 unit = "degrees")
        aspect <- as.numeric(terra::extract(aspect, 
                                            xy[, 2]))
        ha24 <- 0
        for (i in 0:23) {
          har <- microclima::horizonangle(dem, i * 10, 
                                          terra::res(dem)[1])
          ha24[i + 1] <- atan(as.numeric(terra::extract(har, 
                                                        xy)[, 2])) * (180/pi)
        }
        hori <- ha24
      }
    }
    hori <- as.matrix(hori)
    VIEWF <- 1 - sum(sin(as.data.frame(hori) * pi/180))/length(hori)
    SLES <- rep(SLE, timeinterval * nyears)
    if (length(minshade) != timeinterval * nyears) {
      MINSHADES <- rep(0, (timeinterval * nyears)) + minshade[1]
    }
    else {
      MINSHADES <- rep(0, (timeinterval * nyears)) + minshade
    }
    if (length(maxshade) != timeinterval * nyears) {
      MAXSHADES <- rep(0, (timeinterval * nyears)) + maxshade[1]
    }
    else {
      MINSHADES <- rep(0, (timeinterval * nyears)) + minshade
    }
    if (soiltype == 0) {
      BulkDensity <- Density
      cap = 0
      runmoist <- 0
      PE <- rep(CampNormTbl9_1[1, 4], 19)
      KS <- rep(CampNormTbl9_1[1, 6], 19)
      BB <- rep(CampNormTbl9_1[1, 5], 19)
      BD <- rep(BulkDensity, 19)
      DD <- rep(Density, 19)
    }
    else {
      if (soiltype < 12) {
        PE <- rep(CampNormTbl9_1[soiltype, 4], 19)
        KS <- rep(CampNormTbl9_1[soiltype, 6], 19)
        BB <- rep(CampNormTbl9_1[soiltype, 5], 19)
        BD <- rep(BulkDensity, 19)
        DD <- rep(Density, 19)
      }
    }
    if (soilgrids == 1) {
      cat("extracting data from SoilGrids \n")
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("package 'json' is needed to extract data from SoilGrids, please install it.", 
             call. = FALSE)
      }
      require(jsonlite)
      ov <- jsonlite::fromJSON(paste0("https://rest.isric.org/soilgrids/v2.0/properties/query?lon=", 
                                      x[1], "&lat=", x[2], "&property=bdod&property=silt&property=clay&property=sand"), 
                               flatten = TRUE)
      if (length(ov) > 3) {
        soilpro <- cbind(c(0, 5, 15, 30, 60, 100), unlist(ov$properties$layers$depths[[1]]$values.mean)/100, 
                         unlist(ov$properties$layers$depths[[2]]$values.mean)/10, 
                         unlist(ov$properties$layers$depths[[4]]$values.mean)/10, 
                         unlist(ov$properties$layers$depths[[3]]$values.mean)/10)
        soilpro <- rbind(soilpro, soilpro[6, ])
        soilpro[7, 1] <- 200
        colnames(soilpro) <- c("depth", "blkdens", "clay", 
                               "silt", "sand")
        soil.hydro <- pedotransfer(soilpro = as.data.frame(soilpro), 
                                   DEP = DEP)
        PE <- soil.hydro$PE
        BB <- soil.hydro$BB
        BD <- soil.hydro$BD
        KS <- soil.hydro$KS
        BulkDensity <- BD[seq(1, 19, 2)]
      }
      else {
        cat("no SoilGrids data for this site, using user-input soil properties \n")
      }
    }
    gcfolder <- paste(.libPaths()[1], "/gcfolder.rda", sep = "")
    
    if (!is.na(nc_folder) && !is.na(nc_file)) {
      folder <- gsub("/$", "", nc_folder)
      nc_path <- file.path(folder, nc_file)
      if (!file.exists(nc_path)) {
        message(paste0("指定的nc文件不存在: ", nc_path, "\n请检查路径是否正确。\n退出函数 micro_world"))
        opt <- options(show.error.messages = FALSE)
        on.exit(options(opt))
        stop()
      }
      message(paste0("使用自定义气候数据: ", nc_path, "\n"))
    } else if (file.exists(gcfolder) == FALSE) {
      folder <- "c:/globalclimate"
      if (file.exists(paste0(folder, "/global_climate.nc")) == 
          FALSE) {
        message("You don't appear to have the global climate data set - \n run function get.global.climate(folder = 'folder you want to put it in') .....\n exiting function micro_world")
        opt <- options(show.error.messages = FALSE)
        on.exit(options(opt))
        stop()
      }
    }
    else {
      load(gcfolder)
    }
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("package 'terra' is needed. Please install it.", 
           call. = FALSE)
    }
    if (!requireNamespace("ncdf4", quietly = TRUE)) {
      stop("package 'ncdf4' is needed. Please install it.", 
           call. = FALSE)
    }
    message("extracting climate data \n")
    if (!is.na(nc_folder) && !is.na(nc_file)) {
      global_climate <- terra::rast(file.path(folder, nc_file))
    } else {
      global_climate <- terra::rast(paste0(folder, "/global_climate.nc"))
    }
    CLIMATE <- t(as.numeric(terra::extract(global_climate, 
                                           x)))
    ALTT <- as.numeric(CLIMATE[1])
    delta_elev <- 0
    if (is.na(elev) == FALSE) {
      delta_elev <- ALTT - elev
      ALTT <- elev
    }
    adiab_corr_max <- delta_elev * lapse_max
    adiab_corr_min <- delta_elev * lapse_min
    RAINFALL <- CLIMATE[2:13]
    if (is.na(RAINFALL[1])) {
      cat("no climate data for this site, using dummy data so solar is still produced \n")
      CLIMATE <- t(as.numeric(terra::extract(global_climate, 
                                             cbind(140, -35))))
      CLIMATE[2:97] <- 0
      ALTT <- as.numeric(CLIMATE[1])
      delta_elev <- 0
      if (is.na(elev) == FALSE) {
        delta_elev <- ALTT - elev
        ALTT <- elev
      }
      adiab_corr_max <- delta_elev * lapse_max
      adiab_corr_min <- delta_elev * lapse_min
      RAINFALL <- CLIMATE[2:13] * 0
    }
    RAINYDAYS <- CLIMATE[14:25]/10
    WNMAXX <- CLIMATE[26:37]/10 * windfac
    WNMINN <- WNMAXX * 0.1
    TMINN <- CLIMATE[38:49]/10
    TMAXX <- CLIMATE[50:61]/10
    TMAXX <- TMAXX + adiab_corr_max
    TMINN <- TMINN + adiab_corr_min
    ALLMINTEMPS <- TMINN
    ALLMAXTEMPS <- TMAXX
    ALLTEMPS <- cbind(ALLMAXTEMPS, ALLMINTEMPS)
    RHMINN <- CLIMATE[62:73]/10
    RHMAXX <- CLIMATE[74:85]/10
    es <- WETAIR(db = TMAXX, rh = 100)$esat
    e <- WETAIR(db = CLIMATE[50:61]/10, rh = CLIMATE[62:73]/10)$e
    RHMINN <- (e/es) * 100
    RHMINN[RHMINN > 100] <- 100
    RHMINN[RHMINN < 0] <- 0.01
    es <- WETAIR(db = TMINN, rh = 100)$esat
    e <- WETAIR(db = CLIMATE[38:49]/10, rh = CLIMATE[74:85]/10)$e
    RHMAXX <- (e/es) * 100
    RHMAXX[RHMAXX > 100] <- 100
    RHMAXX[RHMAXX < 0] <- 0.01
    CCMINN <- CLIMATE[86:97]/10
    if (clearsky == 1) {
      CCMINN <- CCMINN * 0
    }
    CCMAXX <- CCMINN
    if (runmoist == 0) {
      if (!is.na(nc_folder) && !is.na(soilmoist_file)) {
        soilmoisture <- suppressWarnings(terra::rast(file.path(folder, soilmoist_file)))
      } else {
        soilmoisture <- suppressWarnings(terra::rast(paste(folder, 
                                                           "/soilw.mon.ltm.v2.nc", sep = "")))
      }
      message("extracting soil moisture data")
      SoilMoist <- t(as.numeric(terra::extract(soilmoisture, 
                                               x)))/1000
    }
    if (is.na(max(SoilMoist, ALTT, CLIMATE)) == TRUE) {
      message("Sorry, there is no environmental data for this location")
      SoilMoist <- t(as.numeric(terra::extract(soilmoisture, 
                                               cbind(140, -35))))/1000
    }
    WNMINN <- WNMINN * (1.2/10)^0.15
    WNMAXX <- WNMAXX * (1.2/10)^0.15
    TMAXX <- TMAXX + warm
    TMINN <- TMINN + warm
    if (timeinterval != 12) {
      TMAXX1 <- suppressWarnings(spline(doys12, TMAXX, 
                                        n = timeinterval, xmin = 1, xmax = 365, method = "periodic"))
      TMAXX <- rep(TMAXX1$y, nyears)
      TMINN1 <- suppressWarnings(spline(doys12, TMINN, 
                                        n = timeinterval, xmin = 1, xmax = 365, method = "periodic"))
      TMINN <- rep(TMINN1$y, nyears)
      RHMAXX1 <- suppressWarnings(spline(doys12, RHMAXX, 
                                         n = timeinterval, xmin = 1, xmax = 365, method = "periodic"))
      RHMAXX <- rep(RHMAXX1$y, nyears)
      RHMINN1 <- suppressWarnings(spline(doys12, RHMINN, 
                                         n = timeinterval, xmin = 1, xmax = 365, method = "periodic"))
      RHMINN <- rep(RHMINN1$y, nyears)
      CCMAXX1 <- suppressWarnings(spline(doys12, CCMAXX, 
                                         n = timeinterval, xmin = 1, xmax = 365, method = "periodic"))
      CCMAXX <- rep(CCMAXX1$y, nyears)
      CCMINN <- CCMAXX
      WNMAXX1 <- suppressWarnings(spline(doys12, WNMAXX, 
                                         n = timeinterval, xmin = 1, xmax = 365, method = "periodic"))
      WNMAXX <- rep(WNMAXX1$y, nyears)
      WNMINN1 <- suppressWarnings(spline(doys12, WNMINN, 
                                         n = timeinterval, xmin = 1, xmax = 365, method = "periodic"))
      WNMINN <- rep(WNMINN1$y, nyears)
      if (runmoist == 0) {
        SoilMoist1 <- suppressWarnings(spline(doys12, 
                                              SoilMoist, n = timeinterval, xmin = 1, xmax = 365, 
                                              method = "periodic"))
        SoilMoist <- rep(SoilMoist1$y, nyears)
      }
    }
    if (timeinterval < 365) {
      TMAXX <- rep(TMAXX, nyears)
      TMINN <- rep(TMINN, nyears)
      RHMAXX <- rep(RHMAXX, nyears)
      RHMINN <- rep(RHMINN, nyears)
      CCMAXX <- rep(CCMAXX, nyears)
      CCMINN <- rep(CCMINN, nyears)
      WNMAXX <- rep(WNMAXX, nyears)
      WNMINN <- rep(WNMINN, nyears)
      if (runmoist == 0) {
        SoilMoist <- rep(SoilMoist, nyears)
      }
      RAINFALL <- rep(RAINFALL, nyears)
    }
    orig.RAINFALL <- RAINFALL
    avetemp <- (sum(TMAXX) + sum(TMINN))/(length(TMAXX) * 
                                            2)
    if (is.na(Soil_Init[1])) {
      soilinit <- rep(avetemp, 20)
      spinup <- 1
    }
    else {
      if (snowmodel == 0) {
        soilinit <- c(Soil_Init, rep(avetemp, 10))
      }
      else {
        soilinit <- c(rep(avetemp, 8), Soil_Init[1:10], 
                      rep(avetemp, 2))
      }
      spinup <- 0
    }
    tannul <- mean(unlist(ALLTEMPS))
    tannulrun <- rep(tannul, doynum)
    daymon <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 
                30, 31)
    if (timeinterval == 365) {
      RAINFALL1 <- 1:365
      sort <- matrix(data = 0, nrow = 365, ncol = 2)
      m <- 1
      b <- 0
      for (i in 1:12) {
        ndays = daymon[i]
        for (k in 1:ndays) {
          b <- b + 1
          sort[m, 1] <- i
          sort[m, 2] <- b
          if (k <= RAINYDAYS[i] & rainfrac > 0) {
            if (k == 1) {
              RAINFALL1[m] <- RAINFALL[i] * rainfrac * 
                rainmult
            }
            else {
              RAINFALL1[m] <- (RAINFALL[i] * (1 - rainfrac) * 
                                 rainmult)/RAINYDAYS[i]
            }
          }
          else {
            if (rainfrac == 0) {
              RAINFALL1[m] <- (RAINFALL[i] * rainmult)/RAINYDAYS[i]
            }
            else {
              RAINFALL1[m] <- 0
            }
          }
          m <- m + 1
          if (b > RAINYDAYS[i]) {
            b <- 0
          }
        }
      }
      RAINFALL2 <- as.data.frame(cbind(RAINFALL1, sort))
      RAINFALL <- rep(as.double(RAINFALL2$RAINFALL1), 
                      nyears)
      RAINFALL[!is.finite(RAINFALL)] <- 0
      if (TMINN[1] < snowtemp) {
        RAINFALL[1] <- 0
      }
    }
    else {
      if (timeinterval != 12) {
        RAINFALL <- rep(rep(sum(RAINFALL)/timeinterval, 
                            timeinterval), nyears)
      }
      else {
        RAINFALL <- RAINFALL/rep(daymon, nyears)
      }
    }
    ndays <- length(RAINFALL)
    SOLRhr <- rep(0, 24 * ndays)
    hourly <- 0
    if (microclima > 0 & timeinterval %in% c(12, 365)) {
      cat("using microclima and elevatr to adjust solar for topographic and vegetation effects \n")
      if (!require("microclima", quietly = TRUE)) {
        stop("package 'microclima' is needed. Please install it.", 
             call. = FALSE)
      }
      if (!require("zoo", quietly = TRUE)) {
        stop("package 'zoo' is needed. Please install it.", 
             call. = FALSE)
      }
      cat("Downloading digital elevation data \n")
      lat <- x[2]
      long <- x[1]
      yearlist <- seq(1960, 1960 + nyears - 1)
      tt <- seq(as.POSIXct(paste0("01/01/", yearlist[1]), 
                           format = "%d/%m/%Y", tz = "UTC"), as.POSIXct(paste0("31/12/", 
                                                                               yearlist[nyears]), format = "%d/%m/%Y", tz = "UTC") + 
                  23 * 3600, by = "hours")
      timediff <- x[1]/15
      hour.microclima <- as.numeric(format(tt, "%H")) + 
        timediff - floor(timediff)
      jd <- julday(as.numeric(format(tt, "%Y")), as.numeric(format(tt, 
                                                                   "%m")), as.numeric(format(tt, "%d")))
      if (!is.raster(dem)) {
        dem <- microclima::get_dem(r = NA, lat = lat, 
                                   long = long, resolution = 100, zmin = -20)
      }
      require(terra)
      dem_terra <- terra::rast(dem)
      xy = sf::st_as_sf(data.frame(lon = loc[1], lat = loc[2]), 
                        coords = c("lon", "lat"))
      xy <- sf::st_set_crs(xy, "EPSG:4326")
      xy <- sf::st_transform(xy, sf::st_crs(dem_terra))
      if (class(slope) == "logical") {
        slope <- terra::terrain(dem, v = "slope", unit = "degrees")
        slope <- as.numeric(terra::extract(slope, xy))
      }
      if (class(aspect) == "logical") {
        aspect <- terrain(dem, v = "aspect", unit = "degrees")
        aspect <- as.numeric(terra::extract(aspect, 
                                            xy))
      }
      ha <- 0
      ha36 <- 0
      for (i in 0:35) {
        har <- horizonangle(dem, i * 10, res(dem)[1])
        ha36[i + 1] <- atan(as.numeric(terra::extract(har, 
                                                      xy))) * (180/pi)
      }
      for (i in 1:length(hour.microclima)) {
        saz <- solazi(hour.microclima[i], lat, long, 
                      jd[i], merid = long)
        saz <- round(saz/10, 0) + 1
        saz <- ifelse(saz > 36, 1, saz)
        ha[i] <- ha36[saz]
      }
      cloud <- rep(CCMAXX/2, nyears)
      methspline <- "periodic"
      for (i in 1:nyears) {
        if (yearlist[i] %in% seq(1900, 2060, 4)) {
          xmax <- 366
        }
        else {
          xmax <- 365
        }
        if (i == 1) {
          start <- 1
          end <- 12
          cloud1 <- suppressWarnings(spline(doys12, 
                                            cloud[start:end], n = xmax, xmin = 1, xmax = xmax, 
                                            method = methspline))
          cloud2 <- cloud1$y
        }
        else {
          start <- end + 1
          end <- end + 12
          cloud1 <- suppressWarnings(spline(c(0, doys12), 
                                            c(tail(cloud2, 1), cloud[start:end]), n = xmax, 
                                            xmin = 1, xmax = xmax, method = methspline))
          cloud2 <- c(cloud2, cloud1$y)
        }
      }
      cloudhr <- cbind(rep(seq(1, length(cloud2)), 24), 
                       rep(cloud2, 24))
      cloudhr <- cloudhr[order(cloudhr[, 1]), ]
      cloudhr <- cloudhr[, 2]
      micro_clearsky <- micro_world(loc = c(x[1], x[2]), 
                                    clearsky = 1, TAI = TAI, timeinterval = 365, 
                                    solonly = 1)
      clearskyrad <- micro_clearsky$metout[, c(1, 13)][, 
                                                       2]
      dsw2 <- leapfix(clearskyrad, yearlist, 24) * (0.36 + 
                                                      0.64 * (1 - cloudhr/100))
      si <- microclima::siflat(hour.microclima, lat, long, 
                               jd, merid = long)
      am <- microclima::airmasscoef(hour.microclima, lat, 
                                    long, jd, merid = long)
      dp <- vector(length = length(jd))
      for (i in 1:length(jd)) {
        dp[i] <- microclima:::difprop(dsw2[i], jd[i], 
                                      hour.microclima[i], lat, long, watts = TRUE, 
                                      hourly = TRUE, merid = long)
      }
      dp[dsw2 == 0] <- NA
      dnir <- (dsw2 * (1 - dp))/si
      dnir[si == 0] <- NA
      difr <- (dsw2 * dp)
      edni <- dnir/((4.87/0.0036) * (1 - dp))
      edif <- difr/((4.87/0.0036) * dp)
      bound <- function(x, mn = 0, mx = 1) {
        x[x > mx] <- mx
        x[x < mn] <- mn
        x
      }
      odni <- bound((log(edni)/-am), mn = 0.001, mx = 1.7)
      odif <- bound((log(edif)/-am), mn = 0.001, mx = 1.7)
      nd <- length(odni)
      sel <- which(is.na(am * dp * odni * odif) == F)
      dp[1] <- dp[min(sel)]
      odni[1] <- odni[min(sel)]
      odif[1] <- odif[min(sel)]
      dp[nd] <- dp[max(sel)]
      odni[nd] <- odni[max(sel)]
      odif[nd] <- odif[max(sel)]
      dp[nd] <- dp[max(sel)]
      odni[nd] <- odni[max(sel)]
      odif[nd] <- odif[max(sel)]
      if (!require("terra", quietly = TRUE)) {
        stop("package 'terra' is needed. Please install it.", 
             call. = FALSE)
      }
      dp <- na.approx(dp, na.rm = F)
      odni <- na.approx(odni, na.rm = F)
      odif <- na.approx(odif, na.rm = F)
      h_dp <- bound(dp)
      h_oi <- bound(odni, mn = 0.24, mx = 1.7)
      h_od <- bound(odif, mn = 0.24, mx = 1.7)
      afi <- exp(-am * h_oi)
      afd <- exp(-am * h_od)
      h_dni <- (1 - h_dp) * afi * 4.87/0.0036
      h_dif <- h_dp * afd * 4.87/0.0036
      h_dni[si == 0] <- 0
      h_dif[is.na(h_dif)] <- 0
      diffuse_frac_all <- h_dif/(h_dni + h_dif)
      diffuse_frac_all[is.na(diffuse_frac_all)] <- 1
      radwind2 <- .shortwave.ts(h_dni * 0.0036, h_dif * 
                                  0.0036, jd, hour.microclima, lat, long, slope, 
                                aspect, ha = ha, svv = 1, x = microclima.LOR, 
                                l = mean(microclima.LAI), albr = 0, merid = long, 
                                dst = 0, difani = FALSE)
      SOLRhr_all <- radwind2$swrad/0.0036
      diffuse_frac <- diffuse_frac_all
      if (microclima == 2) {
        hourly <- 2
        VIEWF <- 1
        hori <- rep(0, 24)
      }
    }
    else {
      diffuse_frac <- NA
    }
    if (timeinterval == 12 & microclima > 0) {
      dates_all <- head(seq(as.POSIXct(paste0("01/01/", 
                                              yearlist[1]), format = "%d/%m/%Y", tz = "UTC"), 
                            as.POSIXct(paste0("01/01/", yearlist[nyears] + 
                                                1), format = "%d/%m/%Y ", tz = "UTC"), by = "hours"), 
                        -1)
      dates_15th <- which(format(dates_all, "%d") == "15")
      diffuse_frac <- diffuse_frac_all[dates_15th]
    }
    if (length(TAI) < 111) {
      if (run.gads > 0) {
        relhum <- 1
        if (run.gads == 1) {
          optdep.summer <- as.data.frame(rungads(longlat[2], 
                                                 longlat[1], relhum, 0))
          optdep.winter <- as.data.frame(rungads(longlat[2], 
                                                 longlat[1], relhum, 1))
        }
        else {
          optdep.summer <- as.data.frame(gads.r(longlat[2], 
                                                longlat[1], relhum, 0))
          optdep.winter <- as.data.frame(gads.r(longlat[2], 
                                                longlat[1], relhum, 1))
        }
        optdep <- cbind(optdep.winter[, 1], rowMeans(cbind(optdep.summer[, 
                                                                         2], optdep.winter[, 2])))
        optdep <- as.data.frame(optdep)
        colnames(optdep) <- c("LAMBDA", "OPTDEPTH")
        a <- lm(OPTDEPTH ~ poly(LAMBDA, 6, raw = TRUE), 
                data = optdep)
        LAMBDA <- c(290, 295, 300, 305, 310, 315, 320, 
                    330, 340, 350, 360, 370, 380, 390, 400, 420, 
                    440, 460, 480, 500, 520, 540, 560, 580, 600, 
                    620, 640, 660, 680, 700, 720, 740, 760, 780, 
                    800, 820, 840, 860, 880, 900, 920, 940, 960, 
                    980, 1000, 1020, 1080, 1100, 1120, 1140, 1160, 
                    1180, 1200, 1220, 1240, 1260, 1280, 1300, 
                    1320, 1380, 1400, 1420, 1440, 1460, 1480, 
                    1500, 1540, 1580, 1600, 1620, 1640, 1660, 
                    1700, 1720, 1780, 1800, 1860, 1900, 1950, 
                    2000, 2020, 2050, 2100, 2120, 2150, 2200, 
                    2260, 2300, 2320, 2350, 2380, 2400, 2420, 
                    2450, 2490, 2500, 2600, 2700, 2800, 2900, 
                    3000, 3100, 3200, 3300, 3400, 3500, 3600, 
                    3700, 3800, 3900, 4000)
        TAI <- predict(a, data.frame(LAMBDA))
      }
      else {
        TAI <- c(0.42, 0.415, 0.412, 0.408, 0.404, 0.4, 
                 0.395, 0.388, 0.379, 0.379, 0.379, 0.375, 
                 0.365, 0.345, 0.314, 0.3, 0.288, 0.28, 0.273, 
                 0.264, 0.258, 0.253, 0.248, 0.243, 0.236, 
                 0.232, 0.227, 0.223, 0.217, 0.213, 0.21, 0.208, 
                 0.205, 0.202, 0.201, 0.198, 0.195, 0.193, 
                 0.191, 0.19, 0.188, 0.186, 0.184, 0.183, 0.182, 
                 0.181, 0.178, 0.177, 0.176, 0.175, 0.175, 
                 0.174, 0.173, 0.172, 0.171, 0.17, 0.169, 0.168, 
                 0.167, 0.164, 0.163, 0.163, 0.162, 0.161, 
                 0.161, 0.16, 0.159, 0.157, 0.156, 0.156, 0.155, 
                 0.154, 0.153, 0.152, 0.15, 0.149, 0.146, 0.145, 
                 0.142, 0.14, 0.139, 0.137, 0.135, 0.135, 0.133, 
                 0.132, 0.131, 0.13, 0.13, 0.129, 0.129, 0.128, 
                 0.128, 0.128, 0.127, 0.127, 0.126, 0.125, 
                 0.124, 0.123, 0.121, 0.118, 0.117, 0.115, 
                 0.113, 0.11, 0.108, 0.107, 0.105, 0.103, 0.1)
      }
    }
    Nodes <- matrix(data = 0, nrow = 10, ncol = ndays)
    if (soilgrids == 1) {
      Numtyps <- 10
      Nodes[1:10, ] <- c(1:10)
    }
    else {
      Numtyps <- 2
      Nodes[1, 1:ndays] <- 3
      Nodes[2, 1:ndays] <- 9
    }
    REFLS <- rep(REFL, ndays)
    PCTWET <- rep(PCTWET, ndays)
    if (runmoist == 0) {
      moists2 <- matrix(nrow = 10, ncol = ndays, data = 0)
      moists2[1, ] <- SoilMoist
      moists2[2, ] <- moists2[1, ]
      moists <- moists2
    }
    else {
      moists2 <- matrix(nrow = 10, ncol = ndays, data = 0)
      moists2[1:10, ] <- SoilMoist_Init
      moists2[moists2 > (1 - BulkDensity/Density)] <- (1 - 
                                                         BulkDensity/Density)
      moists <- moists2
    }
    soilprops <- matrix(data = 0, nrow = 10, ncol = 5)
    if (soilgrids == 1) {
      soilprops[, 1] <- BulkDensity
      soilprops[, 2] <- 1 - BulkDensity/Density
      soilprops[soilprops[, 2] < 0.26, 2] <- 0.26
      soilprops[, 3] <- Thcond
      soilprops[, 4] <- SpecHeat
      soilprops[, 5] <- Density
      if (cap == 1) {
        soilprops[1:2, 3] <- 0.2
        soilprops[1:2, 4] <- 1920
      }
      if (cap == 2) {
        soilprops[1:2, 3] <- 0.1
        soilprops[3:4, 3] <- 0.25
        soilprops[1:4, 4] <- 1920
        soilprops[1:4, 5] <- 1.3
        soilprops[1:4, 1] <- 0.7
      }
    }
    else {
      soilprops[1, 1] <- BulkDensity
      soilprops[2, 1] <- BulkDensity
      soilprops[1, 2] <- min(0.26, 1 - BulkDensity/Density)
      soilprops[2, 2] <- min(0.26, 1 - BulkDensity/Density)
      if (cap == 1) {
        soilprops[1, 3] <- 0.2
      }
      else {
        soilprops[1, 3] <- Thcond
      }
      soilprops[2, 3] <- Thcond
      if (cap == 1) {
        soilprops[1, 4] <- 1920
      }
      else {
        soilprops[1, 4] <- SpecHeat
      }
      soilprops[2, 4] <- SpecHeat
      soilprops[1, 5] <- Density
      soilprops[2, 5] <- Density
    }
    hourly <- 0
    rainhourly <- 0
    TAIRhr <- rep(0, 24 * ndays)
    RHhr <- rep(0, 24 * ndays)
    WNhr <- rep(0, 24 * ndays)
    CLDhr <- rep(0, 24 * ndays)
    SOLRhr <- rep(0, 24 * ndays)
    RAINhr <- rep(0, 24 * ndays)
    ZENhr <- rep(-1, 24 * ndays)
    IRDhr <- rep(-1, 24 * ndays)
    microinput <- c(ndays, RUF, ERR, Usrhyt, Refhyt, Numtyps, 
                    Z01, Z02, ZH1, ZH2, idayst, ida, HEMIS, ALAT, AMINUT, 
                    ALONG, ALMINT, ALREF, slope, azmuth, ALTT, CMH2O, 
                    microdaily, tannul, EC, VIEWF, snowtemp, snowdens, 
                    snowmelt, undercatch, rainmult, runshade, runmoist, 
                    maxpool, evenrain, snowmodel, rainmelt, writecsv, 
                    densfun, hourly, rainhourly, lamb, IUV, RW, PC, 
                    RL, SP, R1, IM, MAXCOUNT, IR, message, fail, snowcond, 
                    intercept, grasshade, solonly, ZH, D0, TIMAXS, TIMINS, 
                    spinup, 0, 360, maxsurf)
    if (length(LAI) < ndays) {
      LAI <- rep(LAI[1], ndays)
    }
    if (shore == 0) {
      tides <- matrix(data = 0, nrow = 24 * ndays, ncol = 3)
    }
    micro <- list(tides = tides, microinput = microinput, 
                  doy = doy, SLES = SLES, DEP = DEP, Nodes = Nodes, 
                  MAXSHADES = MAXSHADES, MINSHADES = MINSHADES, TMAXX = TMAXX, 
                  TMINN = TMINN, RHMAXX = RHMAXX, RHMINN = RHMINN, 
                  CCMAXX = CCMAXX, CCMINN = CCMINN, WNMAXX = WNMAXX, 
                  WNMINN = WNMINN, TAIRhr = TAIRhr, RHhr = RHhr, WNhr = WNhr, 
                  CLDhr = CLDhr, SOLRhr = SOLRhr, RAINhr = RAINhr, 
                  ZENhr = ZENhr, IRDhr = IRDhr, REFLS = REFLS, PCTWET = PCTWET, 
                  soilinit = soilinit, hori = hori, TAI = TAI, soilprops = soilprops, 
                  moists = moists, RAINFALL = RAINFALL, tannulrun = tannulrun, 
                  PE = PE, KS = KS, BB = BB, BD = BD, DD = DD, L = L, 
                  LAI = LAI)
    if (write_input == 1) {
      if (dir.exists("micro csv input") == FALSE) {
        dir.create("micro csv input")
      }
      write.table(as.matrix(microinput), file = "micro csv input/microinput.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(doy, file = "micro csv input/doy.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(SLES, file = "micro csv input/SLES.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(DEP, file = "micro csv input/DEP.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(Nodes, file = "micro csv input/Nodes.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(MAXSHADES, file = "micro csv input/Maxshades.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(MINSHADES, file = "micro csv input/Minshades.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(TIMAXS, file = "micro csv input/TIMAXS.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(TIMINS, file = "micro csv input/TIMINS.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(TMAXX, file = "micro csv input/TMAXX.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(TMINN, file = "micro csv input/TMINN.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(RHMAXX, file = "micro csv input/RHMAXX.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(RHMINN, file = "micro csv input/RHMINN.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(CCMAXX, file = "micro csv input/CCMAXX.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(CCMINN, file = "micro csv input/CCMINN.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(WNMAXX, file = "micro csv input/WNMAXX.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(WNMINN, file = "micro csv input/WNMINN.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(REFLS, file = "micro csv input/REFLS.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(PCTWET, file = "micro csv input/PCTWET.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(soilinit, file = "micro csv input/soilinit.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(hori, file = "micro csv input/hori.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(TAI, file = "micro csv input/TAI.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(soilprops, file = "micro csv input/soilprop.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(moists, file = "micro csv input/moists.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(RAINFALL, file = "micro csv input/rain.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(tannulrun, file = "micro csv input/tannulrun.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(PE, file = "micro csv input/PE.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(BD, file = "micro csv input/BD.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(DD, file = "micro csv input/DD.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(BB, file = "micro csv input/BB.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(KS, file = "micro csv input/KS.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(L, file = "micro csv input/L.csv", sep = ",", 
                  col.names = NA, qmethod = "double")
      write.table(LAI, file = "micro csv input/LAI.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(tides, file = "micro csv input/tides.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(TAIRhr, file = "micro csv input/TAIRhr.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(RHhr, file = "micro csv input/RHhr.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(WNhr, file = "micro csv input/WNhr.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(CLDhr, file = "micro csv input/CLDhr.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(SOLRhr, file = "micro csv input/SOLRhr.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(RAINhr, file = "micro csv input/RAINhr.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(ZENhr, file = "micro csv input/ZENhr.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(IRDhr, file = "micro csv input/IRDhr.csv", 
                  sep = ",", col.names = NA, qmethod = "double")
    }
    if (is.numeric(loc[1])) {
      location <- paste("long", loc[1], "lat", loc[2])
    }
    else {
      location <- loc
    }
    message(paste("running microclimate model for", timeinterval, 
                  "days by", nyears, "years at site", location, "\n"))
    message("Note: the output column `SOLR` in metout and shadmet is for unshaded horizontal plane solar radiation \n")
    ptm <- proc.time()
    microut <- microclimate(micro)
    message(paste0("runtime ", (proc.time() - ptm)[3], " seconds"))
    metout <- microut$metout
    shadmet <- microut$shadmet
    soil <- microut$soil
    shadsoil <- microut$shadsoil
    tcond <- microut$tcond
    shadtcond <- microut$shadtcond
    specheat <- microut$specheat
    shadspecheat <- microut$shadspecheat
    densit <- microut$densit
    shaddensit <- microut$shaddensit
    if (runmoist == 1) {
      soilmoist <- microut$soilmoist
      shadmoist <- microut$shadmoist
      humid <- microut$humid
      shadhumid <- microut$shadhumid
      soilpot <- microut$soilpot
      shadpot <- microut$shadpot
      plant <- microut$plant
      shadplant <- microut$shadplant
    }
    else {
      soilpot <- soil
      soilmoist <- soil
      shadpot <- soil
      shadmoist <- soil
      humid <- soil
      shadhumid <- soil
      plant <- cbind(soil, soil[, 3:4])
      shadplant <- cbind(soil, soil[, 3:4])
      soilpot[, 3:12] <- 0
      soilmoist[, 3:12] <- 0.5
      shadpot[, 3:12] <- 0
      shadmoist[, 3:12] <- 0.5
      humid[, 3:12] <- 0.99
      shadhumid[, 3:12] <- 0.99
      plant[, 3:14] <- 0
      shadplant[, 3:14] <- 0
    }
    if (snowmodel == 1) {
      sunsnow <- microut$sunsnow
      shdsnow <- microut$shdsnow
    }
    if (timeinterval == 12) {
      RAINFALL <- orig.RAINFALL
    }
    if (max(metout[, 1] == 0)) {
      cat("ERROR: the model crashed - try a different error tolerance (ERR) or a different spacing in DEP", 
          "\n")
    }
    days <- rep(seq(1, timeinterval * nyears), 24)
    days <- days[order(days)]
    dates <- days + metout[, 2]/60/24 - 1
    dates2 <- seq(1, timeinterval * nyears)
    if (lamb == 1) {
      drlam <- as.data.frame(microut$drlam)
      drrlam <- as.data.frame(microut$drrlam)
      srlam <- as.data.frame(microut$srlam)
      if (snowmodel == 1) {
        return(list(soil = soil, shadsoil = shadsoil, 
                    metout = metout, shadmet = shadmet, soilmoist = soilmoist, 
                    shadmoist = shadmoist, humid = humid, shadhumid = shadhumid, 
                    soilpot = soilpot, shadpot = shadpot, sunsnow = sunsnow, 
                    shdsnow = shdsnow, plant = plant, shadplant = shadplant, 
                    tcond = tcond, shadtcond = shadtcond, specheat = specheat, 
                    shadspecheat = shadspecheat, densit = densit, 
                    shaddensit = shaddensit, RAINFALL = RAINFALL, 
                    ndays = ndays, elev = ALTT, REFL = REFL[1], 
                    longlat = c(x[1], x[2]), nyears = nyears, 
                    timeinterval = timeinterval, minshade = MINSHADES, 
                    maxshade = MAXSHADES, DEP = DEP, drlam = drlam, 
                    drrlam = drrlam, srlam = srlam, dates = dates, 
                    dates2 = dates2, PE = PE, BD = BD, DD = DD, 
                    BB = BB, KS = KS, dem = dem, diffuse_frac = diffuse_frac))
      }
      else {
        return(list(soil = soil, shadsoil = shadsoil, 
                    metout = metout, shadmet = shadmet, soilmoist = soilmoist, 
                    shadmoist = shadmoist, humid = humid, shadhumid = shadhumid, 
                    soilpot = soilpot, shadpot = shadpot, plant = plant, 
                    shadplant = shadplant, tcond = tcond, shadtcond = shadtcond, 
                    specheat = specheat, shadspecheat = shadspecheat, 
                    densit = densit, shaddensit = shaddensit, 
                    RAINFALL = RAINFALL, ndays = ndays, elev = ALTT, 
                    REFL = REFL[1], longlat = c(x[1], x[2]), nyears = nyears, 
                    timeinterval = timeinterval, minshade = MINSHADES, 
                    maxshade = MAXSHADES, DEP = DEP, drlam = drlam, 
                    drrlam = drrlam, srlam = srlam, dates = dates, 
                    dates2 = dates2, PE = PE, BD = BD, DD = DD, 
                    BB = BB, KS = KS, dem = dem, diffuse_frac = diffuse_frac))
      }
    }
    else {
      if (snowmodel == 1) {
        return(list(soil = soil, shadsoil = shadsoil, 
                    metout = metout, shadmet = shadmet, soilmoist = soilmoist, 
                    shadmoist = shadmoist, humid = humid, shadhumid = shadhumid, 
                    soilpot = soilpot, shadpot = shadpot, sunsnow = sunsnow, 
                    shdsnow = shdsnow, plant = plant, shadplant = shadplant, 
                    tcond = tcond, shadtcond = shadtcond, specheat = specheat, 
                    shadspecheat = shadspecheat, densit = densit, 
                    shaddensit = shaddensit, RAINFALL = RAINFALL, 
                    ndays = ndays, elev = ALTT, REFL = REFL[1], 
                    longlat = c(x[1], x[2]), nyears = nyears, 
                    timeinterval = timeinterval, minshade = MINSHADES, 
                    maxshade = MAXSHADES, DEP = DEP, dates = dates, 
                    dates2 = dates2, PE = PE, BD = BD, DD = DD, 
                    BB = BB, KS = KS, dem = dem, diffuse_frac = diffuse_frac))
      }
      else {
        return(list(soil = soil, shadsoil = shadsoil, 
                    metout = metout, shadmet = shadmet, soilmoist = soilmoist, 
                    shadmoist = shadmoist, humid = humid, shadhumid = shadhumid, 
                    soilpot = soilpot, shadpot = shadpot, plant = plant, 
                    shadplant = shadplant, tcond = tcond, shadtcond = shadtcond, 
                    specheat = specheat, shadspecheat = shadspecheat, 
                    densit = densit, shaddensit = shaddensit, 
                    RAINFALL = RAINFALL, ndays = ndays, elev = ALTT, 
                    REFL = REFL[1], longlat = c(x[1], x[2]), nyears = nyears, 
                    timeinterval = timeinterval, minshade = MINSHADES, 
                    maxshade = MAXSHADES, DEP = DEP, dates = dates, 
                    dates2 = dates2, PE = PE, BD = BD, DD = DD, 
                    BB = BB, KS = KS, dem = dem, diffuse_frac = diffuse_frac))
      }
    }
  }
}
