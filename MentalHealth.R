# Setup -------------------------------------------------------------------
rm(list=ls()) # Clears workspace

#system("sudo apt install libgeos-dev libproj-dev libgdal-dev libudunits2-dev -y") # Install linux geospatial dependencies 

# Install/call libraries
install.packages("renv")
#renv::init()

PKG <- c("googledrive","tidyr","purrr", "sf", "tmap", "raster", "rnaturalearth", "rgdal", "exactextractr","ggplot2")

for (p in PKG) {
  if(!require(p,character.only = TRUE)) {  
    install.packages(p)
    require(p,character.only = TRUE)}
}

renv::snapshot()
rm(p,PKG)

## Load data
# Download from google drive to directory "Data"
setwd("~/Github/GZ")
dir.create(file.path('Data'), recursive = TRUE)
folder_url<-"https://drive.google.com/open?id=1yMjM2V09xyrQ2Oekv93IX8YsrNQjBDWg"
folder<-drive_get(as_id(folder_url))
files<-drive_ls(folder)
dl<-function(files){
  walk(files, ~ drive_download(as_id(.x), overwrite = TRUE))
}
setwd("./Data")
system.time(map(files$id,dl))
setwd("..")
rm(files, folder, folder_url, dl)

# Load data to R
NDVI<-raster("./Data/NDVI_mean_Guangzhou.tif")
NDVI<-NDVI/10000 # Rescaling
Pop<-raster("./Data/population.tif") # Doesn't appear to have coverage across the whole serviceshed
Pop2<-raster("./Data/chn_ppp_2020_guangzhou.tif")
hw<-st_read("./Data/haizhu_wetland.gpkg")
aoi<-st_read("./Data/aoi.gpkg")
lulc<-raster("./Data/lulc.tif")
codes<-read.csv("./Data/lulc_codes.csv")
NDVIs<-raster("./Data/ndvi_residential.tif")
NDVIs<-NDVIs/10000 # Rescaling

## Preparing data
prj<-CRS("+init=epsg:32649") # WGS84/UTM49N

# Reprojecting
Pop<-projectRaster(Pop, crs=prj) 
Pop2<-projectRaster(Pop2, crs=prj) 
NDVI<-projectRaster(NDVI, crs=prj)
NDVIs<-projectRaster(NDVIs, crs=prj)

# Removing negative NDVI values (as in Liu et al., 2019)
NDVI[NDVI < 0] <- NA
NDVIs[NDVIs < 0] <- NA

# Buffering the polygon by 1 and 2km - changes in the wetland can affect neighborhoods up to 1km outside of the wetland, and neighborhoods themselves are affected by NDVI within a 1km buffer
hw_1km<-st_buffer(hw, dist = 1000) 
hw_2km<-st_buffer(hw, dist = 2000) 

# Creating a layer that is only the 2km wetland buffer
hw_buff<-st_difference(hw_2km,hw)

## Analysis
# 1. Raster to point for each pop cell inside or within 1 km of wetland
rwl_1km<-mask(Pop2,hw_1km) # mask pop raster to serviceshed - note that "crop" function is different
rwl_p<-rasterToPoints(rwl_1km,spatial = TRUE) # raster to points
rwl_p<-st_as_sf(rwl_p) # SP object to SF (faster/modern spatial format)
# 2. Buffer each point by 1km
rwl_b<-st_buffer(rwl_p,dist=1000)
# 4. Extract values of NDVI for old NDVI raster
rwl_b$NDVI_0<-exact_extract(NDVI,rwl_b)
# Calculate mean 
# https://stackoverflow.com/questions/45317327/apply-function-to-columns-in-a-list-of-data-frames-and-append-results
# https://stackoverflow.com/questions/57834130/pmap-purrr-error-argument-1-must-have-names
rwl_b$meanNDVI_0<-map_dbl(rwl_b$NDVI_0, function(x){
  a<-mean(x$value,na.rm = TRUE)
  return(a)
})
# Calculate count of pixels sampled
rwl_b$countNDVI_0<-map_dbl(rwl_b$NDVI_0, function(x){
  a<-length(na.omit(x$value)) # only counts non-NA cells
  return(a)
})
# 5. Extract values of NDVI for new scenario, using provided "wallpaper" map
rwl_b$NDVI_1<-exact_extract(NDVIs,rwl_b)
# Calculate mean 
# https://stackoverflow.com/questions/45317327/apply-function-to-columns-in-a-list-of-data-frames-and-append-results
# https://stackoverflow.com/questions/57834130/pmap-purrr-error-argument-1-must-have-names
rwl_b$meanNDVI_1<-map_dbl(rwl_b$NDVI_1, function(x){
  a<-mean(x$value,na.rm = TRUE)
  return(a)
})
# Calculate count of pixels sampled
rwl_b$countNDVI_1<-map_dbl(rwl_b$NDVI_1, function(x){
  a<-length(na.omit(x$value)) # only counts non-NA cells
  return(a)
})
# # 5. Extract all NDVI values inside 2km buffer, i.e. NDVI conditions outside the wetland boundary, up to the furthest reach of neighborhoods that could be impacted by a change in wetland conditions
# NDVI_2km<-exact_extract(NDVI,hw_buff)
# NDVI_2km<-NDVI_2km[[1]]$value
# # 6. Extract values of NDVI for new scenario, using values from buffer around wetland
# # Create mask of NDVI to only 2km buffer
# NDVI_2km_crop<-mask(NDVI,hw_buff)
# # Extracting NDVI from the portion of neighborhood outside the wetland
# rwl_b$NDVI_1<-exact_extract(NDVI_2km_crop,rwl_b)
# # Calculate mean for NDVI outside the wetland
# rwl_b$meanNDVI_1<-map_dbl(rwl_b$NDVI_1, function(x){
#   a<-mean(x$value,na.rm = TRUE)
#   return(a)
# })
# # Calculate count of pixels sampled
# rwl_b$countNDVI_1<-map_dbl(rwl_b$NDVI_1, function(x){
#   a<-length(na.omit(x$value)) # only counts non-NA cells
#   return(a)
# })
# # Some NaN for NDVI extract from outside the wetland, for neighborhoods that lie entirely within the wetland, replace with zero
# rwl_b$meanNDVI_1[is.nan(rwl_b$meanNDVI_1)]<-0
# # 7. Weighted mean NDVI from inside/outside wetland boundary for new scenario
# set.seed(5)
# rwl_b$meanNDVI_1<-(mean(sample(NDVI_2km,rwl_b$countNDVI_0-rwl_b$countNDVI_1,replace = FALSE),na.rm = TRUE))*((rwl_b$countNDVI_0-rwl_b$countNDVI_1)/rwl_b$countNDVI_0) + # Portion of mean from inside
#   (1-((rwl_b$countNDVI_0-rwl_b$countNDVI_1)/rwl_b$countNDVI_0))*rwl_b$meanNDVI_1 # Portion of mean from outside wetland    
# 8. Subset to only variables of interest
rwl_b<-rwl_b %>% dplyr::select(chn_ppp_2020_guangzhou, meanNDVI_0, meanNDVI_1, countNDVI_0, countNDVI_1)
# 9. Calculate % change in mental health (WHO-5 score) using function from Liu et al. (2019)
rwl_b$PTcWHO_5<-(rwl_b$meanNDVI_1-rwl_b$meanNDVI_0)/0.1356 # Point change
rwl_b$PCTcWHO_5<-rwl_b$PTcWHO_5/12.081 # Percent change
# 10. Aggregate and calculate net difference
rwl_b$PPcValue<-rwl_b$PCTcWHO_5*356.74 # Per capita change in value from Xu et al. (2016)
rwl_b$cValue<-rwl_b$PPcValue*rwl_b$chn_ppp_2020_guangzhou # Total change in value per cell
sumcValue<-sum(rwl_b$cValue) # Total change in value
# 11. Net present value
# https://stackoverflow.com/questions/47403180/calculate-npv-for-cashflows-at-all-point-in-time
dcf <- function(x, r, t0=FALSE){
  # calculates discounted cash flows (DCF) given cash flow and discount rate
  #
  # x - cash flows vector
  # r - vector or discount rates, in decimals. Single values will be recycled
  # t0 - cash flow starts in year 0, default is FALSE, i.e. discount rate in first period is zero.
  if(length(r)==1){
    r <- rep(r, length(x))
    if(t0==TRUE){r[1]<-0}
  }
  x/cumprod(1+r)
}

npv <- function(x, r, t0=FALSE){
  # calculates net present value (NPV) given cash flow and discount rate
  #
  # x - cash flows vector
  # r - discount rate, in decimals
  # t0 - cash flow starts in year 0, default is FALSE
  sum(dcf(x, r, t0))
}
x<-rep(sumcValue,30) # vector of annual value, for a given number of years
NPVcValue<-npv(x,0.05) # NPV, for a given discount rate and number of years
# 12. Other results
# Population affected
sumpop<-round(sum(rwl_b$chn_ppp_2020_guangzhou))
# Histogram of neighborhood value change
ggplot(rwl_b, aes(x=cValue)) + geom_histogram(color="black", fill="white") +
  theme_minimal() +
  labs(title="Change in value by neighborhood",x="Value (2015 USD)", y = "Count")
# Mean NDVI by scenario
mean(rwl_b$meanNDVI_0)
mean(rwl_b$meanNDVI_1)

## Maps
cval<-rwl_b %>% dplyr::select(cValue,geometry) 
cvalc<-st_centroid(cval)

ggplot() +
  geom_sf(data = hw_2km, fill=NA, color = "black", lty = "dashed") + 
  geom_sf(data = cvalc, aes(color=cValue)) +
  geom_sf(data = hw, fill=NA, color = "black") + 
  theme_void() + 
  labs(colour = "Value Change per \nYear (2020 USD)\n") 

p1<-tm_shape(hw_2km) +
  tm_borders(col = NA) +
  tm_shape(NDVI) +
  tm_raster(palette = "Greens", colorNA = NULL, title = "NDVI") + 
  tm_shape(hw_2km) +
  tm_borders(col = "black", lty = "dashed") +
  tm_shape(hw) +
  tm_borders(col = "black") +  
  tm_layout(legend.outside = TRUE)

p2<-tm_shape(hw_2km) +
  tm_borders(col = NA) +
  tm_shape(NDVIs) +
  tm_raster(palette = "Greens", colorNA = NULL, title = "NDVI") + 
  tm_shape(hw_2km) +
  tm_borders(col = "black", lty = "dashed") +
  tm_shape(hw) +
  tm_borders(col = "black") +  
  tm_layout(legend.outside = TRUE)

tmap_arrange(p1,p2, nrow = 2)  

#tm_raster(style = "fixed", breaks=seq(0, 5, by=1), labels = c("<1","1 - 2","2 - 3","3 - 4",">4"), palette = cols, title = "Water Depth (m)",colorNA = "white",textNA = "")

