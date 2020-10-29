# Setup -------------------------------------------------------------------
rm(list=ls()) # Clears workspace

#system("sudo apt install libgeos-dev libproj-dev libgdal-dev libudunits2-dev -y") # Install linux geospatial dependencies 

# Install/call libraries
install.packages("renv")
renv::init()

PKG <- c("googledrive","tidyr","purrr", "sf", "tmap", "raster", "rnaturalearth", "rgdal")

for (p in PKG) {
  if(!require(p,character.only = TRUE)) {  
    install.packages(p)
    require(p,character.only = TRUE)}
}

renv::snapshot()
rm(p,PKG)

## Load data
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
rm(files, folder, folder_url)

NDVI<-raster("./Data/NDVI_mean_Guangzhou.tif")
Pop<-raster("./Data/chn_ppp_2020_guangzhou.tif")

## Visualize data
cs<-ne_download(scale = 10, type = "countries", returnclass = "sf")
r<-ne_download(scale = 10, type = "rivers_lake_centerlines", category = "physical",  returnclass = "sf")
l<-ne_download(scale = 10, type = "land", category = "physical",  returnclass = "sf")

## Sample wetland polygon 
wl<-st_sf(st_as_sfc(st_bbox(Pop)))

# Function to identify UTM zone of coordinates
lonlat2UTM = function(lonlat) {
  utm = (floor((lonlat[1] + 180) / 6) %% 60) + 1
  if(lonlat[2] > 0) {
    utm + 32600
  } else{
    utm + 32700
  }
}
st_coordinates(wl)

EPSG_2_UTM <- lonlat2UTM(c(112.9246,22.62458))
wl<-st_transform(wl, EPSG_2_UTM)
wl<-st_buffer(wl, dist = -55000)
wl_1km<-st_buffer(wl, dist = 1000)

p1<-tm_shape(l, bbox = st_bbox(Pop), ext = .85) +
  tm_polygons(border.col = "black", col = "gray") +
  tm_shape(Pop) +
  tm_raster() +
  tm_shape(cs) +
  tm_borders(col = "black") +
  tm_shape(r) +
  tm_lines(col = "blue") +
  tm_shape(wl) +
  tm_borders(col = "red") +
  tm_shape(wl_1km) +
  tm_borders(col = "black", lty = "dashed")

p2<-tm_shape(l, bbox = st_bbox(Pop), ext = .85) +
  tm_polygons(border.col = "black", col = "gray") +
  tm_shape(NDVI) +
  tm_raster(midpoint = NA) +
  tm_shape(cs) +
  tm_borders(col = "black") +
  tm_shape(r) +
  tm_lines(col = "blue") +
  tm_shape(wl) +
  tm_borders(col = "red") +
  tm_shape(wl_1km) +
  tm_borders(col = "black", lty = "dashed")

tmap_arrange(p1,p2, nrow = 1)  

## Analysis
# 1. Raster to point for each pop cell within 1 km of wetland
# 2. Buffer each point by 1km
# 3. Extract values of NDVI for old NDVI raster/new NDVI raster
# 3a. Create new raster that changes all cells within the wetland to appropriate NDVI
# 3b. Extract
# 4. Calculate value under both scenarios using function from Liu et al. (2019)
# 5. Aggregate and calculate net difference

  #tm_raster(style = "fixed", breaks=seq(0, 5, by=1), labels = c("<1","1 - 2","2 - 3","3 - 4",">4"), palette = cols, title = "Water Depth (m)",colorNA = "white",textNA = "") +