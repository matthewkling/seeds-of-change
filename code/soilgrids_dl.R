library(terra)
# library(gdalUtilities)
# library(raster)
library(tidyverse)


# SoilGrids v1 (2017) ##############

"https://files.isric.org/soilgrids/former/2017-03-10/data/AWCh1_M_sl1_250m_ll.tif"

vars <- tibble(var = c("AWCh1", "AWCh2", "AWCh3", "AWCtS", "BLDFIE", "CECSOL", "CLYPPT", "CRFVOL",
                       "OCDENS", "OCSTHA", "ORCDRC", "PHIHOX", "PHIKCL", "SLTPPT", "SNDPPT", "TEXMHT", "WWP")) %>%
      expand_grid(depth = c(1, 3, 5, 7))

dl <- function(var, depth){
      x <- paste0(var, "_M_sl", depth, "_250m_ll.tif")
      dest <- paste0("~/data/soilgrids/v1/250m/", basename(x))
      if(file.exists(dest)) return("skipped")
      base <- "https://files.isric.org/soilgrids/former/2017-03-10/data/"
      options(timeout = 100000)
      download.file(paste0(base, x), dest)
      download.file(paste0(base, x, ".xml"), paste(dest, ".xml"))
}

pmap(vars, dl)


# crop & project to CA

r <- rast("assets/climate/ensemble 1981-2010 NA AET.tif")[[1]] # template (chelsa)

f <- list.files("~/data/soilgrids/v1/250m", full.names = T)
f <- f[!grepl("xml", f)]

f %>% map(function(x){
      x %>% 
            rast() %>% 
            crop(r) %>% 
            project(y = r, method = "average") %>%
            mask(r) %>%
            writeRaster(paste0("~/data/soilgrids/v1/250m_California/", basename(x)))
})


      
##### deprecated below ########



# 1000m, v2 ################

vars <- tibble(var = c("bdod", "cec", "clay", "cfvo", "nitrogen", 
                       "ocd", "phh2o", "sand", "silt", "soc"),
               name = c("Bulk density of the fine earth fraction", "Cation Exchange Capacity of the soil", "Proportion of clay particles (< 0.002 mm) in the fine earth fraction", "Volumetric fraction of coarse fragments (> 2 mm)", "Total nitrogen (N)", 
                        "Organic carbon density", "Soil pH", "Proportion of sand particles (> 0.05 mm) in the fine earth fraction", "Proportion of silt particles (>= 0.002 mm and <= 0.05 mm) in the fine earth fraction", "Soil organic carbon content in the fine earth fraction"),
               units = c("cg/cm3", "mmol(c)/kg", "g/kg", "cm3/dm3", "cg/kg",
                         "hg/m3", "ph*10", "g/kg", "g/kg", "dg/kg")) %>%
      expand_grid(depth = c("0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm")) %>%
      mutate(path = paste0(var, "/", var, "_", depth, "_mean"))

# download global 1000m rasters
sg_dl_1000m <- function(x){
      # https://files.isric.org/soilgrids/latest/data_aggregated/1000m/bdod/bdod_0-5cm_mean_1000.tif
      message(x)
      dest <- paste0("~/data/soilgrids/v2/1000m/", basename(x), ".tif")
      if(file.exists(dest)) return("skipped")
      base <- "https://files.isric.org/soilgrids/latest/data_aggregated/1000m/"
      url <- paste0(base, x, "_1000.tif")
      options(timeout = 10000)
      download.file(url, dest)
}
map(vars$path, sg_dl_1000m)


template <- rast("assets/climate/ensemble 1981-2010 NA AET.tif")[[1]]
prj <- function(x, template){
      dest <- paste0("~/data/soilgrids/v2/1000m_ll_california/", basename(x))
      if(file.exists(dest)) return("skipped")
      s <- rast(x)
      sr <- try(terra::project(s, template))
      if(class(sr) == "SpatRaster") writeRaster(sr, dest)
}

list.files("~/data/soilgrids/v2/1000m", full.names = T) %>%
      map(prj, template = template)



# 250m #################


r <- rast("assets/climate/ensemble 1981-2010 NA AET.tif")[[1]] # template (chelsa)
# plot(r)
# e <- draw("extent")
# f <- as.vector(e*1.1) %>% matrix(nrow = 2) %>%
#       as.data.frame()
# coordinates(f) <- c("V1", "V2")
# crs(f) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
# f <- f %>% spTransform('+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs')
# bb <- as.vector(ext(f))[c(1, 4, 2, 3)] # ulx, uly, lrx, lry
bb <- c(-13488471, 4731216, -12254427, 3565411)

igh = '+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs' # proj string for Homolosine projection
sg_url = "/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=https://files.isric.org/soilgrids/latest/data/"

expand_grid(var = c("ocs"),
            name = "organic carbon stock",
            units = "t/ha",
            depth = c("0-30cm"))

vars <- tibble(var = c("bdod", "cec", "clay", "cfvo", "nitrogen", 
                       "ocd", "phh20", "sand", "silt", "soc"),
               name = c("bulk densith", "cation exchange capacity", "clay", "coarse fragments", "nitrogen", 
                        "organic carbon density", "pH water", "sand", "silt", "soil organic carbon"),
               units = c("cg/cm3", "mmol(c)/kg", "g/kg", "cm3/dm3", "cg/kg", 
                         "hg/m3", "ph*10", "g/kg", "g/kg", "dg/kg")) %>%
      expand_grid(depth = c("0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm")) %>%
      mutate(path = paste0(var, "/", var, "_", depth, "_mean"))

# download 250m rasters for california, converting to tif format
sg_dl <- function(x){
      message(x)
      dest <- paste0("soilgrids/igh_proj/", basename(x), ".tif")
      if(file.exists(dest)) return("skipped")
      gdal_translate(paste0(sg_url, x, ".vrt"),
                     dest,
                     tr = c(250, 250),
                     projwin = bb,
                     projwin_srs = igh)
}
map(vars$path, sg_dl)

prj <- function(x, template){
      dest <- paste0("soilgrids/ll_proj/", basename(x))
      if(file.exists(dest)) return("skipped")
      s <- rast(x)
      sr <- terra::project(s, template, method = "average")
      writeRaster(sr, dest)
}

list.files("soilgrids/igh_proj/", full.names = T) %>%
      map(prj, template = r)

