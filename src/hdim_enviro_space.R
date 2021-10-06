library(sp)
library(raster)
library(rgdal)

# ----
# get hdim site locations

# read in data that includes long, lat
x <- read.csv('data/hdimbeat_abund.csv')

# extract long, lat for sites (averaging over subplots)
loc <- lapply(split(x, x$site), function(y) {
    return(c(mean(y$Long1), mean(y$Lat1)))
})

loc <- do.call(rbind, loc)

# convert to spatial object with proper CRS
crs <- '+init=epsg:4269 +proj=longlat +ellps=GRS80'
hdimSites <- SpatialPointsDataFrame(coords = loc,
                                    data = data.frame(name = rownames(loc)),
                                    proj4string = CRS(crs))



# ----
# get precip data

# url for rainfall atlas raw data
rainURL <- 'http://rainfall.geography.hawaii.edu/assets/files/GISLayers/StateASCIIGrids_mm.zip'

# make local folder to store
rainDir <- 'data/precip'
if(!dir.exists(rainDir)) {
    dir.create(rainDir)
}

# download and unzip
rainFile <- file.path(rainDir, 'precip_mm.zip')
download.file(rainURL, destfile = rainFile)
unzip(rainFile, exdir = rainDir)

# read-in annual rainfall as a raster
precipRast <- raster('data/precip/rfgrid_mm_state_ann.txt')

extract(precipRast,
        spTransform(hdimSites, CRS(proj4string(precipRast))))


sum(!is.na(values(precipRast)))

# ----
# read in elevation data


'http://gis.ess.washington.edu/data/raster/tenmeter/hawaii/bigisland.zip'
'http://gis.ess.washington.edu/data/raster/tenmeter/hawaii/kahoolawe.zip'
'http://gis.ess.washington.edu/data/raster/tenmeter/hawaii/kauai.zip'
'http://gis.ess.washington.edu/data/raster/tenmeter/hawaii/lanai.zip'
'http://gis.ess.washington.edu/data/raster/tenmeter/hawaii/maui.zip'
'http://gis.ess.washington.edu/data/raster/tenmeter/hawaii/molokai.zip'
'http://gis.ess.washington.edu/data/raster/tenmeter/hawaii/niihau.zip'
'http://gis.ess.washington.edu/data/raster/tenmeter/hawaii/oahu.zip'
