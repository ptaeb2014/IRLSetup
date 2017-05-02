# download_gefs_data.R
# Bryan Holman || v0.2 || 20170428

# This R script downloads 0.5° GEFS ensemble wind data from the previous day's 
# 18 UTC run for forecasts hours 00 through 81. U and V wind components are
# extracted and interpolated to Melbourne Int'l Airport (KMLB) for each 
# forecast hour, compiled into a dataframe and then saved to disk in the form
# of a .csv file.

# libraries ---------------------------------------------------------------

library(lubridate) # nice date handling

# functions ---------------------------------------------------------------

# grab the url (through ncep's ftp site) of the grib2 we want to download
getGRIBurl <- function(ens.mem, date, run, fcst.hour) {
    return(paste('http://www.ftp.ncep.noaa.gov/data/nccf/com/gens/prod/gefs.', 
                 date, '/', run, '/pgrb2ap5/', ens.mem, '.t', run, 
                 'z.pgrb2a.0p50.f', fcst.hour, sep = ''))
}

# grab the corresponding .idx file to the grib2 we want to download
getIDXurl <- function(ens.mem, date, run, fcst.hour) {
    return(paste(getGRIBurl(ens.mem, date, run, fcst.hour), '.idx', sep = ''))
}

# using get_inv.pl and get_grib.pl, download the grib2 file of interest, return
# the path of the downloaded file
downloadGRIB <- function(get_inv.path, get_grib.path, ens.mem, date, run, 
                         fcst.hour, outpath) {
    outfile <- paste(outpath, '/', ens.mem, '_', date, '_', run, '_', 
                     fcst.hour, '.grb2', sep = '')
    sys.command <- paste(get_inv.path, getIDXurl(ens.mem, date, run, fcst.hour), 
                         '| grep "10 m above" |', get_grib.path, 
                         getGRIBurl(ens.mem, date, run, fcst.hour), outfile)
    system(sys.command)
    return(outfile)
}

# employ the wgrib2 -small_grid command to reduce the downloaded grib2 file
# to only contain the four closest points to KMLB, return the path of the
# trimmed file
trimGRIB <- function(wgrib2.path, file, lats, lons) {
    trimmed.file <- paste(gsub('.grb2', '', file), 'trimmed.grb2', sep = '_')
    sys.command <- paste(wgrib2.path, ' ', file, ' -small_grib ', lons[1], ':', 
                         lons[2], ' ', lats[1], ':', lats[2], ' ', 
                         trimmed.file, sep = '')
    system(sys.command)
    return(trimmed.file)
}

# given a forecast hour (integer), return a three digit string of said forecast
# hour
getFcstHrString <- function(fcst.hour) {
    if (fcst.hour < 10) {
        return(paste('00', fcst.hour, sep = ''))
    }
    else if (fcst.hour >= 10 & fcst.hour < 100) {
        return(paste('0', fcst.hour, sep = ''))
    }
    else {
        return(paste(fcst.hour, sep = ''))
    }
}

# convert from radians to degrees
deg2rad <- function(degrees) {
    return((degrees * pi) / 180)
}

# interpolate some forecast variable (values) to a lat and lon of interest
interpolate <- function(lat, lon, lats, lons, values) {
    # gefs lons are always positive
    lons <- lons - 360
    # calculate distance between each grid cell
    dists <- sqrt(((lats - lat) * 111.32)^2 + ((lons - lon) * 111.32 * 
                                                   cos(deg2rad(lat)))^2)
    # inverse distance weighting
    weights <- 1 / dists
    return(sum(weights * values) / sum(weights))
}

# global variables --------------------------------------------------------

# determine wgrib2 path
wgrib2.path <- system('which wgrib2', intern = TRUE)

# paths for /tmp and /util
tmp.path <- paste(getwd(), '/tmp', sep = '')
data.path <- paste(getwd(), '/data', sep = '')
util.path <- paste(getwd(), '/util', sep = '')

# get get_inv.pl and get_grib.pl paths
get_inv.path <- paste(util.path, '/get_inv.pl', sep = '')
get_grib.path <- paste(util.path, '/get_grib.pl', sep = '')

# latitude and longitude for point of interest, defaults to KMLB
lat <- 28.102778 # KMLB
lon <- -80.645278 # KMLB

# TODO figure out how to automatically get lats on lons given these coordinates
lons <- seq(279, 279.5, by = 0.5)
lats <- seq(28, 28.5, by = 0.5)

# ensemble members
ens.mems <- c('gec00', 'gep01', 'gep02', 'gep03', 'gep04', 'gep05', 'gep06', 
              'gep07', 'gep08', 'gep09', 'gep10', 'gep11', 'gep12', 'gep13', 
              'gep14', 'gep15', 'gep16', 'gep17', 'gep18', 'gep19', 'gep20')

# date/run to grab data (yesterday, 18z)
# date <- format(Sys.Date() - days(1), '%Y%m%d')
# fcst.time <- as.POSIXct(Sys.Date() - days(1) + hours(18), tz = 'GMT')
date <- format(Sys.Date() - days(5), '%Y%m%d')
fcst.time <- as.POSIXct(Sys.Date() - days(5) + hours(18), tz = 'GMT')
run <- '18'

# dataframe to store final information
df.run <- data.frame(runtime = rep(fcst.time, 28), 
                     fcsthour = seq(0, 81, by = 3))
df.run$validtime <- df.run$runtime + hours(df.run$fcsthour)

# loop through all the ensemble members and download all the data
for (ens.mem in ens.mems) {
    # for (ens.mem in 'gep01') { # for testing purposes
    
    # vectors to store the interpolated u and v forecasts for this ensemble
    # member
    mem.u <- NULL
    mem.v <- NULL
    
    # loop through all forecast hours and download data
    for (fcst.hour in seq(0, 81, by = 3)) {
        # for (fcst.hour in 0:0) { # for testing purposes
        
        # download the file for this ensemble member and forecast hour
        gefs.file <- downloadGRIB(get_inv.path, get_grib.path, ens.mem, date, 
                                  run, getFcstHrString(fcst.hour), tmp.path)
        
        # trim the .grb2 file to only contain 4 closest cells to KMLB
        gefs.trimmed <- trimGRIB(wgrib2.path, gefs.file, lats, lons)
        
        # load in u and v information from this trimmed file
        wgrib2.command <- paste(wgrib2.path, ' ', gefs.trimmed, 
                                ' -match UGRD -spread ', getwd(), '/tmp/u.csv', 
                                sep = '')
        system(wgrib2.command)
        wgrib2.command <- paste(wgrib2.path, ' ', gefs.trimmed, 
                                ' -match VGRD -spread ', getwd(), '/tmp/v.csv', 
                                sep = '')
        system(wgrib2.command)
        
        # now load these .csv files in
        df.u <- read.csv(paste(getwd(), '/tmp/u.csv', sep = ''), header = TRUE, 
                         stringsAsFactors = FALSE)
        df.v <- read.csv(paste(getwd(), '/tmp/v.csv', sep = ''), header = TRUE, 
                         stringsAsFactors = FALSE)
        
        # get interpolated u and v values, save them to appropriate vectors
        mem.u <- c(mem.u, interpolate(lat, lon, df.u$lat, df.u$lon, df.u[[3]]))
        mem.v <- c(mem.v, interpolate(lat, lon, df.v$lat, df.v$lon, df.v[[3]]))
        
        # free up some memory and remove files we don't need anymore
        rm(df.u, df.v)
        file.remove(paste(getwd(), '/tmp/u.csv', sep = ''))
        file.remove(paste(getwd(), '/tmp/v.csv', sep = ''))
        file.remove(gefs.file)
        file.remove(gefs.trimmed)
        
    }
    
    # save the completed mem.u and mem.v vectors to df.run
    df.run[[paste(ens.mem, 'u', sep = '.')]] <- mem.u
    df.run[[paste(ens.mem, 'v', sep = '.')]] <- mem.v
}

# save df.run to disk
write.csv(df.run, file = paste(data.path, '/', 'gefs_', date, '.csv', sep = ''), 
          row.names = FALSE)
