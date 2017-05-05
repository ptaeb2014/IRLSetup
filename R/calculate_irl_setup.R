# calculate_irl_setup.R
# Bryan Holman || v0.1 || 20170501

# Given the GEFS ensemble data just downloaded, this R script calculates the
# setup according to the procedure outlined in Colvin et al. (2016).

# v0.1 -> Calculates wind run for each GEFS ensemble member using the
# interpolated data only. No bias correction or ensemble calibration performed
# as of yet.

# libraries ---------------------------------------------------------------

library(rmarkdown) # rendering index.Rmd at the end
library(riem) # access to ASOS data through iowa state
library(lubridate) # awesome date handling
library(xts) # also awesome date handling
library(WindVerification) # wind data handling

# functions ---------------------------------------------------------------

# calculates the 12-hour wind run given a series of u and v forecasts, returns
# the u and v wind run components
getWindRun <- function(times, us, vs, type = 'model') {
    # create a dataframe given input variables
    df <- data.frame(time = times, u = us, v = vs)
    
    # we can't calculate the 12 hour wind run until we have 12 hours of data!
    # So don't start calculating until this starting point
    windRun.start <- df$time[1] + hours(9)
    
    # if we are doing KMLB obs
    if (type == 'asos') windRun.start <- df$time[1] + hours(11)
    
    df.length <- length(df$time)
    
    # by default, all wind runs are NAs
    windRun.u <- rep(NA, df.length)
    windRun.v <- rep(NA, df.length)
    
    # loop through each line and calculate the wind runs
    for (i in 1:df.length) {
        
        # go to next line if we haven't reach windRun.start yet!
        if (df$time[i] < windRun.start) next
        
        # if we have reach windRun.start, then calculate wind run
        windRun.uv <- colMeans(df[df$time <= df$time[i] & 
                                      df$time > df$time[i] - hours(12), 
                                  c('u', 'v')], na.rm = TRUE)
        
        # save current wind run to vectors
        windRun.u[i] <- windRun.uv[1]
        windRun.v[i] <- windRun.uv[2]
    }
    
    # return windruns as a list
    return(list(u = windRun.u, v = windRun.v))
}

# convert from radians to degrees
deg2rad <- function(degrees) {
    return((degrees * pi) / 180)
}

# calculate the setup 
getSetup <- function(times, us, vs, type = 'model') {
    
    # calculate 12 hour wind run
    windRun <- getWindRun(times, us, vs, type)
    
    # get wind speed and wind direction for these wind runs
    windRun.spd.dir <- mapply(getwspdwdir, windRun$u, windRun$v)
    
    # calculate u1 and u2, the irl-oriented wind components for the northern 
    # and southern irl, respectively
    u1 <- unlist(windRun.spd.dir[1,]) * 
        cos(sapply(unlist(windRun.spd.dir[2,]) + 10, deg2rad))
    u2 <- unlist(windRun.spd.dir[1,]) * 
        cos(sapply(unlist(windRun.spd.dir[2,]) + 26, deg2rad))
    
    # now calculate the lagoon relative wind
    u.r <- -1 * (u1 * 28 + u2 * 70) / (28 + 70)
    # u.r <- (u1 * 28 + u2 * 70) / (28 + 70)

    # return the setup given u.r
    return(0.728 * sign(u.r) * abs(u.r)^1.5 - 3)
}

# global variables --------------------------------------------------------

# path for /data
data.path <- paste(getwd(), '/data', sep = '')

# ensemble members
ens.mems <- c('gec00', 'gep01', 'gep02', 'gep03', 'gep04', 'gep05', 'gep06', 
              'gep07', 'gep08', 'gep09', 'gep10', 'gep11', 'gep12', 'gep13', 
              'gep14', 'gep15', 'gep16', 'gep17', 'gep18', 'gep19', 'gep20')

# # most recent data (yesterday)
# date <- format(Sys.Date() - days(1), '%Y%m%d')

# data --------------------------------------------------------------------

# open df.all, we will need a few days of gefs forecasts
df.all <- read.csv(paste(data.path, '/gefs_all.csv', sep = ''), 
                         header = TRUE, stringsAsFactors = FALSE)
df.all$runtime <- as.POSIXct(df.all$runtime, tz = 'UTC')
df.all$validtime <- as.POSIXct(df.all$validtime, tz = 'UTC')

# we will only calculate setup for the last 4 gefs runs
gefs.runs <- tail(unique(df.all$runtime), 4)

# get data frames for only these four runs, and neglect kmlb.u & kmlb.v columns
df.gefs.1 <- df.all[df.all$runtime == gefs.runs[1],-c(46, 47)]
df.gefs.2 <- df.all[df.all$runtime == gefs.runs[2],-c(46, 47)]
df.gefs.3 <- df.all[df.all$runtime == gefs.runs[3],-c(46, 47)]
df.gefs.recent <- df.all[df.all$runtime == gefs.runs[4],-c(46, 47)]

# clear up some memory
rm(df.all)

# data frames to store setup information
gefs.setup.1 <- data.frame(validtime = df.gefs.1$validtime)
gefs.setup.2 <- data.frame(validtime = df.gefs.2$validtime)
gefs.setup.3 <- data.frame(validtime = df.gefs.3$validtime)
gefs.setup.recent <- data.frame(validtime = df.gefs.recent$validtime)

# loop through all the ensemble members and calculate the wind runs
for (ens.mem in ens.mems) {

    # calculate setup for each gefs run and add it to the respective setup
    # data frame
    gefs.setup.1[[paste(ens.mem, 'raw', sep = '.')]] <- 
        getSetup(times = df.gefs.1$validtime, 
                 df.gefs.1[[paste(ens.mem, 'u', sep = '.')]], 
                 df.gefs.1[[paste(ens.mem, 'v', sep = '.')]])
    gefs.setup.2[[paste(ens.mem, 'raw', sep = '.')]] <- 
        getSetup(times = df.gefs.2$validtime, 
                 df.gefs.2[[paste(ens.mem, 'u', sep = '.')]], 
                 df.gefs.2[[paste(ens.mem, 'v', sep = '.')]])
    gefs.setup.3[[paste(ens.mem, 'raw', sep = '.')]] <- 
        getSetup(times = df.gefs.3$validtime, 
                 df.gefs.3[[paste(ens.mem, 'u', sep = '.')]], 
                 df.gefs.3[[paste(ens.mem, 'v', sep = '.')]])
    gefs.setup.recent[[paste(ens.mem, 'raw', sep = '.')]] <- 
        getSetup(times = df.gefs.recent$validtime, 
                 df.gefs.recent[[paste(ens.mem, 'u', sep = '.')]], 
                 df.gefs.recent[[paste(ens.mem, 'v', sep = '.')]])
}
# clear up some memory
rm(df.gefs.1, df.gefs.2, df.gefs.3, df.gefs.recent)

# get hourly KMLB ASOS data for the last few days to plot
df.kmlb <- riem_measures(station = "MLB", 
                         date_start = format(Sys.Date() - days(4), '%Y-%m-%d'), 
                         date_end = format(Sys.Date() + days(1), '%Y-%m-%d'))

# Only keep the hourly updates, which happen to be the only observations with
# MSLP
df.kmlb <- df.kmlb[!is.na(df.kmlb$mslp),]

# round valid times to nearest quarter hour, n is in seconds
df.kmlb$roundvalid <- as.POSIXct(align.time(df.kmlb$valid, n=60*15))
df.kmlb$roundvalid <- as.POSIXct(format(df.kmlb$roundvalid, 
                                        '%Y-%m-%d %H:%M:%S'), tz = 'UTC')

# convert wind speed from knots to m/s
df.kmlb$wspd <- convertunits(df.kmlb$sknt, inunits = 'knots', outunits = 'm/s')

# get u and v
uv <- mapply(getuv, df.kmlb$wspd, df.kmlb$drct)
df.kmlb$u <- uv[1,]
df.kmlb$v <- uv[2,]

# calculate KMLB ASOS setup
kmlb.setup <- getSetup(df.kmlb$roundvalid, df.kmlb$u, df.kmlb$v, 
                       type = 'asos')

# create a dataframe with just this information
asos.setup <- data.frame(roundvalid = df.kmlb$roundvalid, setup = kmlb.setup)
save(gefs.setup.1, gefs.setup.2, gefs.setup.3, gefs.setup.recent, asos.setup, 
     file = 'data/setup.RData')
rmarkdown::render('R/index.Rmd', output_dir = '~/Dropbox/IRLSetup/docs/')
# save df.run to disk
# write.csv(df.run, file = paste(data.path, '/', 'gefs_', date, '.csv', sep = ''), 
#           row.names = FALSE)
