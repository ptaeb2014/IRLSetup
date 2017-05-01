# calculate_irl_setup.R
# Bryan Holman || v0.1 || 20170501

# Given the GEFS ensemble data just downloaded, this R script calculates the
# setup according to the procedure outlined in Colvin et al. (2016).

# v0.1 -> Calculates wind run for each GEFS ensemble member using the
# interpolated data only. No bias correction or ensemble calibration performed
# as of yet.

# libraries ---------------------------------------------------------------

library(lubridate) # nice date handling

# functions ---------------------------------------------------------------

# calculates the 12-hour wind run given a series of u and v forecasts, returns
# the u and v wind run components
getWindRun <- function(times, us, vs) {
    # create a dataframe given input variables
    df <- data.frame(time = times, u = us, v = vs)
    
    # we can't calculate the 12 hour wind run until we have 12 hours of data!
    # So don't start calculating until this starting point
    windRun.start <- df$time[1] + hours(9)
    
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

# get the wind speed and direction from u and v components
getwspdwdir <- function(u, v) {
    
    # If either u or v are missing, we cannot do the calculation!
    if (anyNA(c(u, v))) {return(c(NA, NA))}
    
    # Calculate wind speed and wind direction
    wspd <- sqrt(u^2 + v^2)
    wdir <- atan2(-1 * u, -1 * v) * 57.2957795131 # 57.2957795131 is 180 / pi
    
    # atan2 goes from -pi to pi, so you need to add 360 if negative
    if (wdir < 0) {wdir <- wdir + 360}
    if (round(wdir) >= 360) {wdir <- 0}
    
    # return wind speed first, then wind direction rounded to nearest integer
    return(list(wspd = wspd, wdir = round(wdir, 0)))
}

# global variables --------------------------------------------------------

# path for /data
data.path <- paste(getwd(), '/data', sep = '')

# ensemble members
ens.mems <- c('gec00', 'gep01', 'gep02', 'gep03', 'gep04', 'gep05', 'gep06', 
              'gep07', 'gep08', 'gep09', 'gep10', 'gep11', 'gep12', 'gep13', 
              'gep14', 'gep15', 'gep16', 'gep17', 'gep18', 'gep19', 'gep20')

# most recent data (yesterday)
date <- format(Sys.Date() - days(1), '%Y%m%d')

# data --------------------------------------------------------------------

# in v0.1, all we need is yesterday's data
df.recent <- read.csv(paste(data.path, '/gefs_', date, '.csv', sep = ''), 
                      header = TRUE, stringsAsFactors = FALSE)
df.recent$runtime <- as.POSIXct(df.recent$runtime, tz = 'GMT')
df.recent$validtime <- as.POSIXct(df.recent$validtime, tz = 'GMT')

# dataframe to store setup information
df.setup <- data.frame(validtime = df.recent$validtime)

# loop through all the ensemble members and calculate the wind runs
for (ens.mem in ens.mems) {
# for (ens.mem in 'gec00') { # for testing purposes
    
    # calculate 12 hour wind run given ensemble us and vs
    mem.windruns <- getWindRun(times = df.setup$validtime, 
                               df.recent[[paste(ens.mem, 'u', sep = '.')]], 
                               df.recent[[paste(ens.mem, 'v', sep = '.')]])
    
    # get wind speed and wind direction for these wind runs
    mem.wspd.wdir <- mapply(getwspdwdir, mem.windruns$u, mem.windruns$v)
    
    # calculate u1 and u2, the irl-oriented wind components for the northern and
    # southern irl, respectively
    u1 <- unlist(mem.wspd.wdir[1,]) * cos(sapply(unlist(mem.wspd.wdir[2,]) + 10, 
                                                 deg2rad))
    u2 <- unlist(mem.wspd.wdir[1,]) * cos(sapply(unlist(mem.wspd.wdir[2,]) + 26, 
                                                 deg2rad))
    
    # now calculate the lagoon relative wind
    u.r <- -1 * (u1 * 28 + u2 * 70) / (28 + 70)
    
    # calculate setup and add it to df.setup
    df.setup[[paste(ens.mem, 'raw', sep = '.')]] <- 0.728 * sign(u.r) * 
        abs(u.r)^1.5 - 3
}

# save df.run to disk
# write.csv(df.run, file = paste(data.path, '/', 'gefs_', date, '.csv', sep = ''), 
#           row.names = FALSE)
