# download_asos_data.R
# Bryan Holman || v0.1 || 20170502

# This R script downloads KMLB ASOS wind data to match GEFS ensemble data from
# the previous day's 18 UTC run for bias correction/calibration purposes.


# libraries ---------------------------------------------------------------

library(riem) # access to ASOS data through iowa state
library(lubridate) # awesome date handling
library(xts) # also awesome date handling
library(WindVerification) # wind data handling

# paths for /tmp and /util
tmp.path <- paste(getwd(), '/tmp', sep = '')
data.path <- paste(getwd(), '/data', sep = '')
util.path <- paste(getwd(), '/util', sep = '')

# load gefs_all.csv
df.all <- read.csv(paste(data.path, '/', 'gefs_all.csv', sep = ''), 
                   header = TRUE, stringsAsFactors = FALSE)
df.all$runtime <- as.POSIXct(df.all$runtime, tz = 'UTC')
df.all$validtime <- as.POSIXct(df.all$validtime, tz = 'UTC')

# determine what timeframe we need to look at for KMLB ASOS data
times.needs.kmlb <- df.all$validtime[anyNA(c(df.all$kmlb.u, df.all$kmlb.v))]

# get KMLB ASOS data for this time frame
df.kmlb <- riem_measures(station = "MLB", 
                         date_start = format(times.needs.kmlb[1], '%Y-%m-%d'), 
                         date_end = format(Sys.Date() + days(2), '%Y-%m-%d'))

# Only keep the hourly updates, which happen to be the only observations with
# MSLP
df.kmlb <- df.kmlb[!is.na(df.kmlb$mslp),]

# round valid times to nearest quarter hour, n is in seconds
df.kmlb$roundvalid <- as.POSIXct(align.time(df.kmlb$valid, n=60*15))

# convert wind speed from knots to m/s
df.kmlb$wspd <- convertunits(df.kmlb$sknt, inunits = 'knots', outunits = 'm/s')

# get u and v
uv <- mapply(getuv, df.kmlb$wspd, df.kmlb$drct)
df.kmlb$u <- uv[1,]
df.kmlb$v <- uv[2,]

# A loop can't be the best way to do this, but it will work for now!
for (datetime in times.needs.kmlb) {
    if (datetime %in% df.kmlb$roundvalid) {
        df.all[df.all$validtime == datetime, c('kmlb.u', 'kmlb.v')] <- 
            df.kmlb[df.kmlb$roundvalid == datetime, c('u', 'v')]   
    }
}

# now save df.all
write.csv(df.all, file = paste(data.path, '/gefs_all.csv', sep = ''), 
          row.names = FALSE)