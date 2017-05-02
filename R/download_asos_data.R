# download_asos_data.R
# Bryan Holman || v0.1 || 20170502

# This R script downloads KMLB ASOS wind data to match GEFS ensemble data from
# the previous day's 18 UTC run for bias correction/calibration purposes.


# libraries ---------------------------------------------------------------

library(riem) # access to ASOS data through iowa state
library(lubridate) # awesome date handling
library(xts) # also awesome date handling

df.kmlb <- riem_measures(station = "MLB", date_start = "2017-04-28", 
                          date_end = "2017-05-03")
df.kmlb <- df.kmlb[!is.na(df.kmlb$mslp),]
# round valid times to nearest quarter hour, n is in seconds
df.kmlb$roundvalid <- as.POSIXct(align.time(df.kmlb$valid, n=60*15))
