# analyze_irl_setup.R
# Bryan Holman // v0.2 // 20170503

# This is the testbed for creating the perfect ggplot to go on the site!

library(reshape2)
library(ggplot2)
library(lubridate)

# functions ---------------------------------------------------------------

# given a data frame of ensemble data, add ens average, max, min, and median
addEnsStats <- function(df.ens, cols = 2:22) {
    df.ens$avg <- rowMeans(df.ens[cols], na.rm = FALSE)
    df.ens$min <- apply(df.ens[cols], 1, min)
    df.ens$max <- apply(df.ens[cols], 1, max)
    df.ens$med <- apply(df.ens[cols], 1, median)
    return(df.ens)
}


# df.setup.melt <- melt(df.setup, id.vars = 'validtime')

# # the code below does okay ... not awesome ... too many lines!
# p <- ggplot(df.setup.melt, aes(x = validtime, y = value, color = variable)) + 
#     geom_line() + theme_minimal() + theme(legend.position="bottom")
load('~/Dropbox/IRLSetup/data/setup.RData')
# remove the missing rows in gefs.setup
gefs.setup.1 <- gefs.setup.1[4:28,]
gefs.setup.2 <- gefs.setup.2[4:28,]
gefs.setup.3 <- gefs.setup.3[4:28,]
gefs.setup.recent <- gefs.setup.recent[4:28,]
# match asos times
asos.setup <- asos.setup[asos.setup$roundvalid >= gefs.setup.1$validtime[1],]
gefs.setup.recent.melt <- melt(gefs.setup.recent, id.vars = 'validtime')
gefs.setup.1 <- addEnsStats(gefs.setup.1)
gefs.setup.2 <- addEnsStats(gefs.setup.2)
gefs.setup.3 <- addEnsStats(gefs.setup.3)
gefs.setup.recent <- addEnsStats(gefs.setup.recent)
date.breaks <- seq.POSIXt(asos.setup$roundvalid[1] - hours(3), 
                          gefs.setup.recent$validtime[24], by = '12 hours')
p <- ggplot(gefs.setup.recent, aes(x = validtime)) + 
    geom_hline(aes(yintercept = 0), linetype = 'dashed') + 
    geom_vline(aes(xintercept = as.numeric(gefs.setup.recent$validtime[1])), 
               linetype = 'dashed') +
    geom_line(data = gefs.setup.1, aes(y = max), alpha = 0.5) + 
    # geom_line(data = gefs.setup.1, aes(y = avg), alpha = 0.5) + 
    geom_line(data = gefs.setup.1, aes(y = min), alpha = 0.5) + 
    geom_ribbon(data = gefs.setup.1, aes(ymin = min, ymax = max), 
                alpha = 0.125) + 
    geom_line(data = gefs.setup.2, aes(y = max), alpha = 0.5) + 
    # geom_line(data = gefs.setup.2, aes(y = avg), alpha = 0.5) + 
    geom_line(data = gefs.setup.2, aes(y = min), alpha = 0.5) + 
    geom_ribbon(data = gefs.setup.2, aes(ymin = min, ymax = max), 
                alpha = 0.125) + 
    geom_line(data = gefs.setup.3, aes(y = max), alpha = 0.5) + 
    # geom_line(data = gefs.setup.3, aes(y = avg), alpha = 0.5) + 
    geom_line(data = gefs.setup.3, aes(y = min), alpha = 0.5) + 
    geom_ribbon(data = gefs.setup.3, aes(ymin = min, ymax = max), 
                alpha = 0.125) + 
    geom_line(data = gefs.setup.recent.melt, 
              mapping = aes(x = validtime, y = value, color = variable)) + 
    geom_line(aes(y = gec00.raw, color = 'Ensemble Members')) + 
    geom_line(aes(y = med, color = 'Ensemble Median'), size = 1.5) + 
    geom_line(aes(y = avg, color = 'Ensemble Mean'), size = 1.5) + 
    geom_line(aes(y = max)) + geom_line(aes(y = min)) + 
    geom_ribbon(aes(ymin = min, ymax = max, fill = 'Ensemble Spread'), 
                alpha = 0.25) + 
    geom_point(data = asos.setup, mapping = aes(x = roundvalid, y = setup), 
              size = 3, fill = 'orange', shape = 21) +
    scale_color_manual(breaks = c('Ensemble Median', 'Ensemble Mean', 
                                  'Ensemble Members'), 
                       values = c('red', 'blue', 'grey', rep('grey', 21))) +
    scale_fill_manual(breaks = c('Ensemble Spread'), 
                      values = c('black')) +
    theme_light() + xlab('') + ylab('IRL Setup (cm)') + 
    scale_x_datetime(breaks = date.breaks,
                     date_labels = '%b %d\n %H UTC',
                     limits = c(date.breaks[1], NA)) +
    theme(legend.position="bottom", legend.title = element_blank())
print(p)
