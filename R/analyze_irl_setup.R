library(reshape2)
library(ggplot2)

# df.setup.melt <- melt(df.setup, id.vars = 'validtime')

# # the code below does okay ... not awesome ... too many lines!
# p <- ggplot(df.setup.melt, aes(x = validtime, y = value, color = variable)) + 
#     geom_line() + theme_minimal() + theme(legend.position="bottom")

df.test <- df.setup[4:28,]
df.setup.melt <- melt(df.test, id.vars = 'validtime')
df.test$avg.raw <- rowMeans(df.test[-1], na.rm = FALSE)
df.test$min.raw <- apply(df.test[-c(1, 23)], 1, min)
df.test$max.raw <- apply(df.test[-c(1, 23)], 1, max)
df.test$med.raw <- apply(df.test[-c(1, 23)], 1, median)
p <- ggplot(df.test, aes(x = validtime)) + 
    geom_line(data = df.setup.melt, mapping = aes(x = validtime, y = value, 
                                                  color = variable)) + 
    # scale_color_grey(start = 0.7, end = 0.7, guide = 'none') +
    geom_line(aes(y = gec00.raw, color = 'Ensemble Members')) + 
    geom_line(aes(y = med.raw, color = 'Ensemble Median'), size = 1.5) + 
    geom_line(aes(y = avg.raw, color = 'Ensemble Mean'), size = 1.5) + 
    geom_line(aes(y = max.raw)) + geom_line(aes(y = min.raw)) + 
    geom_ribbon(aes(ymin = min.raw, ymax = max.raw, fill = 'Ensemble Spread'), 
                alpha = 0.25) + 
    scale_color_manual(breaks = c('Ensemble Median', 'Ensemble Mean', 
                                  'Ensemble Members'), 
                       values = c('red', 'blue', 'grey', rep('grey', 21))) +
    scale_fill_manual(breaks = c('Ensemble Spread'), 
                      values = c('black')) +
    geom_hline(aes(yintercept = 0), linetype = 'dashed') + theme_light() + 
    xlab('') + ylab('IRL Setup (cm)') + 
    theme(legend.position="bottom", legend.title = element_blank())
print(p)
ggsave('docs/img/raw_setup.png', width = 8, height = 6, units = 'in', dpi = 150)
