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
    scale_color_grey(start = 0.7, end = 0.7) +
    geom_line(aes(y = med.raw), color = 'red', size = 2) + 
    geom_line(aes(y = avg.raw), color = 'blue', size = 2) + 
    geom_line(aes(y = max.raw)) + geom_line(aes(y = min.raw)) + 
    geom_ribbon(aes(ymin = min.raw, ymax = max.raw), alpha = 0.25) + 
    geom_hline(aes(yintercept = 0), linetype = 'dashed') + theme_light() + 
    xlab('') + ylab('IRL Setup (cm)') + theme(legend.position="none")
print(p)