## 3D plot in R - through plotly, can make an interactive plot!

#install.packages("plotly")
library(plotly)
library(tidyverse)

Sys.setenv("plotly_username"="dasomo")
Sys.setenv("plotly_api_key"="")

## Example from plotly website
#p <- plot_ly(mtcars, x = ~wt, y = ~hp, z = ~qsec, color = ~am, colors = c('#BF382A', '#0C4B8E')) %>%
#  add_markers() %>%
#  layout(scene = list(xaxis = list(title = 'Weight'),
#                      yaxis = list(title = 'Gross horsepower'),
#                      zaxis = list(title = '1/4 mile time')))


pcrit_3d <- plot_ly(md_all_temps, x = ~temp, y = ~avg_smr, z = ~avg_pcrit,
                    color = ~species) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'Temperature'),
                      yaxis = list(title = 'MO2min (umol/g/hr)'),
                      zaxis = list(title = 'Pcrit (torr)')))
pcrit_3d
