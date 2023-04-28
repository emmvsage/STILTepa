#code to convolve emissions with footprints

# Revision from 
# WBB CO2 STILT Tutorial
# for PM2.5 impact 
# Emma V. Sage - Dec 12, 2022

library(dplyr)
library(ggplot2)
library(raster)


start_time <- Sys.time() # program start time

stilt_wd <- file.path('~/myproject')
setwd(stilt_wd)
getwd()
#Part I: read data and save PM2.5 data
# Load emissions inventory, assigning the gridded product to "emissions" and
# extracting the time for each grid to "emissions_time"
emissions <- readRDS('~/myproject/emissions_cal_2018_3kmresol.rds') # umol PM2.5 m-2 hr-1
emissions_time <- getZ(emissions)

# Find all footprint files produced by STILT and 
# select the ones for Mono County control monitor
footprint_all <- list.files("out/footprints", full.names = TRUE, pattern = "*.nc$")
length(footprint_all)
coordsString <- paste0(-119.1203, "_", 37.96207) # for Mono County control monitor
footprint_paths <- grep(coordsString, footprint_all, value = TRUE)
length(footprint_paths)

# For each footprint in "footprint_paths", calculate the PM2.5 contribution 
# from the near-field emisisons

# emissions ug m-2 day-1 to umol m-2 s-1; footprint in umol m-2 s-1 
# 1 day = 24 hours; 1 hr = 3600 seconds; 
# 60% carbon; 10% potassium, chlorine and calcium; 30% hydrogen, oxygen and nitrogen
# PM2.5 mix weight: grams/mole
wght_pm25 <- 0.6 * 12.0107 + 0.1 * (39.0983 + 35.453 + 40.078) / 3.0 + 
  0.3 * (1.00794 + 16.00 + 14.01)/3.0

concentration <- lapply(1:length(footprint_paths), function(i) {
  # Import footprint and extract timestamp
  foot <- brick(footprint_paths[i]) # umol m-2 s-1
  time <- as.POSIXct(getZ(foot), tz = 'UTC', origin = '1970-01-01')
  
  # Convert 3d brick to 2d raster if only a single timestep contains influence
  if (nlayers(foot) == 1)
    foot <- raster(foot, layer = 1)
  
  # Subset emissions to match the footprint timestamps
  band <- findInterval(time, emissions_time)
  emissions_subset <- subset(emissions, band)
  
  # Calculate the near-field PM2.5 contribution by taking the product 
  # of the footprints and the emissions fluxes
  # A mole of a given substance is the number of grams of that substance 
  # that contain 6.022 × 1023 particles (molecules) of that substance
  # emissions ug m-2 day-1 to umol m-2 s-1 to match with the unit of footprint
  emi_ug_to_umol <- 1 * 1/(24 * 60 * 60) * (1/wght_pm25) * 1000000 #1e6 for mole to umol
  data.frame(Time_UTC = max(time) + 3600, PM25 = sum(values(foot * 
                                                              emissions_subset * emi_ug_to_umol), na.rm = T)) #PM2.5 is in ppm

}) %>%
  bind_rows() %>% # ug/m3 = molecular weight x concentration (ppb) ÷ 24.45
  mutate(PM25 = (PM25 * 1000 * wght_pm25 / 24.45)) # ppm to ug m-3

#identify scaling factor
#install.packages("dplyr")
#install.packages("plyr")
library(plyr)

write.table(concentration, file = file.path("Mono", "Mono_control_conc_07_31-08_12.csv"), 
            append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

#method #1
# Plot a timeseries of the modeled concentrations and save the figure

f <- concentration %>%
  ggplot(aes(x = Time_UTC, y = PM25, color = PM25)) +
  geom_line() + geom_point() + 
  scale_color_gradientn(colors = c('blue', 'cyan', 'green', 'yellow', 'orange', 'red'),
                        limits = range(concentration$PM25)) +
  labs(x = 'Time of Day', y = 'PM2.5 Conc [ug/m3]', fill = 'PM2.5') +
  theme_classic()
show(f)
ggsave(filename = file.path("Mono", "timeseries.png"), f,
       width = 8, height = 5, dpi = 300, units = "in", device = "png")


# method #2: Plot a timeseries of the modeled concentrations and save the figure
#library(plotly)
#p <- plot_ly(concentration, x = ~Time_UTC, y = ~PM25,
#             type = 'scatter', mode = 'lines') %>%
#  layout(xaxis = list(title = ''),
#         yaxis = list(title = 'PM<sub>2.5</sub> [ug/m3]'))
#htmlwidgets::saveWidget(p, file = file.path("LosAngeles", "timeseries.html"))

# For each footprint in "footprint_paths", fetch the footprint total into a list
foot_list <- lapply(1:length(footprint_paths), function(i) {
  # Import footprint and extract timestamp
  foot <- brick(footprint_paths[i]) # umol CO2 m-2 s-1
  
  # Convert 3d brick to 2d raster if only a single timestep contains influence
  if (nlayers(foot) == 1) {
    foot <- raster(foot, layer = 1)
  } else {
    foot <- sum(foot)
  }
  
  foot
})
​
# Calculate the total and average footprint from the list of footprints
# method 1
foot_total <- sum(stack(foot_list))
foot_average <- sum(stack(foot_list)) / length(foot_list)
writeRaster(foot_total, file.path("Mono", "foot_total.tif"), 
            format = "GTiff", overwrite = TRUE)
writeRaster(foot_average, file.path("Mono", "foot_average.tif"), 
            format = "GTiff", overwrite = TRUE)

#method 2
#crng <- range(values(foot_average))
#cpal <- colorNumeric('Spectral', domain = crng, reverse = T)

#library(dplyr) # for %>%
#library(leaflet) # for addLegend
#map <- leaflet() %>%
#  addProviderTiles('CartoDB.Positron') %>%
#  addRasterImage(foot_average, opacity = 0.5, colors = cpal) %>%
#  addLegend(position = 'bottomleft',
#            pal = cpal,
#            values = crng,
#            title = paste0('m<sup>2</sup> s ppm<br>',
#                           '<span style="text-decoration:overline">&mu;mol</span>'))
#map
#htmlwidgets::saveWidget(map, file = file.path("LosAngeles", "average_footprint.html"))
​
# For each footprint in "footprint_paths", fetch the convolved flux field using
# emissions * footprints into a list of rasters that represent the contribution
# of fluxes over space
contribution_list <- lapply(1:length(footprint_paths), function(i) {
  # Import footprint and extract timestamp
  foot <- brick(footprint_paths[i]) # umol CO2 m-2 s-1
  time <- as.POSIXct(getZ(foot), tz = 'UTC', origin = '1970-01-01')
  
  # Convert 3d brick to 2d raster if only a single timestep contains influence
  if (nlayers(foot) == 1)
    foot <- raster(foot, layer = 1)
  
  # Subset "emissions" to match the footprint timestamps
  band <- findInterval(time, emissions_time)
  emissions_subset <- subset(emissions, band)
  
  # Calculate the near-field PM2.5 contribution by taking the product of the
  # footprints and the fluxes
  sum(foot * emissions_subset)
})

# Calculate the average contribution from the list of contribution totals
# method 1
contribution_list <- contribution_list[!sapply(contribution_list,is.null)]
contribution_average <- sum(stack(contribution_list)) / length(contribution_list)
contribution_total <- sum(stack(contribution_list)) 
#output hourly individual contribution images
#lapply(1:length(contribution_list), function(i) {
#  writeRaster(contribution_list[[i]], filename=file.path("../data/Temporal_Images", 
#                                                         getZ(contribution_list[[i]])), format = "GTiff", overwrite = TRUE)
#})
#output the total single contribution image
# method 1
writeRaster(contribution_total, file.path("Mono", "contribution_total.tif"), 
            format = "GTiff", overwrite = TRUE)
writeRaster(contribution_average, file.path("Mono", "contribution_average.tif"), 
            format = "GTiff", overwrite = TRUE)

# method 2
#contribution_average <- sum(stack(contribution_list)) / length(contribution_list)
#contribution_average <- log10(contribution_average)
#contribution_average[contribution_average < -5] <- -5

#crng <- range(values(contribution_average))
# crng <- c(-5.001, 0)
#cpal <- colorNumeric('Spectral', domain = crng, reverse = T)

#map <- leaflet() %>%
#  addProviderTiles('CartoDB.Positron') %>%
#  addRasterImage(contribution_average, opacity = 0.5, colors = cpal) %>%
#  addLegend(position = 'bottomleft',
#            pal = cpal,
#            values = crng,
#            title = paste0('log10 <br>',
#                           'm<sup>2</sup> s ppm<br>',
#                           '<span style="text-decoration:overline">&mu;mol</span>'))
#map
#htmlwidgets::saveWidget(map, file = file.path("Mono","average_contribution.html"))

end_time <- Sys.time() # program end time

duration <- end_time - start_time # time elapsed in minutes

print(paste0("Total number of minutes spent: ", duration))
