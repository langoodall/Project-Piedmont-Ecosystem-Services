library(tidyverse)
library(terra)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
job_index <- as.numeric(args[1])

# Set up objects so that we can then run the for loop
mainDir <- "/share/tcsi/lagoodal/R/MAPS"
scenarios <- c("Resistance", "Resilience", "Transition")
rcps <- c("RCP45", "RCP85")
climates <- c("HSHV", "HSLV", "LSHV", "LSLV")
runs <- 1:5
timesteps <- seq(10, 80, 10)

# Determine scenario/rcp/climate combo
combo_df <- expand.grid(Scenario = scenarios, RCP = rcps, Climate = climates)
this_combo <- combo_df[job_index, ]
scenario <- as.character(this_combo$Scenario)
rcp <- as.character(this_combo$RCP)
climate <- as.character(this_combo$Climate)

cat("Processing combo:", scenario, rcp, climate, "\n")

# For loop to calculate the weighte biomass of fruit and nut bearing trees
strList <- list()

for (run in runs) {
  for (year in timesteps) {
    runDir <- file.path(mainDir, scenario, rcp, climate, paste0("run", run))
    filePath <- file.path(runDir, paste0("community-input-file-", year, ".csv"))
    
    if (file.exists(filePath)) {
      strData <- data.table::fread(filePath)
      
      strData$Run <- run
      strData$Time <- year
      
      df <- strData %>%
        filter((SpeciesName == "QUAL" & CohortAge >= 30) |
                 (SpeciesName == "CAAL27" & CohortAge >= 25) |
                 (SpeciesName == "CAGL8" & CohortAge >= 30) |
                 (SpeciesName == "LITU" & CohortAge >= 20) |
                 (SpeciesName == "OXAR" & CohortAge >= 40)) %>%
        group_by(MapCode, SpeciesName, CohortAge, Run, Time) %>%
        summarise(CohortBiomass = sum(CohortBiomass), .groups = "drop") %>%
        mutate(Weight = CohortAge / max(CohortAge, na.rm = TRUE)) %>%
        mutate(WeightedBiomass = CohortBiomass * Weight) %>%
        ungroup() %>%
        group_by(MapCode, Time, Run) %>%
        summarise(Score = sum(WeightedBiomass), .groups = "drop")
      
      df$Scenario <- scenario
      df$RCP <- rcp
      df$Climate <- climate
      
      strList[[length(strList) + 1]] <- df
    }
  }
}

strDf <- bind_rows(strList)

# Make a for loop to assign the fruit/nut tree ecosystem service tree to the individual
# cells within the raster. Check to make sure they all assigned properly
strRasters <- list()

for (run in runs) {
  for (timestep in timesteps) {
    rasterPath <- file.path(mainDir, scenario, rcp, climate, paste0("run", run), paste0("output-community-", timestep, ".img"))
    rasterRun <- rast(rasterPath) %>% flip()
    dfRun <- strDf %>%
      filter(Run == run, Time == timestep) %>%
      mutate(MapCode = as.numeric(MapCode))
    rasterCulture <- subst(rasterRun, from = dfRun$MapCode, to = dfRun$Score)
    
    # medianValue <- median(values(rasterRaster), na.rm = TRUE)
    # rasterH <- ifel(rasterH > 3, medianValue, rasterH)
    
    name <- paste(scenario, rcp, climate, paste0("run", run), timestep, sep = "_")
    strRasters[[name]] <- rasterCulture
  }
}

# Calculate the mean score for each cell across all runs for each raster
meanRasters <- list()
raster_names <- names(strRasters)
group_keys <- gsub("_run[0-9]+", "", raster_names) # Create groupings based on prescription, rcp, climate and timestep
unique_keys <- unique(group_keys)

for (key in unique_keys) {
  matched_rasters <- strRasters[which(group_keys == key)] # Grab the 5 runs for each model scenario
  raster_stack <- rast(matched_rasters)
  mean_raster <- app(raster_stack, mean, na.rm = TRUE)
  meanRasters[[key]] <- mean_raster
}

# Write the raster out
writeRaster(meanRasters[[key]], file = '/share/tcsi/lagoodal/R/Outputs/Cultural_Maps/', paste0(key, '.tif'),
            NAflag = 0, overwrite = TRUE, datatype = "FLT4S")
