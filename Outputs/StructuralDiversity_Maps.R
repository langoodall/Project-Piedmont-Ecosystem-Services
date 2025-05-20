library(tidyverse)
library(terra)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
job_index <- as.numeric(args[1])

# Fixed inputs
mainDir <- "/share/tcsi/lagoodal/R/MAPS"
scenarios <- c("Resistance", "Resilience", "Transition")
rcps <- c("RCP45", "RCP85")
climates <- c("HSHV", "HSLV", "LSHV", "LSLV")
runs <- paste0("run", 1:5)
timesteps <- seq(10, 80, 10)

# Determine scenario/rcp/climate combo
combo_df <- expand.grid(Scenario = scenarios, RCP = rcps, Climate = climates)
this_combo <- combo_df[job_index, ]
scenario <- as.character(this_combo$Scenario)
rcp <- as.character(this_combo$RCP)
climate <- as.character(this_combo$Climate)

cat("Processing combo:", scenario, rcp, climate, "\n")

# STEP 1: Read and summarize community-input CSVs
strList <- list()

for (run in runs) {
  for (year in timesteps) {
    runDir <- file.path(mainDir, scenario, rcp, climate, paste0("run", run))
    filePath <- file.path(runDir, paste0("community-input-file-", year, ".csv"))
    
    if (file.exists(filePath)) {
      strData <- data.table::fread(filePath)
      
      df <- strData %>%
        mutate(Time = year, Run = run) %>%
        group_by(MapCode, Run, Time) %>%
        summarise(
          H = -sum((table(CohortAge) / length(CohortAge)) * log(table(CohortAge) / length(CohortAge))),
          .groups = "drop"
        )
      
      df$Scenario <- scenario
      df$RCP <- rcp
      df$Climate <- climate
      
      strList[[length(strList) + 1]] <- df
    }
  }
}

strDf <- bind_rows(strList)

# STEP 2: Apply H values to rasters
strRasters <- list()

for (run in runs) {
  for (timestep in timesteps) {
    
    rasterPath <- file.path(mainDir, scenario, rcp, climate, paste0("run", run), paste0("output-community-", timestep, ".img"))
    
    rasterRun <- rast(rasterPath) %>% flip()
    
    dfRun <- strDf %>%
      filter(Run == run, Time == timestep) %>%
      mutate(MapCode = as.numeric(MapCode))
    
    rasterH <- subst(rasterRun, from = dfRun$MapCode, to = dfRun$H)
    
    medianValue <- median(values(rasterH), na.rm = TRUE)
    rasterH <- ifel(rasterH > 3, medianValue, rasterH)
    
    name <- paste(scenario, rcp, climate, paste0("run", run), timestep, sep = "_")
    strRasters[[name]] <- rasterH
  }
}

# STEP 3: Calculate cellwise mean across runs
meanRasters <- list()
raster_names <- names(strRasters)
group_keys <- gsub("_run[0-9]+", "", raster_names)
unique_keys <- unique(group_keys)

for (key in unique_keys) {
  matched_rasters <- strRasters[which(group_keys == key)]
  raster_stack <- rast(matched_rasters)
  mean_raster <- app(raster_stack, mean, na.rm = TRUE)
  meanRasters[[key]] <- mean_raster
}

# STEP 4: Save output
writeRaster(meanRasters[[key]], file = '/share/tcsi/lagoodal/R/Outputs/Structural_Diversity_Maps/', paste0(key, '.tif'),
            NAflag = 0, overwrite = TRUE, datatype = "FLT4S")









