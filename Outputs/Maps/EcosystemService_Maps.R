library(terra)
library(tidyterra)
library(tidyverse)
library(vegan)

# Read the index passed in from the job array
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
if (is.na(i)) stop("Error: Invalid value for i (LSB_JOBINDEX).")

# Set map names
mapNames <- c(
  "Resilience_RCP45_HSHV","Resilience_RCP45_HSLV","Resilience_RCP45_LSHV","Resilience_RCP45_LSLV",
  "Resilience_RCP85_HSHV","Resilience_RCP85_HSLV","Resilience_RCP85_LSHV","Resilience_RCP85_LSLV",
  "Resistance_RCP45_HSHV","Resistance_RCP45_HSLV","Resistance_RCP45_LSHV","Resistance_RCP45_LSLV",
  "Resistance_RCP85_HSHV","Resistance_RCP85_HSLV","Resistance_RCP85_LSHV","Resistance_RCP85_LSLV",
  "Transition_RCP45_HSHV","Transition_RCP45_HSLV","Transition_RCP45_LSHV","Transition_RCP45_LSLV",
  "Transition_RCP85_HSHV","Transition_RCP85_HSLV","Transition_RCP85_LSHV","Transition_RCP85_LSLV"
)

if (i < 1 || i > length(mapNames)) stop("Error: i is out of bounds for the mapNames list.")
mapName <- mapNames[i]
cat("Processing map:", mapName, "\n")

# Helper functions
get_rasters_from_subfolder <- function(folder) {
  rasters <- list.files(folder, pattern = "\\.tif$", full.names = TRUE)
  raster_list <- lapply(rasters, rast)
  names(raster_list) <- tools::file_path_sans_ext(basename(rasters))
  raster_list
}

normaliseRasters <- function(rasters){
  lapply(rasters, function(r){
    r <- app(r, function(x) ifelse(x == 0, NA, x))
    r_min <- global(r, "min", na.rm = TRUE)[[1]]
    r_max <- global(r, "max", na.rm = TRUE)[[1]]
    (r - r_min) / (r_max - r_min)
  })
}

makeTilesWithOverlap <- function(r, tile_width, tile_height, overlap = 0,
                                 filename_template = NULL, overwrite = FALSE) {
  ext_r <- ext(r)
  xmin <- ext_r[1]; xmax <- ext_r[2]
  ymin <- ext_r[3]; ymax <- ext_r[4]
  x_tiles <- ceiling((xmax - xmin) / tile_width)
  y_tiles <- ceiling((ymax - ymin) / tile_height)
  count <- 1
  tiles <- list()

  for (i in 0:(x_tiles - 1)) {
    for (j in 0:(y_tiles - 1)) {
      x_start <- max(xmin, xmin + i * tile_width - overlap)
      x_end   <- min(xmax, xmin + (i + 1) * tile_width + overlap)
      y_start <- max(ymin, ymin + j * tile_height - overlap)
      y_end   <- min(ymax, ymin + (j + 1) * tile_height + overlap)
      tile_extent <- ext(x_start, x_end, y_start, y_end)
      tile_rast <- crop(r, tile_extent)
      if (!is.null(filename_template)) {
        filename <- sprintf(filename_template, count)
        writeRaster(tile_rast, filename, overwrite = overwrite)
      }
      tiles[[count]] <- tile_rast
      count <- count + 1
    }
  }
  return(tiles)
}

# Load only relevant rasters for the selected map
mainDir <- "/share/tcsi/lagoodal/R/Data/ES_Rasters"
subfolders <- list.dirs(mainDir, recursive = FALSE)
all_subfolder_rasters <- lapply(subfolders, get_rasters_from_subfolder)

grouped_rasters <- lapply(all_subfolder_rasters, function(subfolder) {
  if (mapName %in% names(subfolder)) subfolder[[mapName]] else NULL
}) %>% compact()

if (length(grouped_rasters) == 0) {
  stop(paste("No rasters found for map:", mapName))
}

# Normalize and stack
normalized_group <- normaliseRasters(grouped_rasters)
currentRast <- tryCatch({
  rast(normalized_group)
}, error = function(e) {
  stop(paste("Failed to stack rasters for", mapName, ":", e$message))
})

# Clean storage and create tiles
storageDir <- file.path("/share/tcsi/lagoodal/R/Data/Storage", mapName)
dir.create(storageDir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(storageDir, "*.tif"))  # delete previous tiles

filename_template <- file.path(storageDir, paste0(mapName, "_tile_%03d.tif"))
makeTilesWithOverlap(currentRast, tile_width = 5000, tile_height = 5000, overlap = 300,
                     filename_template = filename_template, overwrite = TRUE)

# Process each tile
rasterFiles <- list.files(storageDir, pattern = "\\.tif$", full.names = TRUE)
rasterList <- list()

for (j in rasterFiles) {
  r <- rast(j)
  coords <- r %>% as.data.frame(xy = TRUE, na.rm = FALSE) %>% select(x, y)
  valMatrix <- values(r)

  if (all(is.nan(valMatrix))) {
    valMatrix[is.nan(valMatrix)] <- 0
  } else {
    for (k in seq_len(ncol(valMatrix))) {
      colVals <- valMatrix[, k]
      colMedian <- if (all(is.na(colVals))) 0 else median(colVals, na.rm = TRUE)
      colVals[is.na(colVals)] <- colMedian
      valMatrix[, k] <- colVals
    }
  }

  distMatrix <- vegdist(valMatrix, method = "bray")
  distMatrixDf <- as.data.frame(as.matrix(distMatrix))
  distMatrixDf <- cbind(coords, rowSums(distMatrixDf)) %>% rename("Sum" = "rowSums(distMatrixDf)")
  raster <- as_spatraster(distMatrixDf, crs = "epsg:26917")

  if (any(!is.na(values(r[[1]])))) {
    raster <- mask(raster, r[[1]])
  } else {
    warning(paste("Skipping mask for", j, ": mask raster has no values."))
  }

  rasterList[[j]] <- raster
}

# Mosaic all tiles together
mosaicRast <- mosaic(sprc(rasterList), fun = "mean")
mosaicRast <- focal(mosaicRast, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
mosaicRast <- mask(mosaicRast, currentRast[[1]])

# Write final output
outputPath <- file.path("/share/tcsi/lagoodal/R/Data/Finished", paste0(mapName, ".tif"))
dir.create(dirname(outputPath), recursive = TRUE, showWarnings = FALSE)
writeRaster(mosaicRast, filename = outputPath, NAflag = 0, datatype = "FLT4S", overwrite = TRUE)

cat("Finished processing", mapName, "\n")




###############################################
# args <- commandArgs(trailingOnly = TRUE)
# 
# # Ensure 'i' is correctly assigned from the command-line argument
# i <- as.integer(args[1])
# 
# print(i)
# 
# # Validate if 'i' is correctly parsed
# if (is.na(i)) {
#   stop("Error: Invalid value for i (LSB_JOBINDEX).")
# }
# 
# makeTilesWithOverlap <- function(r, tile_width, tile_height, overlap = 0,
#                                  filename_template = NULL, overwrite = FALSE) {
#   ext_r <- ext(r)
#   xmin <- ext_r[1]
#   xmax <- ext_r[2]
#   ymin <- ext_r[3]
#   ymax <- ext_r[4]
#   x_tiles <- ceiling((xmax - xmin) / tile_width)
#   y_tiles <- ceiling((ymax - ymin) / tile_height)
#   tiles <- list()
#   count <- 1
# 
#   for (i in 0:(x_tiles - 1)) {
#     for (j in 0:(y_tiles - 1)) {
#       x_start <- xmin + i * tile_width - overlap
#       x_end   <- xmin + (i + 1) * tile_width + overlap
#       y_start <- ymin + j * tile_height - overlap
#       y_end   <- ymin + (j + 1) * tile_height + overlap
# 
#       x_start <- max(x_start, xmin)
#       x_end   <- min(x_end, xmax)
#       y_start <- max(y_start, ymin)
#       y_end   <- min(y_end, ymax)
# 
#       tile_extent <- ext(x_start, x_end, y_start, y_end)
#       tile_rast <- crop(r, tile_extent)
#       if (!is.null(filename_template)) {
#         filename <- sprintf(filename_template, count)
#         writeRaster(tile_rast, filename, overwrite = overwrite)
#       }
#       tiles[[count]] <- tile_rast
#       count <- count + 1
#     }
#   }
#   return(tiles)
# }
# 
# normaliseRasters <- function(rasters){
#   lapply(rasters, function(r){
#     r <- app(r, function(x) ifelse(x == 0, NA, x))
#     r_min <- global(r, "min", na.rm = TRUE)[[1]]
#     r_max <- global(r, "max", na.rm = TRUE)[[1]]
#     (r - r_min) / (r_max - r_min)
#   })
# }
# 
# get_rasters_from_subfolder <- function(folder) {
#   rasters <- list.files(folder, pattern = "\\.tif$", full.names = TRUE)
#   raster_list <- lapply(rasters, rast)
#   names(raster_list) <- tools::file_path_sans_ext(basename(rasters))
#   return(raster_list)
# }
# 
# # --- Setup inputs and stack ---
# mainDir <- "/share/tcsi/lagoodal/R/Data/ES_Rasters"
# subfolders <- list.dirs(mainDir, recursive = FALSE)
# colNames <- c("x","y","AGB","ANPP","AWB","Carbon","H","NEEC","SOMTC")
# all_subfolder_rasters <- lapply(subfolders, get_rasters_from_subfolder)
# 
# #################################
# #cat("Found", length(subfolders), "subfolders:\n")
# #print(subfolders)
# 
# #cat("\nAll subfolder raster names:\n")
# #print(lapply(all_subfolder_rasters, names))
# #################################
# 
# all_raster_names <- unique(unlist(lapply(all_subfolder_rasters, names)))
# grouped_rasters <- lapply(all_raster_names, function(name) {
#   lapply(all_subfolder_rasters, function(subfolder) {
#     if (name %in% names(subfolder)) subfolder[[name]] else NULL
#   }) %>% compact()
# })
# names(grouped_rasters) <- all_raster_names
# 
# ##################################
# #cat("\nGrouped raster names:\n")
# #print(names(grouped_rasters))
# 
# #cat("Number of rasters per group:\n")
# #print(sapply(grouped_rasters, length))
# ##################################
# 
# normalized_groups <- lapply(grouped_rasters, normaliseRasters)
# 
# ##################################
# #cat("\nNormalized group lengths:\n")
# #print(sapply(normalized_groups, length))
# ##################################
# 
# #stacked_rasters <- lapply(normalized_groups, function(x) rast(x))
# stacked_rasters <- lapply(normalized_groups, function(x) {
#   if (length(x) == 0) return(NULL)
#   tryCatch(
#     rast(x),
#     error = function(e) {
#       message("Stack failed: ", e)
#       return(NULL)
#     }
#   )
# }) %>% compact()
# 
# 
# ##################################
# cat("\nFinal stacked raster count:\n")
# print(length(stacked_rasters))
# cat("Stacked raster names (if any):\n")
# print(names(stacked_rasters))
# #################################
# 
# mapNames <- c("Resilience_RCP45_HSHV","Resilience_RCP45_HSLV","Resilience_RCP45_LSHV","Resilience_RCP45_LSLV",
#               "Resilience_RCP85_HSHV","Resilience_RCP85_HSLV","Resilience_RCP85_LSHV","Resilience_RCP85_LSLV",
#               "Resistance_RCP45_HSHV","Resistance_RCP45_HSLV","Resistance_RCP45_LSHV","Resistance_RCP45_LSLV",
#               "Resistance_RCP85_HSHV","Resistance_RCP85_HSLV","Resistance_RCP85_LSHV","Resistance_RCP85_LSLV",
#               "Transition_RCP45_HSHV","Transition_RCP45_HSLV","Transition_RCP45_LSHV","Transition_RCP45_LSLV",
#               "Transition_RCP85_HSHV","Transition_RCP85_HSLV","Transition_RCP85_LSHV","Transition_RCP85_LSLV")
# 
# # Ensure that the index `i` is valid for the length of `stacked_rasters`
# #if (i < 1 || i > length(stacked_rasters)) {
# #  stop("Error: i is out of bounds for the stacked_rasters list.")
# #}
# 
# # Use the index (i) to pick the right raster stack
# #currentRast <- stacked_rasters[[i]]
# #mapName <- mapNames[i]
# ################################################################################################################
# 
# 
# # -- New: determine map name and load just what we need --
# if (i < 1 || i > length(mapNames)) {
#   stop("Error: i is out of bounds for the mapNames list.")
# }
# mapName <- mapNames[i]
# cat("Processing:", mapName, "\n")
# 
# mainDir <- "R/Data/ES_Rasters"
# subfolders <- list.dirs(mainDir, recursive = FALSE)
# all_subfolder_rasters <- lapply(subfolders, get_rasters_from_subfolder)
# 
# 
# 
# 
# ##############################
# cat("Looking in subfolders for map name:", mapName, "\n")
# cat("Subfolder names:\n")
# print(subfolders)
# 
# cat("Available raster names in each folder:\n")
# print(lapply(all_subfolder_rasters, names))
# ################################
# 
# # Grab only the rasters for this map
# grouped_rasters <- lapply(all_subfolder_rasters, function(subfolder) {
#   if (mapName %in% names(subfolder)) subfolder[[mapName]] else NULL
# }) %>% compact()
# 
# if (length(grouped_rasters) == 0) {
#   stop(paste("No rasters found for map:", mapName))
# }
# 
# normalized_group <- normaliseRasters(grouped_rasters)
# 
# currentRast <- tryCatch({
#   rast(normalized_group)
# }, error = function(e) {
#   stop(paste("Failed to stack rasters for", mapName, ":", e$message))
# })
# 
# 
# names(stacked_rasters) <- mapNames[1:length(stacked_rasters)]
# mapName <- mapNames[i]
# if (!(mapName %in% names(stacked_rasters))) {
# 	stop(paste("No stacked raster found for", mapNames))
# }
# currentRast <- stacked_rasters[[mapName]]
# 
# 
# rasterList <- list()
# 
# storageDir <- paste0("R/Data/Storage/", mapName)
# filesToDelete <- list.files(storageDir, pattern = "\\.tif$", full.names = TRUE)
# if (length(filesToDelete) > 0) file.remove(filesToDelete)
# 
# filename_template <- paste0(storageDir, "/", mapName, "_tile_%03d.tif")
# makeTilesWithOverlap(currentRast, tile_width = 5000, tile_height = 5000, overlap = 300,
#                      filename_template = filename_template, overwrite = TRUE)
# 
# rasterFiles <- list.files(storageDir, pattern = "\\.tif$", full.names = TRUE)
# 
# for (j in rasterFiles) {
#   r <- rast(j)
#   coords <- r %>% as.data.frame(xy = TRUE, na.rm = FALSE) %>% select(x, y)
#   valMatrix <- values(r)
# 
#   if (all(is.nan(valMatrix))) {
#     valMatrix[is.nan(valMatrix)] <- 0
#   } else {
#     for (k in seq_len(ncol(valMatrix))) {
#       colVals <- valMatrix[, k]
#       colMedian <- if (all(is.na(colVals))) 0 else median(colVals, na.rm = TRUE)
#       colVals[is.na(colVals)] <- colMedian
#       valMatrix[, k] <- colVals
#     }
#   }
# 
#   distMatrix <- vegdist(valMatrix, method = "bray")
#   distMatrixDf <- as.data.frame(as.matrix(distMatrix))
#   distMatrixDf <- cbind(coords, rowSums(distMatrixDf)) %>% rename("Sum" = "rowSums(distMatrixDf)")
#   raster <- as_spatraster(distMatrixDf, crs = "epsg:26917")
# 
#   if (any(!is.na(values(r[[1]])))) {
#     raster <- mask(raster, r[[1]])
#   } else {
#     warning(paste("Skipping mask for", j, ": mask raster has no values."))
#   }
# 
#   rasterList[[j]] <- raster
# }
# 
# mosaicRast <- mosaic(sprc(rasterList), fun = "mean")
# mosaicRast <- focal(mosaicRast, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
# mosaicRast <- mask(mosaicRast, currentRast[[1]])
# outputPath <- paste0("R/Outputs/Finished/", mapName, ".tif")
# writeRaster(mosaicRast, filename = outputPath, NAflag = 0, datatype = "FLT4S", overwrite = TRUE)


###############################################

# args <- commandArgs(trailingOnly = TRUE)
# i <- as.integer(args[1])
# 
# makeTilesWithOverlap <- function(r, tile_width, tile_height, overlap = 0, 
#                                  filename_template = NULL, overwrite = FALSE) {
#   ext_r <- ext(r)
#   xmin <- ext_r[1]
#   xmax <- ext_r[2]
#   ymin <- ext_r[3]
#   ymax <- ext_r[4]
#   x_tiles <- ceiling((xmax - xmin) / tile_width)
#   y_tiles <- ceiling((ymax - ymin) / tile_height)
#   tiles <- list()
#   count <- 1
#   
#   for (i in 0:(x_tiles - 1)) {
#     for (j in 0:(y_tiles - 1)) {
#       x_start <- xmin + i * tile_width - overlap
#       x_end   <- xmin + (i + 1) * tile_width + overlap
#       y_start <- ymin + j * tile_height - overlap
#       y_end   <- ymin + (j + 1) * tile_height + overlap
#       
#       x_start <- max(x_start, xmin)
#       x_end   <- min(x_end, xmax)
#       y_start <- max(y_start, ymin)
#       y_end   <- min(y_end, ymax)
#       
#       tile_extent <- ext(x_start, x_end, y_start, y_end)
#       tile_rast <- crop(r, tile_extent)
#       if (!is.null(filename_template)) {
#         filename <- sprintf(filename_template, count)
#         writeRaster(tile_rast, filename, overwrite = overwrite)
#       }
#       tiles[[count]] <- tile_rast
#       count <- count + 1
#     }
#   }
#   return(tiles)
# }
# 
# normaliseRasters <- function(rasters){
#   lapply(rasters, function(r){
#     r <- app(r, function(x) ifelse(x == 0, NA, x))
#     r_min <- global(r, "min", na.rm = TRUE)[[1]]
#     r_max <- global(r, "max", na.rm = TRUE)[[1]]
#     (r - r_min) / (r_max - r_min)
#   })
# }
# 
# get_rasters_from_subfolder <- function(folder) {
#   rasters <- list.files(folder, pattern = "\\.tif$", full.names = TRUE)
#   raster_list <- lapply(rasters, rast)
#   names(raster_list) <- tools::file_path_sans_ext(basename(rasters))
#   return(raster_list)
# }
# 
# # --- Setup inputs and stack ---
# mainDir <- "R/Data/ES_Rasters"
# subfolders <- list.dirs(mainDir, recursive = FALSE)
# colNames <- c("x","y","AGB","ANPP","AWB","Carbon","H","NEEC","SOMTC")
# all_subfolder_rasters <- lapply(subfolders, get_rasters_from_subfolder)
# all_raster_names <- unique(unlist(lapply(all_subfolder_rasters, names)))
# grouped_rasters <- lapply(all_raster_names, function(name) {
#   lapply(all_subfolder_rasters, function(subfolder) {
#     if (name %in% names(subfolder)) subfolder[[name]] else NULL
#   }) %>% compact()
# })
# names(grouped_rasters) <- all_raster_names
# normalized_groups <- lapply(grouped_rasters, normaliseRasters)
# stacked_rasters <- lapply(normalized_groups, function(x) rast(x))
# 
# mapNames <- c("Resilience_RCP45_HSHV","Resilience_RCP45_HSLV","Resilience_RCP45_LSHV","Resilience_RCP45_LSLV",
#               "Resilience_RCP85_HSHV","Resilience_RCP85_HSLV","Resilience_RCP85_LSHV","Resilience_RCP85_LSLV",
#               "Resistance_RCP45_HSHV","Resistance_RCP45_HSLV","Resistance_RCP45_LSHV","Resistance_RCP45_LSLV",
#               "Resistance_RCP85_HSHV","Resistance_RCP85_HSLV","Resistance_RCP85_LSHV","Resistance_RCP85_LSLV",
#               "Transition_RCP45_HSHV","Transition_RCP45_HSLV","Transition_RCP45_LSHV","Transition_RCP45_LSLV",
#               "Transition_RCP85_HSHV","Transition_RCP85_HSLV","Transition_RCP85_LSHV","Transition_RCP85_LSLV")
# 
# currentRast <- stacked_rasters[[i]]
# mapName <- mapNames[i]
# rasterList <- list()
# 
# storageDir <- paste0("R/Data/Storage/", mapName)
# filesToDelete <- list.files(storageDir, pattern = "\\.tif$", full.names = TRUE)
# if (length(filesToDelete) > 0) file.remove(filesToDelete)
# 
# filename_template <- paste0(storageDir, "/", mapName, "_tile_%03d.tif")
# makeTilesWithOverlap(currentRast, tile_width = 5000, tile_height = 5000, overlap = 300,
#                      filename_template = filename_template, overwrite = TRUE)
# 
# rasterFiles <- list.files(storageDir, pattern = "\\.tif$", full.names = TRUE)
# 
# for (j in rasterFiles) {
#   r <- rast(j)
#   coords <- r %>% as.data.frame(xy = TRUE, na.rm = FALSE) %>% select(x, y)
#   valMatrix <- values(r)
#   
#   if (all(is.nan(valMatrix))) {
#     valMatrix[is.nan(valMatrix)] <- 0
#   } else {
#     for (k in seq_len(ncol(valMatrix))) {
#       colVals <- valMatrix[, k]
#       colMedian <- if (all(is.na(colVals))) 0 else median(colVals, na.rm = TRUE)
#       colVals[is.na(colVals)] <- colMedian
#       valMatrix[, k] <- colVals
#     }
#   }
#   
#   distMatrix <- vegdist(valMatrix, method = "bray")
#   distMatrixDf <- as.data.frame(as.matrix(distMatrix))
#   distMatrixDf <- cbind(coords, rowSums(distMatrixDf)) %>% rename("Sum" = "rowSums(distMatrixDf)")
#   raster <- as_spatraster(distMatrixDf, crs = "epsg:26917")
#   
#   if (any(!is.na(values(r[[1]])))) {
#     raster <- mask(raster, r[[1]])
#   } else {
#     warning(paste("Skipping mask for", j, ": mask raster has no values."))
#   }
#   
#   rasterList[[j]] <- raster
# }
# 
# mosaicRast <- mosaic(sprc(rasterList), fun = "mean")
# mosaicRast <- focal(mosaicRast, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
# mosaicRast <- mask(mosaicRast, currentRast[[1]])
# outputPath <- paste0("R/Outputs/Finished/", mapName, ".tif")
# writeRaster(mosaicRast, filename = outputPath, NAflag = 0, datatype = "FLT4S", overwrite = TRUE)


####################################################

# args <- commandArgs(trailingOnly = TRUE)
# task_id <- as.integer(args[1])  # from SLURM_ARRAY_TASK_ID
# 
# mainDir <- "/share/tcsi/lagoodal/R/Data/ES_Rasters/"
# storageDir <- "/share/tcsi/lagoodal/R/Data/Storage/"
# outputDir <- "/share/tcsi/lagoodal/R/Data/Finished/"
# 
# mapNames <- c("Resilience_RCP45_HSHV", "Resilience_RCP45_HSLV", "Resilience_RCP45_LSHV", "Resilience_RCP45_LSLV",
#               "Resilience_RCP85_HSHV", "Resilience_RCP85_HSLV", "Resilience_RCP85_LSHV", "Resilience_RCP85_LSLV",
#               "Resistance_RCP45_HSHV", "Resistance_RCP45_HSLV", "Resistance_RCP45_LSHV", "Resistance_RCP45_LSLV",
#               "Resistance_RCP85_HSHV", "Resistance_RCP85_HSLV", "Resistance_RCP85_LSHV", "Resistance_RCP85_LSLV",
#               "Transition_RCP45_HSHV", "Transition_RCP45_HSLV", "Transition_RCP45_LSHV", "Transition_RCP45_LSLV",
#               "Transition_RCP85_HSHV", "Transition_RCP85_HSLV", "Transition_RCP85_LSHV", "Transition_RCP85_LSLV")
# 
# # Load the raster stacks
# subfolders <- list.dirs(mainDir, recursive = FALSE)
# get_rasters <- function(folder) {
#   rasters <- list.files(folder, pattern = "\\.tif$", full.names = TRUE)
#   raster_list <- lapply(rasters, rast)
#   names(raster_list) <- tools::file_path_sans_ext(basename(rasters))
#   return(raster_list)
# }
# all_rasters <- lapply(subfolders, get_rasters)
# all_names <- unique(unlist(lapply(all_rasters, names)))
# grouped <- lapply(all_names, function(name) {
#   lapply(all_rasters, function(sub) if (name %in% names(sub)) sub[[name]] else NULL) %>% compact()
# })
# names(grouped) <- all_names
# 
# normalize <- function(rasters) {
#   lapply(rasters, function(r) {
#     r <- app(r, function(x) ifelse(x == 0, NA, x))
#     r_min <- global(r, "min", na.rm = TRUE)[[1]]
#     r_max <- global(r, "max", na.rm = TRUE)[[1]]
#     (r - r_min) / (r_max - r_min)
#   })
# }
# 
# normalized_groups <- lapply(grouped, normalize)
# stacked_rasters <- lapply(normalized_groups, function(x) rast(x))
# 
# # Function to create overlapping tiles
# makeTilesWithOverlap <- function(r, tile_width, tile_height, overlap = 0, 
#                                  filename_template = NULL, overwrite = FALSE) {
#   ext_r <- ext(r)
#   xmin <- ext_r[1]; xmax <- ext_r[2]; ymin <- ext_r[3]; ymax <- ext_r[4]
#   x_tiles <- ceiling((xmax - xmin) / tile_width)
#   y_tiles <- ceiling((ymax - ymin) / tile_height)
#   tiles <- list(); count <- 1
#   for (i in 0:(x_tiles - 1)) {
#     for (j in 0:(y_tiles - 1)) {
#       x_start <- max(xmin + i * tile_width - overlap, xmin)
#       x_end   <- min(xmin + (i + 1) * tile_width + overlap, xmax)
#       y_start <- max(ymin + j * tile_height - overlap, ymin)
#       y_end   <- min(ymin + (j + 1) * tile_height + overlap, ymax)
#       tile_extent <- ext(x_start, x_end, y_start, y_end)
#       tile_rast <- crop(r, tile_extent)
#       if (!is.null(filename_template)) {
#         filename <- sprintf(filename_template, count)
#         writeRaster(tile_rast, filename, overwrite = overwrite)
#       }
#       tiles[[count]] <- tile_rast
#       count <- count + 1
#     }
#   }
#   return(tiles)
# }
# 
# # Process current task
# currentRast <- stacked_rasters[[task_id]]
# mapName <- mapNames[task_id]
# 
# # Unique temp folder
# temp_folder <- file.path(storageDir, mapName)
# dir.create(temp_folder, showWarnings = FALSE, recursive = TRUE)
# 
# filename_template <- file.path(temp_folder, paste0(mapName, "_tile_%03d.tif"))
# makeTilesWithOverlap(currentRast, tile_width = 5000, tile_height = 5000, overlap = 300,
#                      filename_template = filename_template, overwrite = TRUE)
# 
# rasterFiles <- list.files(temp_folder, pattern = "\\.tif$", full.names = TRUE)
# rasterList <- list()
# 
# for (j in rasterFiles) {
#   r <- rast(j)
#   coords <- r %>% as.data.frame(xy = TRUE, na.rm = FALSE) %>% select(x, y)
#   valMatrix <- values(r)[, -6]
#   if (all(is.nan(valMatrix))) {
#     valMatrix[is.nan(valMatrix)] <- 0
#   } else {
#     for (k in seq_len(ncol(valMatrix))) {  
#       colVals <- valMatrix[, k]
#       colMedian <- if (all(is.na(colVals))) 0 else median(colVals, na.rm = TRUE)
#       colVals[is.na(colVals)] <- colMedian
#       valMatrix[, k] <- colVals
#     }
#   }
#   distMatrix <- vegdist(valMatrix, method = "bray")
#   distMatrixDf <- as.data.frame(as.matrix(distMatrix))
#   distMatrixDf <- cbind(coords, rowSums(distMatrixDf)) %>% rename("Sum" = "rowSums(distMatrixDf)")
#   raster <- as_spatraster(distMatrixDf, crs = "epsg:26917")
#   raster <- mask(raster, r[[1]])
#   rasterList[[j]] <- raster
# }
# 
# mosaicRast <- mosaic(sprc(rasterList), fun = "mean")
# mosaicRast <- focal(mosaicRast, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
# mosaicRast <- mask(mosaicRast, currentRast)
# 
# output_path <- file.path(outputDir, paste0(mapName, ".tif"))
# writeRaster(mosaicRast, filename = output_path, NAflag = 0, datatype = "FLT4S", overwrite = TRUE)
# 
# # Optional cleanup
# unlink(temp_folder, recursive = TRUE)


###################################

# library(terra)
# library(dplyr)
# library(vegan)
# library(purrr)
# library(doParallel)
# library(foreach)
# 
# 
# args <- commandArgs(trailingOnly = TRUE)
# task_id <- as.integer(args[1])  # from SLURM_ARRAY_TASK_ID
# 
# mainDir <- "/share/tcsi/lagoodal/R/Data/ES_Rasters/"
# storageDir <- "/share/tcsi/lagoodal/R/Data/Storage/"
# outputDir <- "/share/tcsi/lagoodal/R/Data/Finished/"
# 
# mapNames <- c("Resilience_RCP45_HSHV", "Resilience_RCP45_HSLV", "Resilience_RCP45_LSHV", "Resilience_RCP45_LSLV",
#               "Resilience_RCP85_HSHV", "Resilience_RCP85_HSLV", "Resilience_RCP85_LSHV", "Resilience_RCP85_LSLV",
#               "Resistance_RCP45_HSHV", "Resistance_RCP45_HSLV", "Resistance_RCP45_LSHV", "Resistance_RCP45_LSLV",
#               "Resistance_RCP85_HSHV", "Resistance_RCP85_HSLV", "Resistance_RCP85_LSHV", "Resistance_RCP85_LSLV",
#               "Transition_RCP45_HSHV", "Transition_RCP45_HSLV", "Transition_RCP45_LSHV", "Transition_RCP45_LSLV",
#               "Transition_RCP85_HSHV", "Transition_RCP85_HSLV", "Transition_RCP85_LSHV", "Transition_RCP85_LSLV")
# 
# # Load and process rasters only once
# 
# 
# # Load the raster stacks
# subfolders <- list.dirs(mainDir, recursive = FALSE)
# get_rasters <- function(folder) {
#   rasters <- list.files(folder, pattern = "\\.tif$", full.names = TRUE)
#   raster_list <- lapply(rasters, rast)
#   names(raster_list) <- tools::file_path_sans_ext(basename(rasters))
#   return(raster_list)
# }
# all_rasters <- lapply(subfolders, get_rasters)
# all_names <- unique(unlist(lapply(all_rasters, names)))
# grouped <- lapply(all_names, function(name) {
#   lapply(all_rasters, function(sub) if (name %in% names(sub)) sub[[name]] else NULL) %>% compact()
# })
# names(grouped) <- all_names
# 
# normalize <- function(rasters) {
#   lapply(rasters, function(r) {
#     r <- app(r, function(x) ifelse(x == 0, NA, x))
#     r_min <- global(r, "min", na.rm = TRUE)[[1]]
#     r_max <- global(r, "max", na.rm = TRUE)[[1]]
#     (r - r_min) / (r_max - r_min)
#   })
# }
# 
# normalized_groups <- lapply(grouped, normalize)
# stacked_rasters <- lapply(normalized_groups, function(x) rast(x))
# 
# # Function to create overlapping tiles
# makeTilesWithOverlap <- function(r, tile_width, tile_height, overlap = 0, 
#                                  filename_template = NULL, overwrite = FALSE) {
#   ext_r <- ext(r)
#   xmin <- ext_r[1]; xmax <- ext_r[2]; ymin <- ext_r[3]; ymax <- ext_r[4]
#   x_tiles <- ceiling((xmax - xmin) / tile_width)
#   y_tiles <- ceiling((ymax - ymin) / tile_height)
#   tiles <- list(); count <- 1
#   for (i in 0:(x_tiles - 1)) {
#     for (j in 0:(y_tiles - 1)) {
#       x_start <- max(xmin + i * tile_width - overlap, xmin)
#       x_end   <- min(xmin + (i + 1) * tile_width + overlap, xmax)
#       y_start <- max(ymin + j * tile_height - overlap, ymin)
#       y_end   <- min(ymin + (j + 1) * tile_height + overlap, ymax)
#       tile_extent <- ext(x_start, x_end, y_start, y_end)
#       tile_rast <- crop(r, tile_extent)
#       if (!is.null(filename_template)) {
#         filename <- sprintf(filename_template, count)
#         writeRaster(tile_rast, filename, overwrite = overwrite)
#       }
#       tiles[[count]] <- tile_rast
#       count <- count + 1
#     }
#   }
#   return(tiles)
# }
# 
# # Process current task
# currentRast <- stacked_rasters[[task_id]]
# mapName <- mapNames[task_id]
# 
# # Unique temp folder
# temp_folder <- file.path(storageDir, mapName)
# dir.create(temp_folder, showWarnings = FALSE, recursive = TRUE)
# 
# filename_template <- file.path(temp_folder, paste0(mapName, "_tile_%03d.tif"))
# makeTilesWithOverlap(currentRast, tile_width = 5000, tile_height = 5000, overlap = 300,
#                      filename_template = filename_template, overwrite = TRUE)
# 
# rasterFiles <- list.files(temp_folder, pattern = "\\.tif$", full.names = TRUE)
# rasterList <- list()
# 
# for (j in rasterFiles) {
#   r <- rast(j)
#   coords <- r %>% as.data.frame(xy = TRUE, na.rm = FALSE) %>% select(x, y)
#   valMatrix <- values(r)[, -6]
#   if (all(is.nan(valMatrix))) {
#     valMatrix[is.nan(valMatrix)] <- 0
#   } else {
#     for (k in seq_len(ncol(valMatrix))) {
#       colVals <- valMatrix[, k]
#       colMedian <- if (all(is.na(colVals))) 0 else median(colVals, na.rm = TRUE)
#       colVals[is.na(colVals)] <- colMedian
#       valMatrix[, k] <- colVals
#     }
#   }
#   distMatrix <- vegdist(valMatrix, method = "bray")
#   distMatrixDf <- as.data.frame(as.matrix(distMatrix))
#   distMatrixDf <- cbind(coords, rowSums(distMatrixDf)) %>% rename("Sum" = "rowSums(distMatrixDf)")
#   raster <- as_spatraster(distMatrixDf, crs = "epsg:26917")
#   raster <- mask(raster, r[[1]])
#   rasterList[[j]] <- raster
#   mosaicRast <- mosaic(sprc(rasterList), fun = "mean")
#   mosaicRast <- focal(mosaicRast, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
#   mosaicRast <- mask(mosaicRast, currentRast)
# }
# 
# output_path <- file.path(outputDir, paste0(mapName, ".tif"))
# writeRaster(mosaicRast, filename = output_path, NAflag = 0, datatype = "FLT4S", overwrite = TRUE)
# 
# # Optional cleanup
# unlink(temp_folder, recursive = TRUE)
