library(tidyverse)
library(tidyterra)
library(terra)
library(furrr)
library(future)
library(sf)

# Set main directory and common variables
mainDir <- "/share/tcsi/lagoodal/R/Maps/"
colNames <- c("x","y","AGB","ANPP","AWB","Carbon","H","Harvest","NEEC","SOMTC")
subfolders <- list.dirs(mainDir, recursive = FALSE)

get_rasters_from_subfolder <- function(folder) {
  rasters <- list.files(folder, pattern = "\\.tif$", full.names = TRUE)
  raster_list <- lapply(rasters, rast)
  names(raster_list) <- tools::file_path_sans_ext(basename(rasters))
  return(raster_list)
}

all_subfolder_rasters <- lapply(subfolders, get_rasters_from_subfolder)
all_raster_names <- unique(unlist(lapply(all_subfolder_rasters, names)))

grouped_rasters <- lapply(all_raster_names, function(name) {
  lapply(all_subfolder_rasters, function(subfolder) {
    if (name %in% names(subfolder)) subfolder[[name]] else NULL
  }) %>% compact()
})
names(grouped_rasters) <- all_raster_names

normaliseRasters <- function(rasters){
  lapply(rasters, function(r){
    r <- app(r, function(x) ifelse(x == 0, NA, x))
    r_min <- global(r, "min", na.rm = TRUE)[[1]]
    r_max <- global(r, "max", na.rm = TRUE)[[1]]
    (r - r_min) / (r_max - r_min)
  })
}

normalized_groups <- lapply(grouped_rasters, normaliseRasters)
stacked_rasters <- lapply(normalized_groups, function(x) {rast(x)})
df_list <- lapply(stacked_rasters, function(raster_obj) {as.data.frame(raster_obj, xy = TRUE)})
df_list <- lapply(df_list, function(x){
  colnames(x) <- colNames
  return(x)})
df_list <- lapply(df_list, function(x){x %>% mutate(across(-c(x,y), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))})

mapNames <- c(
  "Resilience_RCP45_HSHV","Resilience_RCP45_HSLV","Resilience_RCP45_LSHV","Resilience_RCP45_LSLV",
  "Resilience_RCP85_HSHV","Resilience_RCP85_HSLV","Resilience_RCP85_HSHV","Resilience_RCP85_LSLV",
  "Resistance_RCP45_HSHV","Resistance_RCP45_HSLV","Resistance_RCP45_LSHV","Resistance_RCP45_LSLV",
  "Resistance_RCP85_HSHV","Resistance_RCP85_HSLV","Resistance_RCP85_HSHV","Resistance_RCP85_LSLV",
  "Transition_RCP45_HSHV","Transition_RCP45_HSLV","Transition_RCP45_LSHV","Transition_RCP45_LSLV",
  "Transition_RCP85_HSHV","Transition_RCP85_HSLV","Transition_RCP85_HSHV","Transition_RCP85_LSLV"
)

# Set up parallel plan (this is global and should be done before using furrr)
plan(multisession, workers = 24)

# Define function that processes chunks with block-to-block distance computation
chunked_bray_row_sums <- function(data, chunk_size = 10000) {
  n <- nrow(data)
  feature_data <- data[, -c(1, 2)]
  coords <- data[, c("x", "y")]
  
  chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  
  rowsum_list <- future_map(chunks, function(indexes) {
    chunk <- feature_data[indexes, , drop = FALSE]
    # Compute full Bray-Curtis distance from chunk to ALL rows
    dist_matrix <- ecodist::distance(rbind(chunk, feature_data), method = "bray-curtis")
    dist_matrix <- as.matrix(dist_matrix)
    
    # Only take distances from chunk rows to full data (exclude self-pairings)
    chunk_rows <- seq_len(nrow(chunk))
    full_rows <- (nrow(chunk) + 1):(nrow(chunk) + nrow(feature_data))
    sub_dist <- dist_matrix[chunk_rows, full_rows, drop = FALSE]
    
    row_sums <- rowSums(sub_dist, na.rm = TRUE)
    tibble(index = indexes, dist_sum = row_sums)
  })
  
  full_sums <- bind_rows(rowsum_list)
  full_sums <- full_sums[order(full_sums$index), ]
  final_df <- bind_cols(coords, dist_sum = full_sums$dist_sum)
  return(final_df)
}

outputDir <- "/share/tcsi/lagoodal/R/Outputs"
finalMaps <- list()
for (i in seq_along(df_list)) {
  message("Processing map ", i, " of ", length(df_list), ": ", mapNames[i])
  df <- df_list[[i]]
  
  # Use the new chunked function
  dist_df <- chunked_bray_row_sums(df, chunk_size = 10000)
  
  # Convert the results to a spatial raster
  dist_rast <- as_spatraster(dist_df, crs = "epsg:26917")
  # Apply smoothing via a focal window
  final_map <- focal(dist_rast, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)
  # Mask the raster to a normalized forest area (or similar)
  final_map <- mask(final_map, normalized_groups[[1]][[1]])
  names(final_map) <- mapNames[i]
  finalMaps[[i]] <- final_map
  
  # Write the final map to disk
  file_path <- file.path(outputDir, paste0(mapNames[i], ".tif"))
  writeRaster(final_map, filename = file_path, overwrite = TRUE, NAflag = 0, datatype = "FLT4S")
  gc()
}



# ###################################
# mainDir <- "/share/tcsi/lagoodal/R/Maps/"
# colNames <- c("x","y","AGB","ANPP","AWB","Carbon","H","Harvest","NEEC","SOMTC")
# subfolders <- list.dirs(mainDir, recursive = FALSE)
# get_rasters_from_subfolder <- function(folder) {
#   rasters <- list.files(folder, pattern = "\\.tif$", full.names = TRUE)
#   raster_list <- lapply(rasters, rast)
#   names(raster_list) <- tools::file_path_sans_ext(basename(rasters))
#   return(raster_list)
# }
# all_subfolder_rasters <- lapply(subfolders, get_rasters_from_subfolder)
# all_raster_names <- unique(unlist(lapply(all_subfolder_rasters, names)))
# grouped_rasters <- lapply(all_raster_names, function(name) {
#   lapply(all_subfolder_rasters, function(subfolder) {
#     if (name %in% names(subfolder)) subfolder[[name]] else NULL
#   }) %>% compact()
# })
# names(grouped_rasters) <- all_raster_names
# # agbList <- c(grouped_rasters[[1]][[1]],grouped_rasters[[2]][[1]],grouped_rasters[[3]][[1]],grouped_rasters[[4]][[1]],grouped_rasters[[5]][[1]],grouped_rasters[[6]][[1]],grouped_rasters[[7]][[1]],grouped_rasters[[8]][[1]])
# normaliseRasters <- function(rasters){
#   lapply(rasters, function(r){
#     r <- app(r, function(x) ifelse(x == 0, NA, x))
#     r_min <- global(r, "min", na.rm = TRUE)[[1]]
#     r_max <- global(r, "max", na.rm = TRUE)[[1]]
#     (r - r_min) / (r_max - r_min)
#   })
# }
# 
# normalized_groups <- lapply(grouped_rasters, normaliseRasters)
# stacked_rasters <- lapply(normalized_groups, function(x) {rast(x)})
# df_list <- lapply(stacked_rasters, function(raster_obj) {as.data.frame(raster_obj, xy = TRUE)})
# df_list <- lapply(df_list, function(x){
#   colnames(x) <- colNames
#   return(x)})
# df_list <- lapply(df_list, function(x){x %>% mutate(across(-c(x,y), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))})
# 
# 
# mapNames <- c("Resilience_RCP45_HSHV","Resilience_RCP45_HSLV","Resilience_RCP45_LSHV","Resilience_RCP45_LSLV",
#               "Resilience_RCP85_HSHV","Resilience_RCP85_HSLV","Resilience_RCP85_LSHV","Resilience_RCP85_LSLV",
#               "Resistance_RCP45_HSHV","Resistance_RCP45_HSLV","Resistance_RCP45_LSHV","Resistance_RCP45_LSLV",
#               "Resistance_RCP85_HSHV","Resistance_RCP85_HSLV","Resistance_RCP85_LSHV","Resistance_RCP85_LSLV",
#               "Transition_RCP45_HSHV","Transition_RCP45_HSLV","Transition_RCP45_LSHV","Transition_RCP45_LSLV",
#               "Transition_RCP85_HSHV","Transition_RCP85_HSLV","Transition_RCP85_LSHV","Transition_RCP85_LSLV")
# 
# 
# plan(multisession, workers = 24)
# 
# chunked_bray_row_sums <- function(data, chunk_size = 10000) {
#   n <- nrow(data)
#   chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
#   feature_data <- data[, -c(1, 2)]  # Drop x and y
#   
#   # Process each chunk and compute distances to ALL rows
#   rowsum_list <- future_map(chunks, function(indexes) {
#     chunk <- feature_data[indexes, , drop = FALSE]
#     # Compute Bray-Curtis distance of each row in chunk to all rows
#     mat <- as.matrix(chunk)
#     all_mat <- as.matrix(feature_data)
#     chunk_dists <- vegan::vegdist(rbind(mat, all_mat), method = "bray")
#     # Extract only distances between chunk and full (ignore self)
#     chunk_rows <- 1:nrow(mat)
#     dist_matrix <- as.matrix(chunk_dists)[chunk_rows, -(chunk_rows)]
#     row_sums <- rowSums(dist_matrix, na.rm = TRUE)
#     tibble(index = indexes, dist_sum = row_sums)
#   })
#   
#   # Stitch the rows together
#   full_sums <- bind_rows(rowsum_list)
#   full_sums <- full_sums[order(full_sums$index), ]  # Restore original order
#   
#   coords <- data[, c("x", "y")]
#   final_df <- bind_cols(coords, dist_sum = full_sums$dist_sum)
#   return(final_df)
# }
# 
# # chunk_distance_matrix_parallel <- function(data, chunk_size = 10000, overlap = 500) {
# #   n <- nrow(data)
# #   starts <- seq(1, n, by = chunk_size)
# #   
# #   chunks <- future_map(starts, function(s) {
# #     e <- min(s + chunk_size + overlap - 1, n)
# #     chunk <- data[s:e, ]
# #     dist_chunk <- ecodist::distance(chunk[,-c(1,2)], method = "bray-curtis")
# #     dist_chunk <- rowSums(as.matrix(dist_chunk), na.rm = TRUE)
# #     chunk_coords <- chunk[, c("x", "y")]
# #     chunk_df <- cbind(chunk_coords, dist_sum = dist_chunk)
# #     return(chunk_df)
# #   })
# #   
# #   final_df <- bind_rows(chunks)
# #   final_df <- final_df[!duplicated(final_df[, c("x", "y")]), ]
# #   return(final_df)
# # }
# 
# outputDir <- "/share/tcsi/lagoodal/R/Outputs"
# for (i in seq_along(df_list)) {
#   message("Processing map ", i, " of ", length(df_list), ": ", mapNames[i])
#   df <- df_list[[i]]
#   dist_df <- chunked_bray_row_sums(df, chunk_size = 10000)
#   dist_rast <- as_spatraster(dist_df, crs = "epsg:26917")
#   final_map <- focal(dist_rast, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)
#   final_map <- mask(final_map, normalized_groups[[1]][[1]])
#   names(final_map) <- mapNames[i]
#   finalMaps[[i]] <- final_map
#   file_path <- file.path(outputDir, paste0(mapNames[i], ".tif"))
#   writeRaster(final_map, filename = file_path, overwrite = TRUE, NAflag = 0, datatype = "FLT4S")
#   gc()
# }
# 
# 
# # finalMaps <- list()
# # for (i in seq_along(df_list)) {
# #   x <- df_list[[i]]
# #   dist_mat <- ecodist::distance(x[,-c(1,2)], method = "bray-curtis")
# #   dist_df <- as.data.frame(as.matrix(dist_mat))
# #   dist_sum <- rowSums(dist_df, na.rm = TRUE)
# #   coords <- x[,c("x","y")]
# #   final_df <- cbind(coords, dist_sum)
# #   dist_rast <- as_spatraster(final_df, crs = "epsg:26917")
# #   final_map <- focal(dist_rast, w = matrix(1,5,5), fun = mean, na.rm = TRUE)
# #   final_map <- mask(final_map, normalized_groups[[1]][[1]])
# #   names(final_map) <- mapNames[i]
# #   finalMaps[[i]] <- final_map
# #   gc()
# #   file_path <- file.path(outputDir, paste0(mapNames[i], ".tif"))
# #   writeRaster(final_map, filename = file_path, overwrite = TRUE, NAflag = 0, datatype = "FLT4S")
# }