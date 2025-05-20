# Load libraries
library(tidyverse)
library(data.table)
library(terra)
library(sf)

# STRUCTURAL DIVERSITY
# Read in data
file_path <- "/share/tcsi/lagoodal/R/Data/combinedData.csv"
combinedData <- fread(file_path)

# Calculate H for structural diversity
df <- combinedData %>%
  group_by(MapCode, Time, Run) %>%
  mutate(TotalBiomass = sum(CohortBiomass)) %>%
  group_by(MapCode, CohortAge, Time, Run) %>%
  summarise(Pi = sum(CohortBiomass) / unique(TotalBiomass)) %>%
  ungroup() %>%
  group_by(Time, Run) %>%
  summarise(H = -sum(Pi * log(Pi), na.rm = TRUE))

fwrite(df, file = "/share/tcsi/lagoodal/R/Data/Structural_Diversity.csv")