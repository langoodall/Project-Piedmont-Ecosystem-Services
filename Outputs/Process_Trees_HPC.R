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

# NUT TREES
df <- combinedData %>%
  filter((SpeciesName == "CAAL27" & CohortAge >= 25) |
           (SpeciesName == "CAGL8" & CohortAge >= 30) |
           (SpeciesName == "QUAL" & CohortAge >= 30) |
           (SpeciesName == "QUCO2" & CohortAge >= 30) |
           (SpeciesName == "QUFA" & CohortAge >= 25) |
           (SpeciesName == "QULA3" & CohortAge >= 35) |
           (SpeciesName == "QUNI" & CohortAge >= 35) |
           (SpeciesName == "QUPR2" & CohortAge >= 25) |
           (SpeciesName == "QUPH" & CohortAge >= 35) |
           (SpeciesName == "QURU" & CohortAge >= 25) |
           (SpeciesName == "QUST" & CohortAge >= 35) |
           (SpeciesName == "QUVE" & CohortAge >= 25)) %>%
  group_by(SpeciesName, CohortAge, Run, Time) %>%
  summarise(CohortBiomass = sum(CohortBiomass), .groups = "drop") %>%
  mutate(Weight = CohortAge / max(CohortAge, na.rm = TRUE)) %>%
  mutate(WeightedBiomass = CohortBiomass * Weight) %>%
  ungroup() %>%
  group_by(Time, Run) %>%
  summarise(Score = sum(WeightedBiomass), .groups = "drop")

fwrite(df, file = "/share/tcsi/lagoodal/R/Data/Fruit_Trees.csv")

# CULTURAL TREES
df <- combinedData %>%
  filter((SpeciesName == "QUAL" & CohortAge >= 30) |
           (SpeciesName == "CAAL27" & CohortAge >= 25) |
           (SpeciesName == "CAGL8" & CohortAge >= 30) |
           (SpeciesName == "LITU" & CohortAge >= 20) |
           (SpeciesName == "OXAR" & CohortAge >= 40)) %>%
  group_by(SpeciesName, CohortAge, Run, Time) %>%
  summarise(CohortBiomass = sum(CohortBiomass), .groups = "drop") %>%
  mutate(Weight = CohortAge / max(CohortAge, na.rm = TRUE)) %>%
  mutate(WeightedBiomass = CohortBiomass * Weight) %>%
  ungroup() %>%
  group_by(Time, Run) %>%
  summarise(Score = sum(WeightedBiomass), .groups = "drop")

rm(combinedData)

fwrite(df, file = "/share/tcsi/lagoodal/R/Data/Culture_Trees.csv")