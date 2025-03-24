# Load libraries
library(tidyverse)

# STRUCTURAL DIVERSITY
# Read in data
file_path <- "/share/tcsi/lagoodal/R/Data/combinedData.csv"
combinedData <- data.table::fread(file_path)

# Calculate H for structural diversity
df <- combinedData %>%
  group_by(MapCode, Run, Time) %>%
  summarise(H = -sum((table(CohortAge) / length(CohortAge)) * log(table(CohortAge) / length(CohortAge))), .groups = "drop") %>%
  group_by(Run, Time) %>%
  summarise(H = mean(H), .groups = "drop")

data.table::fwrite(df, file = "/share/tcsi/lagoodal/R/Data/Structural_Diversity.csv")

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

data.table::fwrite(df, file = "/share/tcsi/lagoodal/R/Data/Fruit_Trees.csv")

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

data.table::fwrite(df, file = "/share/tcsi/lagoodal/R/Data/Culture_Trees.csv")
