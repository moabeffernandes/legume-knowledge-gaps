# Analysing Legume Evolutionary Coverage Across Lineages and Geographic Space
# Moabe Fernandes

# Description: 
# This script retrieves DNA sequence data for several genetic markers from NCBI,
# cleans and resolves species names against a taxonomic backbone, and then
# prepares datasets for subsequent phylogenetic analyses.

# Libraries ----
library(rentrez)
library(dplyr)
library(stringr)
library(WorldFlora)
library(rWCVP)
library(tidyr)
library(readr)
library(purrr)

# Load helper functions
source("R/fetch_dna_data.R")
source("R/resolve_wcvp_multi_matches.R")


# Section 1:Load and Process DNA Sequence Data ----

# Load DNA Sequence Data ----

# Define DNA regions and corresponding search strings
dna_regions <- c("matK", "psbA", "rbcL", "trnL", "ITS")
search_terms <- paste0(dna_regions, "[gene]")
search_terms[5] <- "internal transcribed spacer"
search_terms <- paste0(search_terms, " AND Fabaceae[orgn]")

# Loop over each region using the helper function; accumulate results
complete_regions <- purrr::map2_df(dna_regions, search_terms, fetch_dna_data)

# Remove duplicate entries
complete_regions <- complete_regions %>%
  select(species, regions) %>%
  distinct()


# Clean Search Results ----

# Create a new column for rank based on patterns in species names
complete_regions <- complete_regions %>%
  mutate(rank = case_when(
    str_detect(species, " cf\\.") ~ "confer",
    str_detect(species, " aff\\.") ~ "affinis",
    str_detect(species, " sp\\.") ~ "unidentified",
    str_detect(species, " x ") ~ "hybrid",
    str_detect(species, " hybrid ") ~ "hybrid",
    TRUE ~ "complete"
  )) %>%
  # Keep only complete names that can be connected to species
  filter(rank == "complete") %>%
  select(species, regions) %>%
  # Exclude species with round brackets (none accepted) and clean names
  filter(!str_detect(species, "\\("))

# Clean species name to avoid issues in name matching
complete_regions <- complete_regions %>% 
  mutate(species = species %>% 
           str_remove_all("\\[|\\]") %>%
           str_remove_all("[:digit:]+") %>%
           str_squish())


# Name Matching and Resolution ----

# Prepare names for matching
unique_names <- tibble(species = unique(complete_regions$species))
unique_names <- WFO.prepare(unique_names$species)

# Keep only original and edited names
unique_names <- unique_names %>% 
  select(spec.full, spec.name)

# Load the legume taxonomic backbone
fabaceae_taxonomy <- rWCVPdata::wcvp_names %>% 
  filter(family == "Fabaceae")

# Match names using rWCVP package
matches <- wcvp_match_names(unique_names,
                            wcvp_names = fabaceae_taxonomy,
                            name_col = "spec.name")

# Check for multiple matches
count(matches, multiple_matches)

# Resolve multiple matches by sourcing an external helper function
auto_resolved <- matches %>%
  nest_by(spec.full) %>%
  mutate(data = list(resolve_wcvp_multi_matches(data))) %>%
  unnest(col = data) %>%
  ungroup()

# Check match summary
count(auto_resolved, resolved_match_type)

# Exclude unresolved names and join accepted taxon information
auto_resolved <- auto_resolved %>%
  filter(resolved_match_type != "could_not_resolve") %>%
  select(spec.full, spec.name, wcvp_accepted_id, resolved_match_type) %>%
  left_join(select(rWCVPdata::wcvp_names, plant_name_id, taxon_rank, taxon_status, parent_plant_name_id),
            by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  mutate(accepted_species_id = case_when(
    taxon_rank != "Species" ~ parent_plant_name_id,
    TRUE ~ wcvp_accepted_id
  )) %>%
  select(spec.full, spec.name, accepted_species_id) %>%
  left_join(select(rWCVPdata::wcvp_names, plant_name_id, taxon_rank, taxon_status,
                   taxon_name, taxon_authors, parent_plant_name_id),
            by = c("accepted_species_id" = "plant_name_id")) %>%
  filter(taxon_status == "Accepted")

# Merge resolved names into the complete_regions data
complete_regions <- complete_regions %>%
  left_join(auto_resolved, by = c("species" = "spec.full")) %>%
  select(taxon_name, taxon_authors, accepted_species_id, regions) %>%
  distinct() %>%
  filter(!is.na(taxon_name))


# Prepare Data to Create Summaries ----

# Create checklist of accepted legume species
accepted_legumes <- rWCVPdata::wcvp_names %>% 
  filter(family == "Fabaceae",
         taxon_rank == "Species",
         taxon_status == "Accepted",
         !str_detect(taxon_name, "×")) %>% 
  select(plant_name_id, genus, taxon_name, taxon_authors)

# Mark species presence for each DNA region with progress messages
accepted_legumes <- accepted_legumes %>% 
  mutate(
    ITS    = plant_name_id %in% pull(filter(complete_regions, regions == "ITS"), accepted_species_id),
    matK   = plant_name_id %in% pull(filter(complete_regions, regions == "matK"), accepted_species_id),
    rbcL   = plant_name_id %in% pull(filter(complete_regions, regions == "rbcL"), accepted_species_id),
    trnL   = plant_name_id %in% pull(filter(complete_regions, regions == "trnL"), accepted_species_id),
    psbA   = plant_name_id %in% pull(filter(complete_regions, regions == "psbA"), accepted_species_id),
    overall = plant_name_id %in% unique(complete_regions$accepted_species_id)
  )

# Add subfamily information
accepted_legumes <- accepted_legumes %>% 
  left_join(read_csv("data/legume_genera_subfamily.csv"))



# Section 2: Summarise Marker Coverage per Genus and Subfamily ----

# Summaries at the genus level ----
dna_summary <- accepted_legumes %>% 
  group_by(genus) %>% 
  summarise(
    total_spp = n(),
    across(c(ITS, matK, rbcL, trnL, psbA, overall),
           ~ 1 - sum(.x, na.rm = TRUE) / total_spp,
           .names = "ec_{.col}"),
    .groups = "drop"
  )

# Save results per genus
write_csv(dna_summary, "results/evolutionary_coverage_genus.csv")

# Summaries at the subfamily level ----
dna_summary_subfamily <- accepted_legumes %>% 
  group_by(subfamily) %>% 
  summarise(
    total_spp = n(),
    across(c(ITS, matK, rbcL, trnL, psbA, overall),
           ~ 1 - sum(.x, na.rm = TRUE) / total_spp,
           .names = "ec_{.col}"),
    .groups = "drop"
  )


# Save results
write_csv(dna_summary_subfamily, "results/evolutionary_coverage_subfamily.csv")



# Section 3: Summarise Marker Coverage per Botanical Country ----

# Load Legume Occurrence Dataset using wcvp_checklist()
legume_distribution <- wcvp_checklist(
  taxon = "Fabaceae", 
  taxon_rank = "family",
  introduced = FALSE, 
  extinct = FALSE,
  location_doubtful = FALSE, 
  hybrids = FALSE,
  infraspecies = FALSE
) %>% 
  filter(taxon_status == "Accepted") %>%  # Exclude unplaced taxa
  select(accepted_plant_name_id, taxon_name, taxon_authors, area_code_l3)

# Create a geographic summary of total species per botanical country
geographic_summary <- legume_distribution %>% 
  group_by(area_code_l3) %>% 
  summarise(total_spp = n(), .groups = "drop")

# Identify species with available DNA markers from the accepted_legumes dataset
spp_dna <- accepted_legumes %>% 
  filter(ITS | matK | rbcL | trnL | psbA) %>% 
  pull(plant_name_id)

# Summarise species with DNA data per botanical country and calculate incompleteness
geographic_summary <- geographic_summary %>% 
  left_join(
    legume_distribution %>% 
      filter(accepted_plant_name_id %in% spp_dna) %>% 
      group_by(area_code_l3) %>% 
      summarise(spp_dna = n(), .groups = "drop"),
    by = "area_code_l3"
  ) %>% 
  mutate(
    spp_dna = if_else(is.na(spp_dna), 0L, spp_dna),  # Replace NA with 0 for countries with no DNA data
    ec_overall = 1 - (spp_dna / total_spp)
  ) %>% 
  select(area_code_l3, ec_overall)

# Save the geographic summary
write_csv(geographic_summary, "results/geographic_summary.csv")
