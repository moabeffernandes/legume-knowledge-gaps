# Estimating Legume Total Species Richness and Description Rates Growth
# Moabe Fernandes

# Description:
# This script estimates the completeness of the legume species inventory using 
# various non-linear models. It also computes species description rates across 
# geographical regions based on data from the rWCVP package.

# Libraries ----
library(tidyverse)  # Data wrangling
library(rWCVP)  # World Checklist of Vascular Plants (WCVP)
library(nlme) # Generalised nonlinear least squares regression

# load helper functions
source("R/helper_functions.R")

# Section 1: Create Dataset ----

# Generate a list of all accepted legume species
legume_wcvp <- rWCVPdata::wcvp_names %>%
  filter(
    family == "Fabaceae", 
    taxon_rank == "Species", 
    taxon_status == "Accepted",
    !str_detect(taxon_name, "×")  # Exclude hybrids
  ) %>%
  select(
    accepted_name_id = plant_name_id,
    accepted_name = taxon_name,
    accepted_authors = taxon_authors,
    accepted_name_published = first_published,
    basionym_name_id = basionym_plant_name_id,
    genus
  ) %>%
  distinct() %>%
  arrange(accepted_name)

# Load and select necessary columns from the WCVP dataset
complete_wcvp <- rWCVPdata::wcvp_names %>%
  select(
    plant_name_id,
    basionym_name = taxon_name,
    basionym_authors = taxon_authors,
    basionym_published = first_published
  )

# Join basionym information to the legume dataset
legume_wcvp <- legume_wcvp %>%
  left_join(complete_wcvp, by = c("basionym_name_id" = "plant_name_id")) %>%
  mutate(
    is_basionym = is.na(basionym_name),  # Flag whether the record is a basionym
    basionym_published = coalesce(
      basionym_published, 
      accepted_name_published
    ),
    basionym_name = coalesce(
      basionym_name, 
      accepted_name
    ),
    basionym_authors = coalesce(
      basionym_authors, 
      accepted_authors
    )
  )


# Clean and extract publication years
legume_wcvp <- legume_wcvp %>%
  mutate(
    basionym_published = clean_year(basionym_published),
    accepted_name_published = clean_year(accepted_name_published)
  )

# Exclude species with missing publication dates
legume_wcvp <- legume_wcvp %>%
  filter(!is.na(basionym_published))

# Retain only important columns
legume_wcvp <- legume_wcvp %>%
  select(
    accepted_name_id,
    accepted_name,
    accepted_authors,
    accepted_name_published,
    basionym_name,
    basionym_authors,
    basionym_published,
    is_basionym
  )



# Section 2: Estimate Total Number of Species ----

# Prepare Dataset ----

# Define interval lengths (every 1, 5 and 10 years)
interval_lengths <- c(1, 5, 10)

# Initialise an empty dataframe to store results
final_results_model <- tibble(
  model = character(),
  time_range = numeric(),
  lower = numeric(),
  total_species = numeric(),
  upper = numeric(),
  AIC = numeric(),
  interval_length = numeric()
)

# Iterate over interval lengths
for (interval_length in interval_lengths) {
  
  # Define the range for the intervals
  start_year <- 1755
  end_year <- 2020
  interval_breaks <- seq(start_year, end_year, by = interval_length)
  
  # Create intervals and assign them to a new column
  legume_wcvp <- legume_wcvp %>%
    mutate(interval = cut(basionym_published, breaks = interval_breaks,
                          right = TRUE, include.lowest = TRUE))
  
  # Calculate the number of species described in each interval and cumulative species count
  species_by_interval <- legume_wcvp %>%
    group_by(interval) %>%
    summarise(species_described = n(), .groups = 'drop') %>%
    arrange(interval) %>%
    mutate(cumulative_species = cumsum(species_described)) %>%
    mutate(time = as.numeric(sub("^[\\[(](\\d+),.*", "\\1", interval))) %>%
    filter(!is.na(time))
  
  # Calculate the number of authors working in each time interval
  author_by_interval <- legume_wcvp %>% 
    select(basionym_published, interval, basionym_authors) %>%
    mutate(basionym_authors = sub(".* ex ", "", basionym_authors)) %>%
    split_authors() %>%
    mutate(author_names = str_squish(author_names)) %>%
    select(interval, author_names) %>%
    distinct() %>%
    filter(!is.na(interval)) %>%
    group_by(interval) %>%
    summarise(authors = n(), .groups = 'drop') %>%
    arrange(interval)
  
  # Join species and author information
  species_by_interval <- species_by_interval %>%
    left_join(author_by_interval, by = "interval") %>%
    select(time, authors, species_described, cumulative_species)
  
  # Define time ranges for model fitting
  time_range <- seq(from = 1975, to = 2015, by = 10)
  initial_number <- 30000  # Adjust based on the dataset's characteristics
  
  # Fit models
  for (i in 1:length(time_range)) {
    
    # Logistic Model
    logistic_model <- gnls(species_described ~ (a + b * cumulative_species) * (total_species - cumulative_species),
                           data = subset(species_by_interval, time <= time_range[i]),
                           start = list(total_species = initial_number, a = 1e-07, b = 1e-08),
                           weights = varPower(), verbose = TRUE,
                           control = gnlsControl(returnObject = TRUE, minScale = 1e-05,
                                                 tolerance = 1e-04, nlsMaxIter = 3))
    logistic_estimates <- intervals(logistic_model)$coef[1,]
    logistic_AIC <- AIC(logistic_model)
    
    # Joppa Model
    joppa_model <- gnls(species_described ~ authors * (a + b * time) * (total_species - cumulative_species),
                        data = subset(species_by_interval, time <= time_range[i]),
                        start = list(total_species = initial_number, a = -1e-07, b = 1e-07),
                        weights = varPower(), verbose = TRUE,
                        control = gnlsControl(returnObject = TRUE, minScale = 1e-05,
                                              tolerance = 0.001, nlsMaxIter = 3))
    joppa_estimates <- intervals(joppa_model)$coef[1,]
    joppa_AIC <- AIC(joppa_model)
    
    # Lu & He Model
    lu_model <- gnls(species_described ~ authors * (a + b * species_described) * (total_species - cumulative_species), 
                     data = subset(species_by_interval, time <= time_range[i]),
                     start = list(total_species = initial_number, a = 1e-07, b = 1e-07), 
                     weights = varPower(), verbose = TRUE, 
                     control = gnlsControl(returnObject = TRUE, minScale = 1e-05,
                                           tolerance = 0.001, nlsMaxIter = 3))
    lu_estimates <- intervals(lu_model)$coef[1,]
    lu_AIC <- AIC(lu_model)
    
    # Combine results
    results <- tibble(
      model = c("Logistic", "Joppa et al.", "Lu & He"),
      time_range = rep(time_range[i], 3),
      lower = c(logistic_estimates[1], joppa_estimates[1], lu_estimates[1]),
      total_species = c(logistic_estimates[2], joppa_estimates[2], lu_estimates[2]),
      upper = c(logistic_estimates[3], joppa_estimates[3], lu_estimates[3]),
      AIC = c(logistic_AIC, joppa_AIC, lu_AIC),
      interval_length = interval_length
    )
    
    final_results_model <- bind_rows(final_results_model, results)
  }
  
}

# Round numbers and save
final_results_model <- final_results_model %>%
  mutate(
    lower = round(lower, 0),
    total_species = round(total_species, 0),
    upper = round(upper, 0),
    AIC = round(AIC, 2)
  )

# Save model results
write_csv(final_results_model, "results/species_discovery_models.csv")



# Section 3: Compute Regional Growth Rates ----

# Generate geographic dataset
legume_distribution <- wcvp_checklist(taxon = "Fabaceae", taxon_rank = "family") %>%
  filter(taxon_status == "Accepted") %>%
  select(accepted_plant_name_id, area_code_l3) %>%
  left_join(legume_wcvp, by = c("accepted_plant_name_id" = "accepted_name_id")) %>%
  filter(!is.na(basionym_published))

# Annual growth rates
spp_cumulative <- expand_grid(
  area_code_l3 = unique(legume_distribution$area_code_l3),
  year = 1753:2024
) %>%
  left_join(
    legume_distribution %>%
      group_by(area_code_l3, year = basionym_published) %>%
      summarise(spp_number = n(), .groups = "drop"),
    by = c("area_code_l3", "year")
  ) %>%
  mutate(spp_number = replace_na(spp_number, 0)) %>%
  group_by(area_code_l3) %>%
  arrange(year) %>%
  mutate(
    cumulative_species = cumsum(spp_number),
    growth_rate = (cumulative_species - lag(cumulative_species, default = 0)) / 
      lag(cumulative_species, default = 1) * 100
  )

# Compute Average Growth Rates for Selected Intervals
average_growth_rate <- tibble()

for (interval in c(1970, 1980, 1990, 2000, 2010)) {
  message(glue::glue("Calculating growth rates for interval starting: {interval}"))
  
  interval_growth <- spp_cumulative %>%
    filter(year >= interval, year <= 2020) %>%
    group_by(area_code_l3) %>%
    summarise(avg_growth = mean(growth_rate, na.rm = TRUE)) %>%
    mutate(period = interval)
  
  average_growth_rate <- bind_rows(average_growth_rate, interval_growth)
}

# Transform to Wide Format
description_rates <- average_growth_rate %>%
  pivot_wider(names_from = period, values_from = avg_growth) %>%
  rename_with(~ paste0("si_", .), -area_code_l3)

# Merge with WGSRPD Data
description_rates <- rWCVPdata::wgsrpd3 %>%
  sf::st_drop_geometry() %>%
  select(wgsrpd_level3_name = LEVEL3_NAM, wgsrpd_level3_code = LEVEL3_COD) %>%
  left_join(description_rates, by = c("wgsrpd_level3_code" = "area_code_l3")) %>%
  mutate(across(starts_with("si_"), ~ replace_na(round(., 3), 0)))

# Save Results
write_csv(description_rates, "results/new_species_description_rates.csv")
