# Helper functions
# Moabe Fernandes

# Function to extract publication years
clean_year <- function(year_string) {
  year_string %>%
    str_remove_all("\\(|\\)") %>%  # Remove parentheses
    str_extract("[[:digit:]]+$") %>%  # Extract the last numeric part
    as.numeric()  # Convert to numeric
}


# Function to split authors into individual entries
split_authors <- function(data) {
  tibble(
    basionym_published = rep(data$basionym_published, lengths(strsplit(data$basionym_authors, ", | & "))),
    interval = rep(data$interval, lengths(strsplit(data$basionym_authors, ", | & "))),
    author_names = unlist(strsplit(data$basionym_authors, ", | & "))
  )
}

