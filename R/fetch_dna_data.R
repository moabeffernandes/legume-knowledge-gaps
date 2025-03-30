# Helper functions (Fetch DNA data)
# Moabe Fernandes

# Function to fetch data for a given DNA region using entrez_search & entrez_summary
fetch_dna_data <- function(region, search_term, batch_size = 500) {
  message(glue::glue("Fetching data for region: {region}"))
  
  # Initial search to get record count and web history
  initial_search <- rentrez::entrez_search(db = "nuccore",
                                           term = search_term,
                                           use_history = TRUE)
  total_count <- as.numeric(initial_search$count)
  web_history <- initial_search$web_history
  
  # Prepare to collect summaries in batches
  batches <- ceiling(total_count / batch_size)
  all_summaries <- vector("list", batches)
  
  # Iterate over batches
  for (j in seq(1, total_count, by = batch_size)) {
    batch_summaries <- rentrez::entrez_summary(db = "nuccore", 
                                               web_history = web_history,
                                               retmax = batch_size, 
                                               retstart = j - 1)
    all_summaries[[ceiling(j / batch_size)]] <- batch_summaries
  }
  
  # Combine summaries and extract species information
  all_summaries <- do.call(c, all_summaries)
  species_info <- sapply(all_summaries, function(x) c(x$organism, x$title))
  
  species_df <- as.data.frame(t(species_info), stringsAsFactors = FALSE)
  colnames(species_df) <- c("species", "title")
  species_df <- tibble::as_tibble(species_df)
  species_df <- species_df %>% dplyr::mutate(regions = region)
  
  return(species_df)
}

