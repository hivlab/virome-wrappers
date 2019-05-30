

stopifnot(length(snakemake@output) == 1)

message("Loading libraries.")
library(dplyr)
library(readr)
library(purrr)
library(tidyr)

parse_taxonomy <- function(input, output) {

  message("Importing Blast results\n")
  all_files <- tibble(file = unlist(input))
  queries <- all_files %>%
    mutate(data = map(file, read_csv, col_types = "cccdiiiiiiiddccccccccccccc")) %>%
    unnest()

  message("Hits per query\n
  Queries can have multple blast hits, we need to identify most plausble hits for each query.\n
  1) By evalue,\n
  2) Ranking by pident has second best impact.\n
  Also dropping gi column to select rows with unique query tax_id combinations.\n")
  top_queries <- queries %>%
    group_by(query) %>%
    filter(dplyr::min_rank(evalue) == 1, dplyr::min_rank(desc(pident)) == 1)
  
  message("Selecting rows with unique query tax_id combinations.\n")
  unique_top_queries <- top_queries %>%
    select(query, tax_id, parent_tax_id) %>% 
    ungroup %>%
    distinct()
  
  message("Filter by taxon abundance\n
  Calculating weighted taxon counts for group tax-per-query")
  abun_top_queries <- unique_top_queries %>%
    add_count(tax_id) %>%
    group_by(query) %>%
    top_n(1, n)
  
  message("Summarising number of taxons by query.\n")
  per_query <- abun_top_queries %>%
    select(query, tax_id) %>%
    distinct() %>%
    count(query)
  
  message("List of queries with multiple tax_ids but unique parent_taxi_id.\n")
  distinct_parent_taxid <- per_query %>%
    filter(n > 1) %>%
    left_join(abun_top_queries) %>%
    select(query, parent_tax_id) %>%
    distinct() %>%
    count(query) %>%
    filter(n == 1) %>%
    pull(query)
  
  message("Assigning parent_tax_id to tax_id for queries with multiple tax_ids but unique parent_taxi_id\n
  Dropping ambiguous queries.\n")
  unique_queries <- abun_top_queries %>%
    ungroup() %>%
    mutate(tax_id = case_when(
      query %in% distinct_parent_taxid ~ parent_tax_id,
      TRUE ~ tax_id
    )) %>%
    select(query, tax_id) %>%
    distinct() %>%
    add_count(query) %>%
    filter(n == 1) %>%
    select(query, tax_id)

  write_csv(unique_queries, output)
}

parse_taxonomy(input = snakemake@input, output = snakemake@output)
