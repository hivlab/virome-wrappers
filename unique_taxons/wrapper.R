

stopifnot(length(snakemake@output) == 1)

#' Loading libraries
library(dplyr)
library(readr)
library(purrr)
library(tidyr)

parse_taxonomy <- function(input, output) {

  all_files <- tibble(file = input)
  queries <- all_files %>%
    mutate(data = map(file, read_csv, col_types = "cccdiiiiiiiddccccccccccccc")) %>%
    unnest()

  #' Hits per query
  #' Queries can have multple blast hits, we need to identify most plausble hits for each query.
  #' 1) By evalue,
  #' 2) Ranking by pident has second best impact.
  #' Also dropping gi column to select rows with unique query tax_id combinations.
  top_queries <- queries %>%
    group_by(query) %>%
    filter(dplyr::min_rank(evalue) == 1, dplyr::min_rank(desc(pident)) == 1)

  #' Selecting only query and assigned taxons
  unamb_top_queries <- top_queries %>%
    select(query, tax_id, parent_tax_id)

  #' Select rows with unique query tax_id combinations.
  unique_top_queries <- unamb_top_queries %>%
    ungroup %>%
    distinct()

  #' Filter by taxon abundance
  #' Calculate weighted taxon counts for group
  #+ tax-per-query
  abun_top_queries <- unique_top_queries %>%
    add_count(query, tax_id) %>%
    group_by(tax_id) %>%
    mutate(tax_per_sample = sum(n)) %>%
    group_by(query) %>%
    top_n(1, tax_per_sample)

  #' Summarise number of taxons by query.
  per_query <- abun_top_queries %>%
    select(query, tax_id) %>%
    distinct() %>%
    count(query)

  #' List of queries with multiple tax_ids but unique parent_taxi_id
  distinct_parent_taxid <- per_query %>%
    filter(n > 1) %>%
    left_join(unamb_top_queries) %>%
    select(query, parent_tax_id) %>%
    distinct() %>%
    count(query) %>%
    arrange(desc(n)) %>%
    filter(n == 1) %>%
    pull(query)

  #' Assigning parent_tax_id to tax_id for queries with multiple tax_ids but unique parent_taxi_id
  #' Dropping ambiguous queries
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
