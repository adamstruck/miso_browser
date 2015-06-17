## UI
library(shiny)
library(DT)

## Plotting
library(ggplot2)

## Data manipulation
## library(GenomicRanges)
library(RSQLite)
library(dplyr)
library(tidyr)
library(lazyeval)

## Gene info for lookups
## load("data/gene_lookup.RData")
## gene_list <- sort(unique(gene_lookup$symbol))

## Miso event lookup
## load("data/AltMISOevent_lookup.RData")
## load("data/AltMISOevent_lookup_granges.RData")

## sample metadata
sample_metadata <- read.table("data/sample_metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_metadata_columns <- colnames(sample_metadata)

## connect to sqlite db containing the MISO data
my_db <- src_sqlite("data/miso_summaries.sqlite", create = FALSE)
miso_sqlite <- tbl(my_db, "consolidatedMISOSummaries")
MISOdata_columns <- colnames(miso_sqlite)
