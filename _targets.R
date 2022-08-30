# Main target script for running everything.
#
# Author: Yuriy Sverchkov

src_folder <- here::here()
data_folder <- here::here("data")

metadata_folder <- file.path(data_folder, "metadata")

library(targets)

source(file.path(src_folder, "functions.R"))

options(tidyverse.quiet = TRUE)
tar_option_set(packages= c("magrittr", "tidyverse", "ggplot2", "patchwork", "qvalue"))

# Targets list
list(

  # Metadata files

  tar_target(
    mxp_genes_file,
    file.path(metadata_folder, "mxp_genes.csv"),
    format = "file"
  ),

  tar_target(
    metadata_xlsx,
    file.path(metadata_folder, 'H3K Project Master Lists.xlsx'),
    format = "file"
  ),

  # Data files

  tar_target(
    data_xlsx,
    file.path(data_folder, 'Supplementary Table 3.xlsx'),
    format = "file"
  ),

  # QC

  tar_target(
    excluded_lines,
    c(
      'FUNDC2-KO2', 'HDHD3-KO1', 'HDHD3-KO2', 'MIPEP-KO', 'MPC1-KO2', 'MTRES1-KO', 'PISD-KO1', 'PISD-KO2',
      'SFXN3-KO1', 'SFXN3-KO2', 'SPATA20-KO1', 'SPATA20-KO2', 'TOMM7-KO2'
    )
  ),

  # Load data
  tar_target(lipids, read_lipids(data_xlsx, "Lipidomics log2FCsPvaluesSDs")),
  tar_target(metabolites, read_metabolites(data_xlsx, "Metabolomics log2FCsPvaluesSDs")),
  tar_target(proteins, read_proteins(data_xlsx, "Proteomics log2FCsPvaluesSDs")),

  # Make a combined molecule type dataframe
  tar_target(
    molecule_types,
    bind_rows(
      lipids$metadata %>% select(ID) %>% mutate(`molecule type` = 'Lipid'),
      metabolites$metadata %>% select(ID) %>% mutate(`molecule type` = 'Metabolite'),
      proteins$metadata %>% select(ID) %>% mutate(`molecule type` = 'Protein')
    )
  ),

  # Make combined molecules frame, and compute q-values
  tar_target(
    molecules_df,
    bind_rows(lipids$data, metabolites$data, proteins$data) %>%
      filter(
        !str_detect(`cell line`, '^P?WT'),
        !(`cell line` %in% excluded_lines)) %>%
      mutate(`Q Value` = qvalue(`P Value`)$qvalues)
  )
)

