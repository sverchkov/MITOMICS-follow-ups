# Main target script for running everything.
#
# Author: Yuriy Sverchkov

src_folder <- here::here()
data_folder <- here::here("data")
results_folder <- here::here("results")

metadata_folder <- file.path(data_folder, "metadata")

library(targets)

source(file.path(src_folder, "functions.R"))

options(tidyverse.quiet = TRUE)
tar_option_set(packages= c("magrittr", "tidyverse", "ggplot2", "qvalue", "processx"))

# Targets list
list(

  # Templates
  tar_target(
    embedding_html_template,
    file.path(src_folder, "templates", "embedding_plot.html"),
    format = "file"
  ),

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
  ),

  # Write metabolites files
  tar_target(
    metabolites_data_file,
    {
      fp <- file.path(data_folder, 'processed', 'metabolites_data.csv')
      metabolites$data %>%
        mutate(`Q Value` = qvalue(`P Value`)$qvalue) %>%
        write_csv(fp)
      fp
    },
    format = 'file'
  ),

  tar_target(
    metabolites_metadata_file,
    {
      fp <- file.path(data_folder, 'processed', 'metabolites_metadata.csv')
      metabolites$metadata %>%
        write_csv(fp)
      fp
    },
    format = 'file'
  ),

  # Make a matrix of q-adjusted relative differences    
  tar_target(
    qardm,
    derive_QARD_matrix(molecules_df)
  ),

  # Save QARD matrix to file
  tar_target(
    qardm_file,
    {
      fp <- file.path(data_folder, 'processed', 'qardm.csv')
      qardm %>% as_tibble(rownames='ID') %>% write_csv(fp)
      fp
    },
    format = 'file'
  ),

  # Run PaCMAP
  tar_target(
    pacmap_euclidean_embedding_file,
    {
      fp <- file.path(data_folder, 'processed', 'pacmap_euclidean_embedding.csv')
      run('python', c('-m', 'pyscripts.pacmap', qardm_file, fp, 'euclidean'), echo = TRUE)
      fp
    },
    format = 'file'
  ),

  tar_target(
    pacmap_angular_embedding_file,
    {
      fp <- file.path(data_folder, 'processed', 'pacmap_angular_embedding.csv')
      run('python', c('-m', 'pyscripts.pacmap', qardm_file, fp, 'angular'), echo = TRUE)
      fp
    },
    format = 'file'
  ),

  # Read PACMAP embedding from file
  tar_target(pacmap_angular_embedding, read_csv(pacmap_angular_embedding_file)),
  tar_target(pacmap_euclidean_embedding, read_csv(pacmap_euclidean_embedding_file)),

  # Create visualization
  tar_target(
    pacmap_angular_plot,
    plot_embedding2(
      x = pacmap_angular_embedding$`Component 1`,
      y = pacmap_angular_embedding$`Component 2`,
      id = pacmap_angular_embedding$ID,
      metadata = metadata_for_plot(lipids$metadata, metabolites$metadata, proteins$metadata),
      colors = c(
          "Lipid" = "#2cb781",
          "Metabolite" = "#e3b74a",
          "Protein" = "#99b0b5",
          "Protein (Mitochondria)" = "#003947",
          "Protein (Mitoribosome)" = "#9d2063",
          "Protein (OxPhos)" = "#ed2c51"),
      html_template = embedding_html_template,
      html_out = file.path(results_folder, "angular_embedding.html")
    ),
    format = "file"
  ),

  tar_target(
    pacmap_euclidean_plot,
    plot_embedding2(
      x = pacmap_euclidean_embedding$`Component 1`,
      y = pacmap_euclidean_embedding$`Component 2`,
      id = pacmap_euclidean_embedding$ID,
      metadata = metadata_for_plot(lipids$metadata, metabolites$metadata, proteins$metadata),
      colors = c(
          "Lipid" = "#2cb781",
          "Metabolite" = "#e3b74a",
          "Protein" = "#99b0b5",
          "Protein (Mitochondria)" = "#003947",
          "Protein (Mitoribosome)" = "#9d2063",
          "Protein (OxPhos)" = "#ed2c51"),
      html_template = embedding_html_template,
      html_out = file.path(results_folder, "eucledian_embedding.html")
    ),
    format = "file"
  )

)

