library(tidyverse)
library(DESeq2)
library(here)
library(readxl)
library(fs)
library("data.table")

dds <- get(load(here("data", "counts.RData")))
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)

clinical_data <- read_excel(here("data", "clinical_data.xlsx")) |>
  filter(!is.na(`ID PATIENT`)) |>
  mutate(`ID CRF PRESERVE` = gsub("-", ".", `ID CRF PRESERVE`))

neoplastic_data <- tribble(
  ~SAMPLE, ~ANNOTATION_1, ~ANNOTATION_2, ~ANNOTATION_3,
  "24086",              80,            70,            70,
  "15897",              70,            NA,            NA,
  "052",                70,            50,            NA,
  "046",                50,            60,            50,
  "5456",               70,            80,            70,
  "103382",              75,            NA,            NA
)

counts_data <-  counts(dds, normalized = TRUE) |>
  as_tibble(rownames = "GENE") |>
  pivot_longer(cols = -GENE, names_to = "Sample", values_to = "Expression") |>
  pivot_wider(names_from = GENE, values_from = Expression) |>
  filter(Sample %in% clinical_data$`ID CRF PRESERVE`)

process_patient <- function(csv_path) {

  patient_id <- path_ext_remove(path_file(csv_path))   # patient1.csv -> patient1

  message(sprintf("LOADING PATIENT: %s", patient_id))

  ihc_data <- fread(csv_path) |> as_tibble()

  ihc_clean_data <- ihc_data |>
    mutate(
      phenotype_clean = str_extract(phenotype, "(?<=\\().*?(?=\\))"),
      patient_id      = patient_id
    )

  ihc_clean_data
}

load_ihc_data <- function() {

  csv_paths  <- dir_ls("~/ihc_method/data/flowpath", glob = "*.csv")   # point this at folder X
  safe_process <- possibly(process_patient, otherwise = NULL, quiet = FALSE)

  map(csv_paths, safe_process) |> list_rbind()
}

ihc_data <- load_ihc_data()

colnames(clinical_data)
colnames(neoplastic_data)
colnames(counts_data)
colnames(ihc_data)

