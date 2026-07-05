library(tidyverse)
library(DESeq2)
library(here)
library(readxl)

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




