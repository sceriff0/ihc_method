library(tidyverse)
library(DESeq2)
library(here)
library(readxl)
library(fs)
library("data.table")

dds <- get(load(here("data", "counts.RData")))
dds <- DESeq(dds)

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
  "10338",               75,            NA,            NA
)

counts_data <-  counts(dds, normalized = TRUE) |>
  as_tibble(rownames = "GENE") |>
  pivot_longer(cols = -GENE, names_to = "Sample", values_to = "Expression") |>
  pivot_wider(names_from = GENE, values_from = Expression) |>
  filter(Sample %in% clinical_data$`ID CRF PRESERVE`)

# FlowPath now exports ONE csv per (patient, annotation): "<patient>_a<k>.csv".
# The cell rows are identical across a patient's files — only `Out_of_annotation`
# differs, precomputed for ANNOTATION_<k>. So membership is read straight from the
# flag (no sf point-in-polygon) and whole-slide metrics must use ONE copy of the
# cells, not all k files. `.cell_row` is a within-file row index used as a fallback
# join key when `cell_id` is absent.
process_patient <- function(csv_path) {

  stem  <- path_ext_remove(path_file(csv_path))          # "046_a1"
  parts <- str_match(stem, "^(.*)_a(\\d+)$")             # -> patient, k
  patient_id <- parts[, 2]
  ann_num    <- as.integer(parts[, 3])
  if (is.na(patient_id) || is.na(ann_num))
    stop(sprintf("unexpected flowpath filename (want <patient>_a<k>.csv): %s", stem))

  message(sprintf("LOADING PATIENT %s / ANNOTATION_%d", patient_id, ann_num))

  fread(csv_path) |>
    as_tibble() |>
    mutate(
      phenotype_clean = str_extract(phenotype, "(?<=\\().*?(?=\\))"),
      patient_id      = patient_id,
      annotation      = paste0("ANNOTATION_", ann_num),
      ann_num         = ann_num,
      .cell_row       = row_number()
    )
}

# Long table: one row per cell PER annotation file (carries Out_of_annotation).
# Feeds annotation-membership metrics.
load_ihc_annot <- function() {
  csv_paths    <- dir_ls("~/ihc_method/data/flowpath", glob = "*.csv")  # top level only; skips old/
  safe_process <- possibly(process_patient, otherwise = NULL, quiet = FALSE)
  map(csv_paths, safe_process) |> list_rbind()
}

# One deduped copy of each patient's cells (the lowest-index annotation file),
# stripped of the per-annotation columns. Feeds every whole-slide metric so cells
# are never counted k times.
dedupe_whole_slide <- function(annot) {
  annot |>
    group_by(patient_id) |>
    filter(ann_num == min(ann_num)) |>
    ungroup() |>
    select(-any_of(c("annotation", "ann_num", ".cell_row", "Out_of_annotation")))
}

ihc_annot <- load_ihc_annot()
ihc_data  <- dedupe_whole_slide(ihc_annot)

colnames(clinical_data)
colnames(neoplastic_data)
colnames(counts_data)
colnames(ihc_data)

