# ============================================================================
# per_annotation_flag.R — flag-based per-annotation IHC metrics
# ============================================================================
# Builds the two driver frames that clinical_data_per_annotation.Rmd consumes
# (ihc_tumor_per_ann, ihc_tumor_union) from the per-annotation FlowPath flag files
#   data/flowpath/per_annotation/<patient>_a<k>.csv
# — one file per (patient, annotation), each carrying the cell centroids + the
# precomputed Out_of_annotation flag (and the usual FlowPath phenotype /
# <marker>_sign columns). Membership here is the FLAG (inside = Out_of_annotation
# is falsey), NOT sf point-in-polygon.
#
# Depends on validation_helpers.R (slide_key, .is_outside, region_ratios_area) —
# source THAT first. The returned frames match the schema of
# ihc_annotation_metrics() so the whole downstream analysis works unchanged; with
# no polygon the area_mm2 / dens_*_per_mm2 columns come back NA
# (region_ratios_area(..., 0, ...)).
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
  library(readr)
  library(stringr)
  library(purrr)
})

PER_ANNOT_DIR <- here::here("data", "flowpath", "per_annotation")

# Read every <patient>_a<k>.csv (top-level glob, so old/ is skipped) into one long
# cells table keyed by (patient_id, annotation). patient_id = slide_key of the stem
# before "_a" (shared numeric slide space); annotation = "ANNOTATION_<k>" to line up
# with neoplastic_data's ANNOTATION_1..3 after pivot. phenotype_clean is the
# parenthetical label in `phenotype` (same extraction as test.R's process_patient).
load_per_annotation_cells <- function(dir = PER_ANNOT_DIR) {
  files <- fs::dir_ls(dir, glob = "*.csv")            # non-recursive: old/ excluded
  purrr::map_dfr(as.character(files), function(path) {
    stem <- fs::path_ext_remove(fs::path_file(path))  # "046_a1"
    m    <- stringr::str_match(stem, "^(.*)_a(\\d+)$") # [ , patient, k ]
    if (is.na(m[1, 1])) {
      warning("per_annotation: skipping ", path, " — name is not <patient>_a<k>.csv")
      return(tibble::tibble())
    }
    cells <- readr::read_csv(path, show_col_types = FALSE)
    if (!all(c("Out_of_annotation", "phenotype") %in% names(cells))) {
      warning("per_annotation: skipping ", path, " — missing Out_of_annotation / phenotype")
      return(tibble::tibble())
    }
    cells |>
      dplyr::mutate(
        patient_id      = slide_key(m[1, 2]),
        annotation      = paste0("ANNOTATION_", as.integer(m[1, 3])),
        phenotype_clean = stringr::str_extract(phenotype, "(?<=\\().*?(?=\\))")
      )
  })
}

# inside = the FlowPath flag says the cell is NOT out of the annotation.
.flag_inside <- function(cells) !.is_outside(cells$Out_of_annotation)

# Cell identity for union-dedup: cell_id if the export carries it, else the centroid.
.cell_key_cols <- function(cells) {
  if ("cell_id" %in% names(cells)) "cell_id"
  else intersect(c("centroid_x", "centroid_y"), names(cells))
}

# PER-ANNOTATION metrics: one row per (patient_id, annotation). Inside cells only,
# no polygon area (area_mm2 / dens_* = NA); source = "flag". Same column schema as
# ihc_annotation_metrics(..., scope = "per_annotation").
ihc_flag_metrics_per_annotation <- function(ann_cells, um_per_px = 0.325) {
  ann_cells |>
    dplyr::group_by(patient_id, annotation) |>
    dplyr::group_split() |>
    purrr::map_dfr(function(cells) {
      region_ratios_area(cells[.flag_inside(cells), , drop = FALSE], 0, um_per_px) |>
        dplyr::mutate(patient_id = cells$patient_id[1],
                      annotation = cells$annotation[1],
                      source     = "flag", .before = 1)
    })
}

# UNION metrics: one row per patient. Pool inside cells across the patient's
# annotation files and dedup to one row per cell (so a cell inside two annotations
# counts once), then the same area-free ratios. Same schema as scope = "union".
ihc_flag_metrics_union <- function(ann_cells, um_per_px = 0.325) {
  key_cols <- .cell_key_cols(ann_cells)
  ann_cells |>
    dplyr::group_by(patient_id) |>
    dplyr::group_split() |>
    purrr::map_dfr(function(cells) {
      inside <- cells[.flag_inside(cells), , drop = FALSE]
      if (length(key_cols))
        inside <- dplyr::distinct(inside, dplyr::across(dplyr::all_of(key_cols)),
                                  .keep_all = TRUE)
      region_ratios_area(inside, 0, um_per_px) |>
        dplyr::mutate(patient_id = cells$patient_id[1],
                      annotation = "union", source = "flag", .before = 1)
    })
}

# Per-patient UNION of inside cells (deduped) as a plain cells table — the
# in-annotation analogue of whole-slide ihc_data for ihc_lineage_fraction() /
# ihc_marker_fraction() in the biological-validity section.
union_inside_cells <- function(ann_cells) {
  key_cols <- .cell_key_cols(ann_cells)
  inside <- ann_cells[.flag_inside(ann_cells), , drop = FALSE]
  if (length(key_cols))
    inside <- inside |>
      dplyr::group_by(patient_id) |>
      dplyr::distinct(dplyr::across(dplyr::all_of(key_cols)), .keep_all = TRUE) |>
      dplyr::ungroup()
  inside
}
