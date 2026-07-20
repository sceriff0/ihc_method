# ============================================================================
# per_annotation_flag.R — flag-based per-annotation IHC metrics
# ============================================================================
# Builds the driver frames that the clinical_data per-annotation Rmds consume
# (ann_cells, ihc_tumor_per_ann, ihc_tumor_union) from the FlowPath exports under
#   data/flowpath/per_annotation/<patient>_a<k>.csv        (one file per annotation)
#   data/flowpath/per_annotation/old/<patient>_a<k>.csv    (superseding re-exports)
#   data/flowpath/<patient>.csv                            (whole-slide fallback)
# Membership is ALWAYS the precomputed FlowPath flag (inside = Out_of_annotation is
# falsey), NEVER sf point-in-polygon — so no geojson is read here and every
# area_mm2 / dens_*_per_mm2 column comes back NA (region_ratios_area(..., 0, ...)).
#
# THREE CELL SOURCES, one long table
#   per-annotation  annotation = "ANNOTATION_<k>", source = "flag"
#   old/ overlay    same, but REPLACES the top-level file for the same
#                   (patient, annotation) when an overlay dir is passed
#   whole-slide     patients with NO per-annotation csv at all; the slide's single
#                   Out_of_annotation flag is one region: annotation = "csv",
#                   source = "csv" — mirroring ihc_annotation_metrics()'s fallback
#                   convention in validation_helpers.R so the two analyses line up.
# Because the fallback rows live in the SAME cells table, the metric functions below
# pick them up with no branching. Note that annotation "csv" never joins to the
# pathologist's ANNOTATION_k, so fallback patients enter the UNION concordance but
# not the per-annotation one — same behaviour as clinical_data.Rmd.
#
# Depends on validation_helpers.R (slide_key, .is_outside, region_ratios_area) —
# source THAT first. The returned frames match the schema of ihc_annotation_metrics()
# so the whole downstream analysis works unchanged.
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
  library(readr)
  library(stringr)
  library(purrr)
})

PER_ANNOT_DIR     <- here::here("data", "flowpath", "per_annotation")
PER_ANNOT_OLD_DIR <- here::here("data", "flowpath", "per_annotation", "old")
FLOWPATH_DIR      <- here::here("data", "flowpath")

# Read every <patient>_a<k>.csv directly in `dir` (non-recursive, so old/ and new/
# are only read when passed explicitly) into one long cells table keyed by
# (patient_id, annotation). patient_id = slide_key of the stem before "_a" (shared
# numeric slide space); annotation = "ANNOTATION_<k>" to line up with
# neoplastic_data's ANNOTATION_1..3 after pivot. phenotype_clean is the parenthetical
# label in `phenotype` (same extraction as test.R's process_patient). `file_origin`
# records which directory the rows came from, for the provenance table in the Rmd.
.read_per_annotation_dir <- function(dir, origin) {
  if (!fs::dir_exists(dir)) {
    warning("per_annotation: directory not found, skipping: ", dir)
    return(tibble::tibble())
  }
  files <- fs::dir_ls(dir, glob = "*.csv")            # non-recursive
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
        patient_id         = slide_key(m[1, 2]),
        annotation         = paste0("ANNOTATION_", as.integer(m[1, 3])),
        phenotype_clean    = stringr::str_extract(phenotype, "(?<=\\().*?(?=\\))"),
        file_origin        = origin,
        .membership_source = "flag"
      )
  })
}

# Per-annotation cells, optionally with an OVERLAY directory layered on top.
# `overlay` (e.g. PER_ANNOT_OLD_DIR) WINS on collision: for every (patient_id,
# annotation) the overlay provides, the top-level rows are dropped entirely and
# replaced — the two files are alternative exports of the same annotation, so mixing
# their cells would double-count. Annotations the overlay does not mention are kept
# from `dir` untouched. overlay = NULL (default) reproduces the top-level-only load.
load_per_annotation_cells <- function(dir = PER_ANNOT_DIR, overlay = NULL) {
  base <- .read_per_annotation_dir(dir, origin = "per_annotation")
  if (is.null(overlay)) return(base)

  over <- .read_per_annotation_dir(overlay, origin = "old")
  if (nrow(over) == 0) return(base)
  if (nrow(base) == 0) return(over)

  superseded <- dplyr::distinct(over, patient_id, annotation)
  base |>
    dplyr::anti_join(superseded, by = c("patient_id", "annotation")) |>
    dplyr::bind_rows(over) |>
    dplyr::arrange(patient_id, annotation)
}

# Which ANNOTATION_k columns of neoplastic_data a patient actually has a score for.
# Returns a named list: slide_key(patient) -> character vector of annotation names.
.scored_annotations <- function(neoplastic_data) {
  ann_cols <- grep("^ANNOTATION_", names(neoplastic_data), value = TRUE)
  stats::setNames(
    lapply(seq_len(nrow(neoplastic_data)), function(i)
      ann_cols[!is.na(unlist(neoplastic_data[i, ann_cols]))]),
    slide_key(neoplastic_data$SAMPLE)
  )
}

# Whole-slide FALLBACK cells for patients with NO per-annotation file. The slide csv
# carries a SINGLE Out_of_annotation flag, computed over all of that patient's
# annotations at once, so it yields ONE region per patient.
#
# How that region is LABELLED decides whether the patient reaches the per-annotation
# panels, because those join to the pathologist on (patient_id, annotation):
#   patient has exactly ONE scored annotation -> label it that annotation
#       (e.g. "ANNOTATION_1"). Sound: the flag's region and the pathologist's single
#       scored annotation are the same region, so the join compares like with like.
#   patient has TWO OR MORE scored annotations -> label "csv" and warn.
#       The one flag region spans ALL of them, so it cannot be attributed to any
#       single ANNOTATION_k; joining it to one score would compare a whole-slide
#       region against a part-region score. Such a patient reaches the union panels
#       only (the union score is the mean over its annotations, which does match).
# `source` stays "csv" either way, so the provenance is still visible in every table
# and the relabelling never disguises a fallback patient as a per-annotation one.
#
# `have_ids` are the slide_key ids already covered per-annotation; those patients are
# skipped so a patient is never counted twice.
load_flag_fallback_cells <- function(ihc_data, have_ids = character(),
                                     neoplastic_data = NULL) {
  if (!"Out_of_annotation" %in% names(ihc_data)) {
    warning("fallback: ihc_data has no Out_of_annotation column — no fallback patients")
    return(tibble::tibble())
  }
  out <- ihc_data |>
    dplyr::mutate(patient_id = slide_key(patient_id)) |>
    dplyr::filter(!patient_id %in% have_ids)
  if (nrow(out) == 0) return(tibble::tibble())

  scored <- if (is.null(neoplastic_data)) list() else .scored_annotations(neoplastic_data)

  # One label per fallback patient, decided by its scored-annotation count.
  labels <- vapply(unique(out$patient_id), function(pid) {
    anns <- scored[[as.character(pid)]]
    if (length(anns) == 1L) return(anns)
    if (length(anns) > 1L)
      warning("fallback: patient ", pid, " has ", length(anns), " scored annotations (",
              paste(anns, collapse = ", "), ") but only one whole-slide flag region; ",
              "labelling it \"csv\" — it will appear in the union panels only.")
    "csv"
  }, character(1))

  out |>
    dplyr::mutate(
      annotation         = unname(labels[as.character(patient_id)]),
      file_origin        = "whole_slide",
      .membership_source = "csv"
    )
}

# The full cells table both Rmds run on: per-annotation (optionally overlaid with
# old/) plus the whole-slide fallback for any patient the per-annotation exports do
# not cover. The ONLY difference between the two analyses is `overlay`.
# Pass `neoplastic_data` so a fallback patient with a single scored annotation is
# labelled with that annotation and therefore reaches the per-annotation panels too
# (see load_flag_fallback_cells); without it every fallback patient stays "csv".
build_flag_cells <- function(ihc_data, dir = PER_ANNOT_DIR, overlay = NULL,
                             neoplastic_data = NULL) {
  ann <- load_per_annotation_cells(dir, overlay = overlay)
  fb  <- load_flag_fallback_cells(ihc_data, have_ids = unique(ann$patient_id),
                                  neoplastic_data = neoplastic_data)
  dplyr::bind_rows(ann, fb)
}

# inside = the FlowPath flag says the cell is NOT out of the annotation.
.flag_inside <- function(cells) !.is_outside(cells$Out_of_annotation)

# Provenance label for a metrics row: "flag" for per-annotation cells, "csv" for the
# whole-slide fallback. Collapsed with "+" on the (by construction impossible) mix so
# a mistake shows up in the table instead of being silently hidden.
.membership_label <- function(cells) {
  if (!".membership_source" %in% names(cells)) return("flag")
  paste(sort(unique(cells$.membership_source)), collapse = "+")
}

# Cell identity for union-dedup: cell_id if the export carries it, else the centroid.
.cell_key_cols <- function(cells) {
  if ("cell_id" %in% names(cells)) "cell_id"
  else intersect(c("centroid_x", "centroid_y"), names(cells))
}

# PER-ANNOTATION metrics: one row per (patient_id, annotation). Inside cells only,
# no polygon area (area_mm2 / dens_* = NA). Same column schema as
# ihc_annotation_metrics(..., scope = "per_annotation"). Fallback patients contribute
# their single annotation = "csv" row here.
ihc_flag_metrics_per_annotation <- function(ann_cells, um_per_px = 0.325) {
  ann_cells |>
    dplyr::group_by(patient_id, annotation) |>
    dplyr::group_split() |>
    purrr::map_dfr(function(cells) {
      region_ratios_area(cells[.flag_inside(cells), , drop = FALSE], 0, um_per_px) |>
        dplyr::mutate(patient_id = cells$patient_id[1],
                      annotation = cells$annotation[1],
                      source     = .membership_label(cells), .before = 1)
    })
}

# UNION metrics: one row per patient. Pool inside cells across the patient's
# annotation files and dedup to one row per cell (so a cell inside two annotations
# counts once), then the same area-free ratios. Same schema as scope = "union". For a
# fallback patient the "union" is just its single whole-slide region.
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
                      annotation = "union",
                      source     = .membership_label(cells), .before = 1)
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

# Provenance inventory: which file fed each (patient, annotation), how many cells,
# how many the flag puts inside. `file_origin` is "per_annotation" / "old" /
# "whole_slide", so the overlay Rmd can show at a glance which annotations old/
# superseded.
flag_cells_inventory <- function(ann_cells) {
  ann_cells |>
    dplyr::group_by(patient_id, annotation, file_origin, membership = .membership_source) |>
    dplyr::summarise(n_cells  = dplyr::n(),
                     n_inside = sum(!.is_outside(Out_of_annotation)),
                     .groups  = "drop") |>
    dplyr::arrange(patient_id, annotation)
}
