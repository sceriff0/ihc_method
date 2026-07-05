# =============================================================================
# validation_helpers.R  —  shared logic for the three IHC-validation reports
#   (clinical_data.Rmd, deconvolution.Rmd, bulkRna.Rmd).
#
# The FlowPath single-cell IHC table is the "measurement under test". These
# helpers derive the quantities each report validates against an independent
# reference (pathologist / bulk-RNA deconvolution / bulk-RNA marker genes).
#
# Deliberate divergence from ATTEND's code/attend_ihc.R: "inside the annotation"
# is recomputed here from the raw pathologist geojson via sf point-in-polygon,
# NOT read from the CSV `Out_of_annotation` flag (which is not trusted for this
# project). Multiple annotations per patient are reported both per-annotation and
# as their dissolved union.
#
# renv: needs sf, jsonlite (+ tidyverse, here, fs already in the lockfile).
#   renv::install(c("sf", "jsonlite")); renv::snapshot()
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(fs)
  library(sf)
})

# --- ID handling ------------------------------------------------------------
# All datasets share one patient/slide ID but differ in punctuation (clinical
# CRF uses dots, neoplastic uses bare digits, IHC comes from a filename). Strip
# non-alphanumerics and upper-case so the shared key joins; leading zeros are
# kept (they are significant, e.g. "052").
norm_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

# TRUE where a FlowPath `<marker>_sign` column marks a positive cell. The
# FlowPath export uses "+"; the alternatives guard against encoding drift.
is_pos <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  tolower(trimws(as.character(x))) %in% c("+", "pos", "positive", "yes", "true", "1")
}

# --- Phenotype / cell-type vocabulary ---------------------------------------
# `phenotype_clean` (the parenthetical label in `phenotype`) is the cell type.
# Collapse the 12 observed labels into lineages that have a bulk-RNA counterpart.
# NOTE: there is no macrophage/B-cell phenotype in this panel, so those deconv
# cell types have no phenotype-level IHC match (see marker fractions instead).
phenotype_lineage <- tibble::tribble(
  ~phenotype_clean,            ~lineage,
  "PANCK+Tumor",               "Tumor",
  "VIM+Tumor",                 "Tumor",
  "T helper",                  "CD4T",
  "T cytotoxic",               "CD8T",
  "Activated T cytotoxic",     "CD8T",
  "CD8+ T reg",                "Treg",
  "CD4+ Treg",                 "Treg",
  "Natural Killer",            "NK",
  "Activated Natural Killer",  "NK",
  "Immune",                    "Immune_other",
  "Stroma",                    "Stroma",
  "Unknown",                   "Unknown"
)

# Lineages with a clean deconvolution counterpart, used for the method comparison.
comparable_lineages <- c("CD8T", "CD4T", "Treg", "NK")

# IHC protein marker -> canonical gene symbol(s). Multi-gene markers are averaged
# downstream. `wrongL1CAM` and DAPI are intentionally excluded.
marker_gene_map <- tibble::tribble(
  ~marker,     ~genes,
  "CD45",      "PTPRC",
  "CD3",       "CD3D,CD3E,CD3G",
  "CD8",       "CD8A,CD8B",
  "CD4",       "CD4",
  "GZMB",      "GZMB",
  "FOXP3",     "FOXP3",
  "CD14",      "CD14",
  "CD163",     "CD163",
  "CD56",      "NCAM1",
  "SMA",       "ACTA2",
  "PANCK",     "EPCAM,KRT5,KRT8,KRT14,KRT18,KRT19",
  "VIMENTIN",  "VIM",
  "CD74",      "CD74",
  "L1CAM",     "L1CAM",
  "P53",       "TP53",
  "PDL1",      "CD274",
  "PD1",       "PDCD1",
  "ARID1A",    "ARID1A",
  "FSP1",      "S100A4"
)

# Markers carrying a `<marker>_sign` column in ihc_data (wrongL1CAM excluded).
ihc_markers <- marker_gene_map$marker

# --- GeoJSON reading (sf) ----------------------------------------------------
# Parse a QuPath-style annotation geojson into an sf polygon. Handles a bare
# Feature, a FeatureCollection, or a raw geometry; closes open rings; dissolves
# multiple features in one file into a single (multi)polygon. jsonlite is used
# (simplifyVector = FALSE) so nested coordinate arrays parse predictably across
# GDAL versions.
.ring_matrix <- function(ring) {
  m <- do.call(rbind, lapply(ring, function(pt) c(as.numeric(pt[[1]]), as.numeric(pt[[2]]))))
  if (!isTRUE(all.equal(m[1, ], m[nrow(m), ]))) m <- rbind(m, m[1, ])  # close ring
  m
}

read_polygon_geojson <- function(path) {
  j <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  feats <- if (identical(j$type, "FeatureCollection")) j$features
           else if (identical(j$type, "Feature"))       list(j)
           else                                          list(list(geometry = j))

  geoms <- lapply(feats, function(f) {
    g <- f$geometry
    if (identical(g$type, "Polygon")) {
      sf::st_polygon(lapply(g$coordinates, .ring_matrix))
    } else if (identical(g$type, "MultiPolygon")) {
      sf::st_multipolygon(lapply(g$coordinates,
        function(poly) lapply(poly, .ring_matrix)))
    } else {
      stop(sprintf("unsupported geometry type: %s", g$type))
    }
  })
  # planar image pixel coordinates -> no CRS; dissolve to one geometry
  sf::st_union(sf::st_sfc(geoms))
}

# Load every annotation polygon under `dir`, keyed by patient and annotation index
# from the filename `<patient_id>_a<k>.geojson` -> annotation "ANNOTATION_<k>"
# (matching neoplastic_data's ANNOTATION_1..3 columns).
load_annotations <- function(dir = here::here("data", "annotation"),
                             patient_ids = NULL) {
  files <- fs::dir_ls(dir, glob = "*.geojson")
  if (length(files) == 0) stop(sprintf("no geojson files found in %s", dir))

  meta <- tibble::tibble(
    path = as.character(files),
    stem = fs::path_ext_remove(fs::path_file(files))
  ) |>
    tidyr::extract(stem, c("patient_id", "ann_num"), "^(.*)_a(\\d+)$", remove = FALSE)

  if (!is.null(patient_ids)) {
    keep <- norm_id(meta$patient_id) %in% norm_id(patient_ids)
    meta <- meta[keep, , drop = FALSE]
  }

  polys <- purrr::pmap(meta, function(path, stem, patient_id, ann_num) {
    sf::st_sf(
      patient_id = patient_id,
      annotation = paste0("ANNOTATION_", ann_num),
      geometry   = read_polygon_geojson(path)
    )
  })
  do.call(rbind, polys)  # sf provides rbind(); preserves geometry across versions
}

# --- Cell-in-annotation metrics ---------------------------------------------
# Ratio metrics for one set of cells that fall inside a region (mirrors ATTEND's
# tumor_over_all / cd45_over_inside plus lineage composition). Returns one row.
region_ratios <- function(cells) {
  n_inside <- nrow(cells)
  if (n_inside == 0) {
    return(tibble::tibble(n_inside = 0L, tumor_over_inside = NA_real_,
                          cd45_over_inside = NA_real_, frac_CD8T = NA_real_,
                          frac_CD4T = NA_real_, frac_Treg = NA_real_,
                          frac_NK = NA_real_, frac_Stroma = NA_real_))
  }
  pc       <- tidyr::replace_na(cells$phenotype_clean, "")
  is_tumor <- stringr::str_detect(pc, "Tumor")
  is_cd45  <- is_pos(cells$CD45_sign)
  lin      <- dplyr::left_join(tibble::tibble(phenotype_clean = cells$phenotype_clean),
                               phenotype_lineage, by = "phenotype_clean")$lineage
  frac <- function(l) sum(lin == l, na.rm = TRUE) / n_inside

  tibble::tibble(
    n_inside          = n_inside,
    tumor_over_inside = sum(is_tumor, na.rm = TRUE) / n_inside,
    cd45_over_inside  = sum(is_cd45,  na.rm = TRUE) / n_inside,
    frac_CD8T = frac("CD8T"), frac_CD4T = frac("CD4T"),
    frac_Treg = frac("Treg"), frac_NK   = frac("NK"),
    frac_Stroma = frac("Stroma")
  )
}

# Point-in-polygon assignment + region metrics, per patient. `scope`:
#   "per_annotation" -> one row per (patient, annotation polygon)
#   "union"          -> one row per patient over the dissolved annotations
# Cells are matched to their own patient's polygons only. Uses sf::st_within on
# the raw geojson; the CSV Out_of_annotation flag is ignored by design.
ihc_annotation_metrics <- function(ihc_data, annots,
                                    scope = c("per_annotation", "union")) {
  scope <- match.arg(scope)
  ihc_data <- dplyr::mutate(ihc_data, .pid = norm_id(patient_id))
  annots   <- dplyr::mutate(annots,   .pid = norm_id(patient_id))
  common   <- intersect(unique(ihc_data$.pid), unique(annots$.pid))
  if (length(common) == 0) {
    warning("no patient IDs shared between ihc_data and annotations")
    return(tibble::tibble())
  }

  purrr::map_dfr(common, function(pid) {
    cells_p <- dplyr::filter(ihc_data, .pid == pid)
    polys_p <- dplyr::filter(annots,   .pid == pid)
    pts <- sf::st_as_sf(cells_p, coords = c("centroid_x", "centroid_y"),
                        remove = FALSE, crs = sf::st_crs(polys_p))
    within <- sf::st_within(pts, polys_p)  # per-cell list of polygon indices

    if (scope == "union") {
      inside <- lengths(within) > 0
      region_ratios(cells_p[inside, , drop = FALSE]) |>
        dplyr::mutate(patient_id = pid, annotation = "union", .before = 1)
    } else {
      purrr::map_dfr(seq_len(nrow(polys_p)), function(i) {
        inside_i <- vapply(within, function(idx) i %in% idx, logical(1))
        region_ratios(cells_p[inside_i, , drop = FALSE]) |>
          dplyr::mutate(patient_id = pid, annotation = polys_p$annotation[i], .before = 1)
      })
    }
  })
}

# --- Whole-slide IHC quantities (for the RNA comparisons) -------------------
# Per-patient fraction of cells positive for each marker (`<marker>_sign == "+"`)
# plus the mean z-score as a sensitivity readout. Whole slide (bulk RNA has no
# annotation boundary). Returns wide: patient_id + <marker>_posfrac + <marker>_z.
ihc_marker_fraction <- function(ihc_data, markers = ihc_markers) {
  sign_cols <- paste0(markers, "_sign")
  z_cols    <- paste0(markers, "_zscore")
  ihc_data |>
    dplyr::mutate(patient_id = norm_id(patient_id)) |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(sign_cols), ~ mean(is_pos(.x), na.rm = TRUE),
                    .names = "{.col}"),
      dplyr::across(dplyr::any_of(z_cols),   ~ mean(.x, na.rm = TRUE),
                    .names = "{.col}"),
      .groups = "drop"
    ) |>
    dplyr::rename_with(~ sub("_sign$",   "_posfrac", .x), dplyr::ends_with("_sign")) |>
    dplyr::rename_with(~ sub("_zscore$", "_z",       .x), dplyr::ends_with("_zscore"))
}

# Per-patient lineage composition over ALL cells (whole slide), for deconvolution
# comparison. Long: patient_id, lineage, n, frac.
ihc_lineage_fraction <- function(ihc_data) {
  ihc_data |>
    dplyr::mutate(patient_id = norm_id(patient_id)) |>
    dplyr::left_join(phenotype_lineage, by = "phenotype_clean") |>
    dplyr::count(patient_id, lineage, name = "n") |>
    dplyr::group_by(patient_id) |>
    dplyr::mutate(frac = n / sum(n)) |>
    dplyr::ungroup()
}

# Map an immunedeconv cell_type string to one of `comparable_lineages` (or NA).
# Pattern-based because the exact strings differ across methods.
deconv_to_lineage <- function(cell_type) {
  ct <- tolower(cell_type)
  dplyr::case_when(
    stringr::str_detect(ct, "regulatory|treg")        ~ "Treg",
    stringr::str_detect(ct, "cd8")                    ~ "CD8T",
    stringr::str_detect(ct, "cd4")                    ~ "CD4T",
    stringr::str_detect(ct, "^nk|nk cell|natural killer") ~ "NK",
    TRUE                                              ~ NA_character_
  )
}

# Spearman rho + p + n for two paired numeric vectors, as a one-row tibble.
# Guards the small-n / zero-variance cases that abort cor.test().
paired_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3 || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(tibble::tibble(n = length(x), rho = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman"))
  tibble::tibble(n = length(x), rho = unname(ct$estimate), p = ct$p.value)
}
