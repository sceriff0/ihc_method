# =============================================================================
# validation_helpers.R  —  shared logic for the three IHC-validation reports
#   (clinical_data.Rmd, deconvolution.Rmd, bulkRna.Rmd).
#
# The FlowPath single-cell IHC table is the "measurement under test". These
# helpers derive the quantities each report validates against an independent
# reference (pathologist / bulk-RNA deconvolution / bulk-RNA marker genes).
#
# "Inside the annotation" is read from the FlowPath `Out_of_annotation` flag: the
# export now writes one csv per (patient, annotation), so the flag is precomputed
# per ANNOTATION_<k> and point-in-polygon is no longer needed for membership.
# Multiple annotations per patient are reported both per-annotation and as their
# union (a cell counts if inside ANY of the patient's annotations).
#
# sf/geojson survive ONLY for the invasive-margin (periphery) metrics, where a
# distance band must be buffered from the real polygon geometry. There the cell
# centroids (microns) are mapped into the geojson pixel frame by centroid / 0.325.
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

# --- House plot style (adapted from mirage benchmarks/analysis/plots.R) ------
# Okabe-Ito colourblind-safe categorical palette.
oi <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
        "#E69F00", "#56B4E9", "#F0E442", "#000000")

# Semantic colours for the immune "hot/cold" phenotype: HOT -> red, COLD -> light
# blue, intermediate/other -> orange/grey. Robust to case/spelling. Returns a
# named vector keyed by the given levels, for scale_colour_manual/scale_fill_manual.
hotcold_cols <- function(levels) {
  lv  <- as.character(levels)
  key <- toupper(trimws(lv))
  col <- dplyr::case_when(
    grepl("HOT|INFLAM", key)               ~ "#D7191C",  # red
    grepl("COLD|DESERT", key)              ~ "#74ADD1",  # light blue
    grepl("INTERMED|VARI|MIX|EXCLUD", key) ~ "#FDAE61",  # orange
    TRUE                                   ~ "grey65"
  )
  stats::setNames(col, lv)
}

# Order immune-phenotype levels cold -> intermediate -> hot (axis/legend order).
hotcold_order <- function(x) {
  lv   <- unique(as.character(x[!is.na(x)]))
  key  <- toupper(trimws(lv))
  rank <- ifelse(grepl("COLD|DESERT", key), 1L,
          ifelse(grepl("INTERMED|VARI|MIX|EXCLUD", key), 2L,
          ifelse(grepl("HOT|INFLAM", key), 3L, 4L)))
  factor(x, levels = lv[order(rank)])
}

# Publication theme: generous type, restrained gridlines, bold titles, grey
# subtitles/captions, top-left legend. Apply per-Rmd with theme_set(theme_paper).
theme_paper <- ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.05)),
    plot.subtitle    = ggplot2::element_text(colour = "grey35",
                                             margin = ggplot2::margin(b = 8)),
    plot.caption     = ggplot2::element_text(colour = "grey55",
                                             size = ggplot2::rel(.7), hjust = 1),
    plot.title.position = "plot", plot.caption.position = "plot",
    axis.title       = ggplot2::element_text(colour = "grey20"),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(linewidth = .3, colour = "grey90"),
    strip.text       = ggplot2::element_text(face = "bold"),
    legend.position  = "top", legend.justification = "left",
    plot.margin      = ggplot2::margin(12, 16, 8, 12))

# --- ID handling ------------------------------------------------------------
# All datasets share one patient/slide ID but differ in punctuation (clinical
# CRF uses dots, neoplastic uses bare digits, IHC comes from a filename). Strip
# non-alphanumerics and upper-case so the shared key joins; leading zeros are
# kept (they are significant, e.g. "052").
norm_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

# Normalise a SLIDE-space id (IHC patient_id, neoplastic SAMPLE, clinical
# `ID PATIENT`). These are numeric slide codes that sometimes carry an alpha
# prefix (e.g. "EPM - 052"); keep only the digits so "EPM - 052" -> "052" matches
# the IHC "052". Falls back to norm_id() for any purely non-numeric id.
norm_slide_id <- function(x) {
  x <- as.character(x)
  d <- gsub("[^0-9]", "", x)
  ifelse(d == "" | is.na(x), norm_id(x), d)
}

# Known clinical `ID PATIENT` data-entry errors: the value on the LEFT is what the
# clinical sheet reads, the RIGHT is the true slide/neoplastic id it should match.
# Default handles the verified transposition (clinical "15879" == slide "15897").
SLIDE_ID_FIXES <- c("15879" = "15897")

# THE canonical slide-space join key: normalise + apply the typo fixes. Every
# join on a slide id (IHC, neoplastic, clinical `ID PATIENT`) must use this so no
# call site silently omits the correction and drops a patient.
slide_key <- function(x, fixes = SLIDE_ID_FIXES) {
  s   <- norm_slide_id(x)
  hit <- s %in% names(fixes)
  s[hit] <- unname(fixes[s[hit]])
  s
}

# Crosswalk between the two id systems, read from clinical_data:
#   slide_id (norm_slide_id of `ID PATIENT`)  <->  crf_id (norm_id of `ID CRF PRESERVE`)
# `ID CRF PRESERVE` also keys counts_data$Sample, so this bridges IHC/neoplastic
# (slide space) to bulk RNA (CRF space). `slide_id_fixes` corrects known clinical
# data-entry errors: the default handles the verified transposition where clinical
# reads "15879" but the slide/neoplastic id is "15897". Extend or clear as needed.
id_crosswalk <- function(clinical_data,
                         patient_col = "ID PATIENT",
                         crf_col     = "ID CRF PRESERVE",
                         slide_id_fixes = SLIDE_ID_FIXES) {
  slide <- slide_key(clinical_data[[patient_col]], slide_id_fixes)
  crf   <- norm_id(clinical_data[[crf_col]])
  keep  <- slide != "" & crf != "" & !is.na(slide) & !is.na(crf)
  unique(tibble::tibble(slide_id = slide[keep], crf_id = crf[keep]))
}

# Named vector  crf_id -> slide_id  built DIRECTLY from clinical_data (base R: no
# tibble, no join). Map bulk-RNA (CRF space) sample ids to IHC slide ids with
# `slide[crf_id]`. Prefer this over id_crosswalk() + join in the RNA reports —
# dplyr join column-propagation is unreliable on this project's dplyr build.
crf_to_slide_map <- function(clinical_data,
                             patient_col = "ID PATIENT",
                             crf_col     = "ID CRF PRESERVE",
                             slide_id_fixes = SLIDE_ID_FIXES) {
  slide <- slide_key(clinical_data[[patient_col]], slide_id_fixes)
  crf   <- norm_id(clinical_data[[crf_col]])
  keep  <- !is.na(slide) & !is.na(crf) & slide != "" & crf != ""
  stats::setNames(slide[keep], crf[keep])
}

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
  # planar image pixel coordinates -> no CRS; dissolve to one geometry.
  # st_make_valid: QuPath annotation rings are often self-intersecting, and
  # st_within() silently returns NO matches on an invalid polygon (every cell
  # reads as outside). Repairing the geometry is what makes point-in-polygon work.
  sf::st_make_valid(sf::st_union(sf::st_sfc(geoms)))
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
# Lineages tracked in every region (immune subsets + stroma).
region_lineages <- c("CD8T", "CD4T", "Treg", "NK", "Immune_other", "Stroma")

# Per-region counts + ratios for one set of cells inside a region. Emits the three
# ATTEND denominators as raw counts (n_inside / n_tumor_inside / n_cd45_inside) and
# per-lineage counts (n_<lineage>), so composition can be normalised any of ATTEND's
# three ways downstream (see region_composition()). `frac_<lineage>` is the default
# normalisation (per all cells inside). Returns one row.
region_ratios <- function(cells) {
  n_inside <- nrow(cells)
  is_tumor <- if (n_inside) stringr::str_detect(tidyr::replace_na(cells$phenotype_clean, ""), "Tumor") else logical(0)
  is_cd45  <- if (n_inside) is_pos(cells$CD45_sign) else logical(0)
  lin      <- if (n_inside)
    dplyr::left_join(tibble::tibble(phenotype_clean = cells$phenotype_clean),
                     phenotype_lineage, by = "phenotype_clean")$lineage else character(0)
  n_tumor  <- sum(is_tumor, na.rm = TRUE)
  n_cd45   <- sum(is_cd45,  na.rm = TRUE)
  ncount   <- function(l) sum(lin == l, na.rm = TRUE)
  safe     <- function(num, den) if (den > 0) num / den else NA_real_

  out <- tibble::tibble(
    n_inside          = n_inside,
    n_tumor_inside    = n_tumor,
    n_cd45_inside     = n_cd45,
    tumor_over_inside = safe(n_tumor, n_inside),
    cd45_over_inside  = safe(n_cd45,  n_inside)
  )
  for (l in region_lineages) out[[paste0("n_", l)]]    <- ncount(l)
  for (l in region_lineages) out[[paste0("frac_", l)]] <- safe(ncount(l), n_inside)
  out
}

# ATTEND-style multi-normalisation composition (mirrors code/attend_ihc.R
# ihc_celltype_metrics). From a region-metrics table (rows from region_ratios,
# carrying patient_id/n_* counts), returns LONG per lineage with the same three
# denominators ATTEND uses:
#   frac_inside = lineage / all cells inside   (overall composition)
#   frac_tumor  = lineage / tumour cells inside (immune-to-tumour density)
#   frac_cd45   = lineage / CD45+ cells inside  (composition of the immune compartment)
region_composition <- function(region_metrics,
                               lineages = c("CD8T", "CD4T", "Treg", "NK", "Immune_other")) {
  id_cols <- intersect(c("patient_id", "annotation", "source"), names(region_metrics))
  region_metrics |>
    dplyr::select(dplyr::all_of(c(id_cols, "n_inside", "n_tumor_inside", "n_cd45_inside")),
                  dplyr::all_of(paste0("n_", lineages))) |>
    tidyr::pivot_longer(dplyr::all_of(paste0("n_", lineages)),
                        names_to = "lineage", names_prefix = "n_", values_to = "n_cell") |>
    dplyr::mutate(
      frac_inside = dplyr::if_else(n_inside       > 0, n_cell / n_inside,       NA_real_),
      frac_tumor  = dplyr::if_else(n_tumor_inside > 0, n_cell / n_tumor_inside, NA_real_),
      frac_cd45   = dplyr::if_else(n_cd45_inside  > 0, n_cell / n_cd45_inside,  NA_real_)
    )
}

# TRUE where a FlowPath `Out_of_annotation` value marks a cell OUTSIDE the region.
.is_outside <- function(x) tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "t")

# Stable per-cell key within a patient, for OR-ing membership across a patient's
# annotation files. Prefers the real `cell_id`; falls back to the loader's
# within-file row index `.cell_row` (the files are row-aligned, identical cells).
.cell_key <- function(df) {
  if ("cell_id" %in% names(df)) df$cell_id else df$.cell_row
}

# Region metrics per patient from the FlowPath `Out_of_annotation` flag (no sf).
# `ihc_annot` is the LONG table: one row per cell per annotation file, carrying
# `patient_id`, `annotation` (ANNOTATION_<k>) and `Out_of_annotation`.
#   scope = "per_annotation": one row per (patient, annotation); a cell is inside
#           iff !Out_of_annotation in that annotation's file.
#   scope = "union":          one row per patient; a cell is inside iff it is
#           inside ANY of the patient's annotations (OR over `cell_id`).
# `source` is always "csv" (the flag is the sole membership source now).
ihc_annotation_metrics <- function(ihc_annot,
                                    scope = c("per_annotation", "union")) {
  scope    <- match.arg(scope)
  ihc_annot <- dplyr::mutate(ihc_annot, .pid = slide_key(patient_id),
                             .inside = !.is_outside(Out_of_annotation),
                             .key = .cell_key(ihc_annot))

  purrr::map_dfr(unique(ihc_annot$.pid), function(pid) {
    cells_p <- dplyr::filter(ihc_annot, .pid == pid)

    if (scope == "per_annotation") {
      purrr::map_dfr(sort(unique(cells_p$annotation)), function(ann) {
        ca <- dplyr::filter(cells_p, annotation == ann)
        region_ratios(ca[ca$.inside, , drop = FALSE]) |>
          dplyr::mutate(patient_id = pid, annotation = ann, source = "csv", .before = 1)
      })
    } else {
      # union: a cell is inside if inside in any annotation; take ONE copy of the
      # cell (lowest-index file) so counts are not multiplied by the file count.
      inside_key <- unique(cells_p$.key[cells_p$.inside])
      base       <- dplyr::filter(cells_p, ann_num == min(ann_num))
      region_ratios(base[base$.key %in% inside_key, , drop = FALSE]) |>
        dplyr::mutate(patient_id = pid, annotation = "union", source = "csv", .before = 1)
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
    dplyr::mutate(patient_id = slide_key(patient_id)) |>
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

# Long form of the per-cell marker table: one row per (cell, marker) with the
# raw / zscore / sign triple side by side. Used by the internal-QC benchmarks
# (phenotype-marker concordance, per-channel usability). wrongL1CAM/DAPI excluded
# via the default marker list.
ihc_marker_long <- function(ihc_data, markers = ihc_markers) {
  keep <- intersect(c("cell_id", "patient_id", "phenotype_clean"), names(ihc_data))
  trip <- c(paste0(markers, "_raw"), paste0(markers, "_zscore"), paste0(markers, "_sign"))
  ihc_data |>
    dplyr::select(dplyr::any_of(keep), dplyr::any_of(trip)) |>
    tidyr::pivot_longer(
      -dplyr::any_of(keep),
      names_to     = c("marker", ".value"),
      names_pattern = "(.*)_(raw|zscore|sign)"
    ) |>
    dplyr::mutate(patient_id = slide_key(patient_id), pos = is_pos(sign))
}

# Per-patient lineage composition over ALL cells (whole slide), for deconvolution
# comparison. Long: patient_id, lineage, n, frac.
ihc_lineage_fraction <- function(ihc_data) {
  ihc_data |>
    dplyr::mutate(patient_id = slide_key(patient_id)) |>
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

# =============================================================================
# Tumour border / periphery (invasive-margin) metrics
# -----------------------------------------------------------------------------
# The tumour annotation gives an "inside vs outside" partition; the biology at
# the tumour-host interface needs a third region — the invasive margin (IM), a
# band centred on the annotation boundary. Literature grounding for the band
# half-width `d` (microns each side of the border):
#   * Consensus Immunoscore (Galon 2014 J Pathol; Pages 2018 Lancet) scores CD3/
#     CD8 densities in the tumour core (CT) and a 500 um invasive margin (IM) —
#     the de-facto standard, and the default in QuPath IM workflows.
#   * Reviews put the IM width in a 200-500 um range, with some invasive-margin
#     detection algorithms using up to 1000 um.
#   * Immune spatial phenotypes (Chen & Mellman; Hegde 2016) are defined by WHERE
#     CD8 sits relative to this border: inflamed/hot = throughout the core,
#     excluded/cold = trapped at the margin, desert/cold = sparse everywhere.
# So we report several thresholds (default 100 / 250 / 500 um) and, per threshold,
# both the margin band and the eroded "deep core" so a margin-vs-core contrast can
# recover the inflamed/excluded/desert axis.
#
# All geometry is done in the polygon's PIXEL coordinate frame. Cell centroids are
# in microns and the QuPath geojson is in pixels, so centroids are mapped in by
# dividing by `um_per_px` (= 0.325): x_px = centroid_x / 0.325. One pixel is then
# `um_per_px` microns, so a micron threshold `d` becomes `d / um_per_px` pixel
# units, and a pixel area is scaled by `um_per_px^2` to reach microns^2.

# From a dissolved core polygon (sfc) and a signed buffer distance `d` (polygon
# units), return the invasive-margin band (+/- d around the boundary) and the
# eroded interior "deep core" (further than d inside the boundary). If d exceeds
# the core's inradius the eroded core is empty (sfc of length 0) — the caller
# treats that as "no deep core at this threshold".
margin_regions <- function(core, d) {
  outer <- suppressWarnings(sf::st_buffer(core, d))
  inner <- suppressWarnings(sf::st_buffer(core, -d))
  band  <- if (length(inner) == 0 || all(sf::st_is_empty(inner))) outer
           else suppressWarnings(sf::st_difference(outer, inner))
  list(margin = band, core = inner)
}

# region_ratios() augmented with region AREA and area-normalised DENSITIES
# (cells / mm^2), the native Immunoscore unit. `area_units2` is the region area in
# squared polygon units; `um_per_unit` converts it to mm^2. Densities are NA when
# the region is empty/degenerate so they never masquerade as a real zero.
region_ratios_area <- function(cells, area_units2, um_per_unit) {
  rr       <- region_ratios(cells)
  area_mm2 <- area_units2 * (um_per_unit^2) / 1e6         # units^2 -> um^2 -> mm^2
  dens     <- function(n) if (is.finite(area_mm2) && area_mm2 > 0) n / area_mm2 else NA_real_
  rr$area_mm2          <- if (area_mm2 > 0) area_mm2 else NA_real_
  rr$dens_all_per_mm2  <- dens(rr$n_inside)
  rr$dens_cd45_per_mm2 <- dens(rr$n_cd45_inside)
  for (l in c("CD8T", "CD4T", "Treg", "NK"))
    rr[[paste0("dens_", l, "_per_mm2")]] <- dens(rr[[paste0("n_", l)]])
  rr
}

# Invasive-margin metrics per patient. For every requested `thresholds_um` and,
# per `scope`, either the dissolved union (one core) or each single annotation
# (one core each), returns TWO rows — region = "margin" (the +/- d band) and
# region = "core" (interior beyond d) — carrying every region_ratios_area column.
# Only patients WITH a geojson polygon are covered: a margin band cannot be built
# from the binary CSV Out_of_annotation flag, so there is no CSV fallback here.
# `um_per_unit` and `source` are echoed so the physical band width is auditable.
ihc_periphery_metrics <- function(ihc_data, annots,
                                  scope = c("union", "per_annotation"),
                                  thresholds_um = c(100, 250, 500),
                                  um_per_px = 0.325) {
  scope    <- match.arg(scope)
  ihc_data <- dplyr::mutate(ihc_data, .pid = slide_key(patient_id))
  annots   <- dplyr::mutate(annots,   .pid = slide_key(patient_id))
  ann_ids  <- unique(annots$.pid)

  purrr::map_dfr(intersect(unique(ihc_data$.pid), ann_ids), function(pid) {
    cells_p <- dplyr::filter(ihc_data, .pid == pid,
                             is.finite(centroid_x), is.finite(centroid_y))
    if (nrow(cells_p) == 0) return(tibble::tibble())
    polys_p <- dplyr::filter(annots, .pid == pid)
    crs     <- sf::st_crs(polys_p)
    # Fixed µm -> pixel mapping (centroids in microns, geojson in pixels): / 0.325.
    upu     <- um_per_px    # microns per polygon (pixel) unit
    pts     <- sf::st_as_sf(data.frame(x = cells_p$centroid_x / um_per_px,
                                       y = cells_p$centroid_y / um_per_px),
                            coords = c("x", "y"), crs = crs)

    cores <- if (scope == "union") {
      list(union = sf::st_union(sf::st_geometry(polys_p)))
    } else {
      stats::setNames(lapply(seq_len(nrow(polys_p)),
                             function(i) sf::st_geometry(polys_p)[i]),
                      polys_p$annotation)
    }

    purrr::imap_dfr(cores, function(core, ann_label) {
      purrr::map_dfr(thresholds_um, function(d_um) {
        rg <- margin_regions(core, d_um / upu)
        row_for <- function(region_geom, region_label) {
          empty  <- length(region_geom) == 0 || all(sf::st_is_empty(region_geom))
          inside <- if (empty) rep(FALSE, nrow(cells_p))
                    else lengths(sf::st_within(pts, region_geom)) > 0
          area_u2 <- if (empty) 0 else as.numeric(sum(sf::st_area(region_geom)))
          region_ratios_area(cells_p[inside, , drop = FALSE], area_u2, upu) |>
            dplyr::mutate(region = region_label, .before = 1)
        }
        dplyr::bind_rows(row_for(rg$margin, "margin"),
                         row_for(rg$core,   "core")) |>
          dplyr::mutate(patient_id = pid, annotation = ann_label,
                        threshold_um = d_um, um_per_unit = upu,
                        source = "sf:um2px(/0.325)", .before = 1)
      })
    })
  })
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

# Pairwise cross-source agreement, for comparing several deconvolution methods (and
# IHC) against each other rather than only against one reference. `df` is long with
# columns source / patient_id / lineage / value (one value per source x lineage x
# patient). For every ordered source pair it correlates each lineage across the
# shared patients (Spearman, rank-based so differing score scales are fine) and
# averages the per-lineage rho, returning a symmetric long table (a == b -> rho 1).
# Different methods put scores on different scales; only rank agreement is meaningful.
pairwise_agreement <- function(df, sources = sort(unique(df$source))) {
  grid <- expand.grid(a = sources, b = sources, stringsAsFactors = FALSE)
  purrr::pmap_dfr(grid, function(a, b) {
    if (a == b)
      return(tibble::tibble(a = a, b = b, mean_rho = 1, n_lineages = NA_integer_))
    j <- dplyr::inner_join(
      dplyr::filter(df, source == a),
      dplyr::filter(df, source == b),
      by = c("patient_id", "lineage"), suffix = c("_a", "_b"))
    if (nrow(j) == 0)
      return(tibble::tibble(a = a, b = b, mean_rho = NA_real_, n_lineages = 0L))
    per <- j |>
      dplyr::group_by(lineage) |>
      dplyr::group_modify(~ paired_spearman(.x$value_a, .x$value_b)) |>
      dplyr::ungroup()
    tibble::tibble(a = a, b = b,
                   mean_rho   = mean(per$rho, na.rm = TRUE),
                   n_lineages = sum(is.finite(per$rho)))
  })
}
