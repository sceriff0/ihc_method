# ============================================================================
# aggregation_compare.R — how the annotation-aggregation choice moves the
# IHC-vs-pathologist concordance
# ============================================================================
# A patient can carry SEVERAL annotations, so both sides of the neoplastic-
# cellularity comparison have to be collapsed to one number per patient before they
# can be paired — and the collapse is a free choice that changes the answer.
#
# THE WEIGHTING MISMATCH this module exists to expose
#   The Rmds' existing "union" pair collapses the two sides INCONSISTENTLY:
#     pathologist  path_union = mean(ANNOTATION_k / 100)      <- UNWEIGHTED mean
#     IHC          ihc_tumor_union = tumour cells / all cells
#                    over the POOLED cells of every annotation <- CELL-WEIGHTED mean
#   For a patient whose annotations differ in size those are different summaries of
#   the same slide, so part of the union gap is aggregation artefact rather than
#   IHC-vs-pathologist disagreement. The grid below pairs every IHC aggregator with
#   every pathologist aggregator so that artefact is visible and quantified.
#
# MATCHED pairs — the two combinations that weight both sides the same way, and so
# are the only ones whose gap is interpretable as pure disagreement:
#   (pooled,   wmean_cells)  both cell-weighted
#   (mean_ann, mean)         both annotation-unweighted
# `.is_matched()` flags them; the Rmds mark them in every table and figure.
#
# CONSISTENCY CHECK built into the grid: `wmean_ann` (cell-weighted mean of the
# per-annotation fractions) should reproduce `pooled` (the union) almost exactly,
# because pooling cells IS cell-weighting. The only thing separating them is the
# union's dedup of cells shared by overlapping annotations — so a visible gap
# between those two rows measures annotation OVERLAP, not an aggregation error.
#
# Depends on validation_helpers.R (paired_cor3) — source THAT first.
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# How each side collapses a patient's annotations. Names are the codes used in the
# frames; values are the labels shown on axes, facets and tables.
IHC_AGG_LABELS <- c(
  pooled     = "pooled cells (union)",
  wmean_ann  = "cell-weighted mean of annotations",
  mean_ann   = "mean of annotations",
  median_ann = "median of annotations"
)

PATH_AGG_LABELS <- c(
  wmean_cells = "cell-weighted mean of ANNOTATION_k",
  mean        = "mean of ANNOTATION_k",
  median      = "median of ANNOTATION_k",
  max         = "max of ANNOTATION_k",
  min         = "min of ANNOTATION_k"
)

# The two combinations that weight both sides identically (see header).
.is_matched <- function(ihc_agg, path_agg) {
  (ihc_agg == "pooled"   & path_agg == "wmean_cells") |
  (ihc_agg == "mean_ann" & path_agg == "mean")
}

# The cell count that acts as the weight for a given metric: a metric computed over
# the CLEAN cell set must be weighted by the CLEAN denominator, otherwise annotations
# whose cells were mostly dropped as outliers/Unknown would be over-weighted.
.weight_col <- function(metric) {
  if (identical(metric, "tumor_over_inside_clean")) "n_inside_clean" else "n_inside"
}

# Weighted mean over the finite (value, weight) pairs only; NA when nothing is left
# or the weights sum to zero (an annotation set with no cells carries no information).
.wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

# ---------------------------------------------------------------------------
# IHC side: one value per (patient, aggregator) for one metric.
#   pooled      taken from the union frame (cells pooled and deduped, then the ratio
#               recomputed) — this is what the Rmds' union panels already use
#   wmean_ann   per-annotation fractions averaged with the cell counts as weights
#   mean_ann    per-annotation fractions averaged, every annotation equal
#   median_ann  median of the per-annotation fractions (outlier-resistant; at 2
#               annotations it equals mean_ann, so it only bites at 3+)
# ---------------------------------------------------------------------------
aggregate_ihc <- function(per_ann, union, metric) {
  w <- .weight_col(metric)
  stopifnot(metric %in% names(per_ann), w %in% names(per_ann))

  from_ann <- per_ann |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(
      n_ann      = sum(is.finite(.data[[metric]])),
      wmean_ann  = .wmean(.data[[metric]], .data[[w]]),
      mean_ann   = mean(.data[[metric]][is.finite(.data[[metric]])]),
      median_ann = stats::median(.data[[metric]][is.finite(.data[[metric]])]),
      .groups    = "drop"
    ) |>
    # mean()/median() of an empty vector give NaN; keep the frame honest with NA.
    dplyr::mutate(dplyr::across(c(mean_ann, median_ann),
                                ~ dplyr::if_else(is.finite(.x), .x, NA_real_)))

  from_union <- union |>
    dplyr::select(patient_id, pooled = dplyr::all_of(metric)) |>
    dplyr::distinct(patient_id, .keep_all = TRUE)

  from_ann |>
    dplyr::full_join(from_union, by = "patient_id") |>
    tidyr::pivot_longer(dplyr::any_of(names(IHC_AGG_LABELS)),
                        names_to = "ihc_agg", values_to = "ihc_val")
}

# ---------------------------------------------------------------------------
# Pathologist side: one value per (patient, aggregator).
#   wmean_cells  ANNOTATION_k scores averaged with that annotation's IHC cell count
#                as the weight — the pathologist analogue of pooling the cells, and
#                the only pathologist summary that matches `pooled`'s weighting
#   mean         plain mean over the patient's annotations (what the Rmds use today)
#   median / max / min  order statistics; max is the "worst area governs" reading a
#                pathologist may intend, min the most conservative
# `per_ann` supplies the weights, joined on (patient_id, annotation) — so an
# annotation the IHC never measured contributes no weight.
# ---------------------------------------------------------------------------
aggregate_pathologist <- function(path_long, per_ann, metric) {
  w <- .weight_col(metric)

  weighted <- path_long |>
    dplyr::left_join(dplyr::select(per_ann, patient_id, annotation,
                                   .w = dplyr::all_of(w)),
                     by = c("patient_id", "annotation")) |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(wmean_cells = .wmean(path_frac, .w), .groups = "drop")

  path_long |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(
      n_path = sum(is.finite(path_frac)),
      mean   = mean(path_frac,          na.rm = TRUE),
      median = stats::median(path_frac, na.rm = TRUE),
      max    = max(path_frac,           na.rm = TRUE),
      min    = min(path_frac,           na.rm = TRUE),
      .groups = "drop"
    ) |>
    # max()/min() of an all-NA vector warn and return +/-Inf; normalise to NA.
    dplyr::mutate(dplyr::across(c(mean, median, max, min),
                                ~ dplyr::if_else(is.finite(.x), .x, NA_real_))) |>
    dplyr::left_join(weighted, by = "patient_id") |>
    tidyr::pivot_longer(dplyr::any_of(names(PATH_AGG_LABELS)),
                        names_to = "path_agg", values_to = "path_frac")
}

# ---------------------------------------------------------------------------
# The full grid: every IHC aggregator x every pathologist aggregator, for every
# metric. Returns
#   $pairs  one row per (metric, ihc_agg, path_agg, patient) — the points plotted
#   $stats  one row per (metric, ihc_agg, path_agg) — paired_cor3 + n + bias, where
#           bias = mean(IHC - pathologist), the systematic offset the correlation
#           coefficients are deliberately blind to
# ---------------------------------------------------------------------------
aggregation_grid <- function(per_ann, union, path_long,
                             metrics = c("tumor_over_inside", "tumor_over_inside_clean")) {
  pairs <- purrr::map_dfr(metrics, function(the_metric) {
    ihc  <- aggregate_ihc(per_ann, union, the_metric)
    path <- aggregate_pathologist(path_long, per_ann, the_metric)

    # The grid IS a cartesian product of the two aggregator sets, but it is built
    # one combination at a time so each join stays ONE-TO-ONE on patient_id. Joining
    # the two long frames directly would be a many-to-many join, which errors on
    # dplyr < 1.1 (no `relationship` argument) and warns on >= 1.1 — and dplyr is not
    # pinned in renv.lock, so the version on the machine is not knowable from here.
    combos <- tidyr::expand_grid(
      ia = intersect(names(IHC_AGG_LABELS),  unique(ihc$ihc_agg)),
      pa = intersect(names(PATH_AGG_LABELS), unique(path$path_agg))
    )
    purrr::pmap_dfr(combos, function(ia, pa) {
      i <- ihc  |> dplyr::filter(.data$ihc_agg  == ia) |>
        dplyr::select(patient_id, ihc_val)
      p <- path |> dplyr::filter(.data$path_agg == pa) |>
        dplyr::select(patient_id, path_frac)
      dplyr::inner_join(i, p, by = "patient_id") |>
        dplyr::filter(is.finite(ihc_val), is.finite(path_frac)) |>
        dplyr::mutate(metric = the_metric, ihc_agg = ia, path_agg = pa, .before = 1)
    })
  })
  if (nrow(pairs) == 0) return(list(pairs = pairs, stats = tibble::tibble()))

  stats <- pairs |>
    dplyr::group_by(metric, ihc_agg, path_agg) |>
    dplyr::group_modify(~ {
      cc <- paired_cor3(.x$ihc_val, .x$path_frac)
      dplyr::mutate(cc, bias = mean(.x$ihc_val - .x$path_frac))
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(matched = .is_matched(ihc_agg, path_agg)) |>
    dplyr::arrange(metric, dplyr::desc(spearman))

  list(pairs = .label_aggs(pairs), stats = .label_aggs(stats))
}

# Turn the aggregator codes into ordered, human-readable factors so facets and tiles
# come out in a deliberate order (most-defensible weighting first) rather than
# alphabetically.
.label_aggs <- function(df) {
  df |>
    dplyr::mutate(
      ihc_lab  = factor(IHC_AGG_LABELS[ihc_agg],   levels = IHC_AGG_LABELS),
      path_lab = factor(PATH_AGG_LABELS[path_agg], levels = PATH_AGG_LABELS)
    )
}

# The per-annotation baseline: NO aggregation at all — every (patient, annotation)
# is its own point. It answers a different question from the grid (does the method
# agree annotation-by-annotation, rather than patient-by-patient) and is the honest
# reference the aggregated correlations should be read against, since it neither
# gains the smoothing nor loses the within-patient signal that aggregation does.
per_annotation_baseline <- function(per_ann, path_long,
                                    metrics = c("tumor_over_inside", "tumor_over_inside_clean")) {
  purrr::map_dfr(metrics, function(metric) {
    d <- per_ann |>
      dplyr::inner_join(path_long, by = c("patient_id", "annotation")) |>
      dplyr::filter(is.finite(.data[[metric]]), is.finite(path_frac))
    if (nrow(d) == 0) return(tibble::tibble())
    paired_cor3(d[[metric]], d$path_frac) |>
      dplyr::mutate(metric = metric, bias = mean(d[[metric]] - d$path_frac),
                    .before = 1)
  })
}

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

# Scatter grid: IHC aggregator (rows) x pathologist aggregator (columns), one figure
# per metric. Both axes are 0..1 fractions of the same quantity, so the dashed x = y
# line is the target and vertical distance from it is the disagreement. Matched-
# weighting panels are outlined so the interpretable cells stand out.
plot_aggregation_grid <- function(pairs, which_metric, title = NULL) {
  # `metric` is a COLUMN of `pairs`, so the argument is deliberately named
  # differently — `filter(pairs, metric == metric)` would be a tautology matching
  # every row and silently overlay all metrics in one panel.
  d <- dplyr::filter(pairs, metric == which_metric)
  if (nrow(d) == 0) return(invisible())

  rho <- d |>
    dplyr::group_by(ihc_lab, path_lab, matched = .is_matched(ihc_agg, path_agg)) |>
    dplyr::summarise(lab = sprintf("rho = %.2f (n = %d)",
                                   suppressWarnings(stats::cor(ihc_val, path_frac,
                                                               method = "spearman")),
                                   dplyr::n()),
                     .groups = "drop")

  ggplot(d, aes(ihc_val, path_frac)) +
    # Shading FIRST: a later geom_rect would paint over the x = y line and the points.
    geom_rect(data = dplyr::filter(rho, matched),
              aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              fill = "grey92", inherit.aes = FALSE) +
    geom_abline(slope = 1, linetype = "dashed", colour = "red") +
    geom_point(aes(colour = patient_id), alpha = 0.8, size = 2) +
    geom_text(data = rho, aes(x = -Inf, y = Inf, label = lab),
              hjust = -0.05, vjust = 1.3, size = 2.7, inherit.aes = FALSE) +
    facet_grid(ihc_lab ~ path_lab, labeller = label_wrap_gen(18)) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title %||% paste0("Aggregation sensitivity: ", metric),
         subtitle = paste("shaded panels weight both sides the same way;",
                          "dashed line is x = y"),
         x = "IHC tumour fraction (aggregated per patient)",
         y = "Pathologist tumour fraction (aggregated per patient)",
         colour = "Patient") +
    theme_classic(base_size = 10) +
    theme(strip.text.y = element_text(angle = 0))
}

# Heatmap of one correlation statistic over the grid, faceted by metric. Reading a
# ROW tells you how much the pathologist aggregator matters; a COLUMN, how much the
# IHC one does. A grid that is flat means the choice is immaterial for this cohort —
# which at n of a handful of patients is the outcome to hope for.
plot_aggregation_heatmap <- function(stats, stat = "spearman") {
  if (nrow(stats) == 0) return(invisible())
  d <- dplyr::mutate(stats, value = .data[[stat]])
  ggplot(d, aes(path_lab, ihc_lab, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    # Ring the matched-weighting cells so they read as the reference points.
    geom_tile(data = dplyr::filter(d, matched),
              colour = "black", linewidth = 0.9, fill = NA) +
    geom_text(aes(label = ifelse(is.finite(value), sprintf("%.2f\n(n=%d)", value, n), "—")),
              size = 3) +
    facet_wrap(~ metric, ncol = 1) +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                         midpoint = 0, limits = c(-1, 1), na.value = "grey90",
                         name = stat) +
    labs(title = paste0("Aggregation sensitivity of the IHC-vs-pathologist ", stat),
         subtitle = "boxed cells weight both sides the same way",
         x = "pathologist aggregator", y = "IHC aggregator") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank())
}

# How far apart the aggregators put a single patient: per patient, the spread
# (max - min) of its value across the aggregators on each side. A patient with one
# annotation has spread 0 by construction; large spreads identify the patients whose
# annotations disagree most, which is where the aggregation choice actually bites.
aggregation_spread <- function(pairs) {
  pairs |>
    dplyr::group_by(metric, patient_id) |>
    dplyr::summarise(
      n_pairs     = dplyr::n(),
      ihc_min     = min(ihc_val),   ihc_max   = max(ihc_val),
      ihc_spread  = max(ihc_val)   - min(ihc_val),
      path_min    = min(path_frac), path_max  = max(path_frac),
      path_spread = max(path_frac) - min(path_frac),
      .groups = "drop"
    ) |>
    dplyr::arrange(metric, dplyr::desc(ihc_spread))
}

