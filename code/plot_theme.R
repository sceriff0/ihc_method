# =============================================================================
# plot_theme.R  —  THE house figure style for this project.
#
# Every figure in analysis/*.Rmd and code/*.R goes through this file, so a plot
# rendered by clinical_data.Rmd and one rendered by benchmarks.Rmd are visually
# indistinguishable apart from their content. It is modelled on the journal
# figure the PI supplied as the target (Nature-style multi-panel):
#
#   * white panel, NO gridlines — the axes carry the reading, not a grid
#   * thin black axis lines with short OUTWARD ticks
#   * black text (not grey) in a small humanist sans (Helvetica/Arial)
#   * no strip background — facet labels are plain small bold text
#   * compact legends, colourbars horizontal and short
#   * colourblind-safe categorical palette (Okabe-Ito), blue-white-red diverging
#     and single-hue sequential continuous ramps
#
# HOW IT IS APPLIED. Sourcing this file has three side effects, on purpose:
#   1. theme_set(theme_paper())          -> every plot inherits the theme
#   2. paper_geom_defaults()             -> geom_point/line/text/boxplot defaults
#   3. options(ggplot2.discrete.*)       -> unspecified discrete scales use `oi`
# (2) and (3) are what make a NEW plot publication-ready without the author
# remembering anything. (1) only works if the plot does not append its own
# `theme_*()`: an inline `theme_classic()` REPLACES the active theme wholesale.
# So the rule for this repo is: never call `theme_classic()`/`theme_bw()`/
# `theme_minimal()` in an analysis. Add bare `theme(...)` for per-plot tweaks,
# which layers on top of the house theme instead of discarding it.
#
# SIZING FOR PRINT. The Rmds render at fig.width = 9in so the workflowr site is
# readable. A 9in-wide PDF placed in a 183mm (7.2in) double-column slot shrinks
# by 0.8x, so the 10pt base type lands at ~8pt — inside the 5-8pt journal range.
# For a single 89mm column, re-knit that chunk at fig.width = 4.5 rather than
# shrinking a wide figure, or the type drops below 5pt.
#
# Dependencies: ggplot2 + grid only (no tidyverse/here), so benchmark_plots.R,
# validation_helpers.R and any bare Rscript can all source it.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

# --- unit helpers ------------------------------------------------------------
# ggplot2 measures line widths and geom text in mm, journals specify them in pt.
# A `linewidth` of 1 draws at ~2.13pt; `size` in geom_text() is pt/2.845. These
# two converters let the specs below (and call sites) be written in real points.
pt_line <- function(pt) pt / 2.13          # pt -> ggplot2 `linewidth`
pt_text <- function(pt) pt / 2.845276      # pt -> geom_text/geom_label `size`

# --- font --------------------------------------------------------------------
# The style calls for a humanist sans (Helvetica/Arial), and "" is how you ask for
# it SAFELY. Do not name the font here.
#
# WHY. A family name has to be registered in the DEVICE's font database, not
# merely installed on the machine. grDevices::pdf() ships the Type 1 base-14 set
# (Helvetica, Times, Courier, ...) and knows nothing about an OS face like "Arial",
# "Nimbus Sans" or "DejaVu Sans" — those are what systemfonts::system_fonts()
# reports, which is a DIFFERENT database. Naming an OS-only face makes every knit
# with `dev = c("png", "pdf")` die at the first text grob with
#   Error in grid.Call.graphics(C_text, ...) : invalid font type
# and the PNG (cairo) pass gives no warning of it, because cairo resolves OS fonts
# happily. That failure mode cost a build; hence this comment.
#
# "" means "the device's own default family" — which IS Helvetica on pdf() and a
# Helvetica-metric sans on cairo/ragg. Identical look, no font database to satisfy,
# works unchanged on a headless cluster node.
#
# To force a specific face (a journal insisting on Arial, say): register it with
# the device first (see ?pdfFonts / ?grDevices::Type1Font), then set
#   options(ihc.plot.family = "Arial")
# BEFORE sourcing this file. paper_family() validates the request against the
# device databases and degrades to "" with a warning rather than letting a knit
# fail three chunks in.
paper_family <- function(family = getOption("ihc.plot.family", "")) {
  if (!nzchar(family)) return("")
  known <- tryCatch(unique(c(names(grDevices::pdfFonts()),
                             names(grDevices::postscriptFonts()))),
                    error = function(e) character(0))
  if (family %in% known) return(family)
  warning("plot_theme: font family ", sQuote(family), " is not registered with the ",
          "pdf/postscript device, so it would abort any knit that renders PDFs. ",
          "Falling back to the device default. Register it with grDevices::pdfFonts() ",
          "first if you need it.", call. = FALSE)
  ""
}

# --- the theme ---------------------------------------------------------------
# Built on theme_classic() because that is the only built-in with the right
# grammar (axis lines, no grid, no panel border); everything below re-specifies
# the parts theme_classic gets wrong for print — grey-ish text, chunky lines,
# inward ticks, a grey strip background, and oversized legend keys.
#
# Args let a caller deviate deliberately without hand-rolling a theme:
#   base_size   type size in pt (10 = the house default; see SIZING FOR PRINT)
#   grid        "none" (default), "y", "x", or "both" — a faint reference grid,
#               for the rare panel (e.g. a wide dot plot) that is unreadable
#               without one
#   axis_lines  FALSE drops the L-shaped axis, for heatmaps/tile plots
theme_paper <- function(base_size = 10, base_family = paper_family(),
                        grid = c("none", "y", "x", "both"), axis_lines = TRUE) {
  grid <- match.arg(grid)
  line_col <- "black"
  grid_line <- element_line(linewidth = pt_line(0.25), colour = "grey92")

  th <- theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # -- text: black, tight, left-aligned titles over the WHOLE plot so a
      # title and its y-axis label do not fight for the same column.
      text             = element_text(colour = "black"),
      plot.title       = element_text(face = "bold", size = rel(1.0),
                                      margin = margin(b = 4)),
      plot.subtitle    = element_text(colour = "grey30", size = rel(0.85),
                                      margin = margin(b = 6)),
      plot.caption     = element_text(colour = "grey45", size = rel(0.65),
                                      hjust = 1, margin = margin(t = 6)),
      plot.title.position = "plot", plot.caption.position = "plot",
      # Panel letters (a, b, c...) for patchwork::plot_annotation(tag_levels = "a").
      # NOTE: a tag and a plot.title share the top-left slot and WILL overlap. That
      # is the correct trade-off, because a panel inside an assembled figure should
      # not carry its own title anyway — the caption does that work. When composing
      # the final figure, set labs(title = NULL, subtitle = NULL) on each panel and
      # let patchwork place the tags.
      plot.tag          = element_text(face = "bold", size = rel(1.3), hjust = 0),
      plot.tag.position = "topleft",

      # -- axes: hairline black rules, short OUTWARD ticks, black labels
      axis.line   = if (axis_lines) element_line(linewidth = pt_line(0.75),
                                                 colour = line_col,
                                                 lineend = "square")
                    else element_blank(),
      axis.ticks  = if (axis_lines) element_line(linewidth = pt_line(0.75),
                                                 colour = line_col)
                    else element_blank(),
      axis.ticks.length = unit(2, "pt"),
      axis.text   = element_text(colour = "black", size = rel(0.85)),
      axis.title  = element_text(colour = "black", size = rel(0.95)),

      # -- panel/grid
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = if (grid %in% c("y", "both")) grid_line else element_blank(),
      panel.grid.major.x = if (grid %in% c("x", "both")) grid_line else element_blank(),
      # 9pt, not less: with no panel border the only thing keeping one panel's
      # last x tick label off its neighbour's first one is this gap.
      panel.spacing    = unit(9, "pt"),

      # -- facets: no grey box, just small bold text sitting on the panel
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = rel(0.85),
                                      colour = "black",
                                      margin = margin(3, 3, 3, 3)),

      # -- legends: compact, top-left, no box. Colourbars are short and thin;
      # `legend.key.*` is also what sizes guide_colourbar() in ggplot2 >= 3.5.
      legend.position   = "top",
      legend.justification = "left",
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.title      = element_text(size = rel(0.85)),
      legend.text       = element_text(size = rel(0.8)),
      legend.key.height = unit(8, "pt"),
      legend.key.width  = unit(14, "pt"),
      legend.margin     = margin(0, 0, 2, 0),

      plot.margin = margin(6, 8, 4, 6)
    )
  th
}

# Heatmap / tile variant: cells ARE the panel, so the L-shaped axis and the ticks
# only add clutter. Keeps the type, legend and title styling identical.
theme_paper_tile <- function(base_size = 10, ...) {
  theme_paper(base_size = base_size, axis_lines = FALSE, ...) +
    theme(axis.ticks.length = unit(0, "pt"))
}

# --- palettes ----------------------------------------------------------------
# Okabe-Ito: the standard 8-colour colourblind-safe categorical palette. Ordered
# so the first two (blue, vermillion) are the maximally separable pair, which is
# what a 2-level scale gets.
oi <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
        "#E69F00", "#56B4E9", "#F0E442", "#000000")

# `oi` runs out at 8, and several figures here colour by `patient_id` with more
# patients than that. Extended set: the first 8 entries ARE `oi`, so a 3-patient
# and a 14-patient panel still open with the same blue; the tail adds darker,
# still-distinguishable hues from the Okabe-Ito barrier-free set. Past ~12
# categories colour stops separating anything — facet, or label points directly,
# rather than reaching for more colours.
oi_ext <- c(oi,
            "#004949", "#920000", "#490092", "#B66DFF",
            "#6DB6FF", "#924900", "#009292", "#FF6DB6")

# Diverging (ColorBrewer RdBu endpoints): for signed quantities around a real
# zero — log2 fold changes, correlations, differences. Blue = low/negative.
PAL_DIV <- c(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B")

# Sequential single-hue: for magnitudes with no meaningful midpoint (counts,
# densities, fractions). Perceptually ordered and safe in greyscale.
PAL_SEQ <- c(low = "#F7FBFF", high = "#08519C")

# Semantic colours for the two recurring ANNOTATION marks, so "the dashed x = y
# line" and "the fitted trend" look the same in every report. Reference lines are
# deliberately achromatic — a red identity line competes with the data for
# attention, and in these plots the data is what carries the finding.
REF_LINE <- "grey35"       # identity line, threshold, bias / limits of agreement
FIT_LINE <- "#0072B2"      # geom_smooth trend (= oi[1])

# Semantic colours for the immune "hot/cold" phenotype: HOT -> red, COLD -> light
# blue, intermediate/other -> orange/grey. Robust to case/spelling. Returns a
# named vector keyed by the given levels, for scale_colour_manual/scale_fill_manual.
hotcold_cols <- function(levels) {
  lv  <- as.character(levels)
  key <- toupper(trimws(lv))
  col <- ifelse(grepl("HOT|INFLAM", key),               "#D7191C",   # red
         ifelse(grepl("COLD|DESERT", key),              "#74ADD1",   # light blue
         ifelse(grepl("INTERMED|VARI|MIX|EXCLUD", key), "#FDAE61",   # orange
                                                        "grey65")))
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

# --- scale shorthands --------------------------------------------------------
# Thin wrappers so a call site names the SEMANTICS ("this is diverging") rather
# than repeating hex codes. Every `...` passes through to the underlying scale,
# so limits/name/labels/trans all still work.
scale_fill_div <- function(midpoint = 0, ...)
  scale_fill_gradient2(low = PAL_DIV[["low"]], mid = PAL_DIV[["mid"]],
                       high = PAL_DIV[["high"]], midpoint = midpoint,
                       na.value = "grey92", ...)

scale_colour_div <- function(midpoint = 0, ...)
  scale_colour_gradient2(low = PAL_DIV[["low"]], mid = PAL_DIV[["mid"]],
                         high = PAL_DIV[["high"]], midpoint = midpoint,
                         na.value = "grey92", ...)

scale_fill_seq <- function(...)
  scale_fill_gradient(low = PAL_SEQ[["low"]], high = PAL_SEQ[["high"]],
                      na.value = "grey92", ...)

scale_colour_seq <- function(...)
  scale_colour_gradient(low = PAL_SEQ[["low"]], high = PAL_SEQ[["high"]],
                        na.value = "grey92", ...)

scale_colour_oi <- function(...) scale_colour_manual(values = unname(oi), ...)
scale_fill_oi   <- function(...) scale_fill_manual(values = unname(oi), ...)

# ORDINAL discrete: levels that have an order (image size, tile count) or simply
# outnumber the 8 Okabe-Ito colours. `oi` is wrong for both — an unordered hue set
# hides the ordering, and scale_*_manual errors outright when it runs out of
# colours. One ramp (viridis D) for every such case, so ordered discrete looks the
# same everywhere; the repo previously mixed viridis options B, C and D.
scale_colour_ordinal <- function(...) scale_colour_viridis_d(option = "D", ...)
scale_fill_ordinal   <- function(...) scale_fill_viridis_d(option = "D", ...)
scale_color_ordinal  <- scale_colour_ordinal

# American/British spelling aliases, so a call site can use either.
scale_fill_diverging <- scale_fill_div
scale_color_div      <- scale_colour_div
scale_color_seq      <- scale_colour_seq
scale_color_oi       <- scale_colour_oi

# Short horizontal colourbar with the title above it, matching the reference
# figure's under-panel bars. ggplot2 3.5 moved bar sizing out of the guide and
# into the theme (and deprecated barwidth/barheight), so detect which API this
# installation has instead of pinning one.
# `title` defaults to waiver(), NOT NULL: in ggplot2 guides waiver() means
# "inherit the scale's name" while NULL means "draw no title at all", so a NULL
# default would silently strip the label off every colourbar it is used on.
guide_cbar <- function(title = waiver(), width = 72, height = 6, ...) {
  args <- list(title = title, title.position = "top", direction = "horizontal",
               ticks = FALSE, ...)
  if ("theme" %in% names(formals(guide_colourbar))) {
    args$theme <- theme(legend.key.width  = unit(width, "pt"),
                        legend.key.height = unit(height, "pt"))
    # 3.5+ renamed these to `theme` entries; drop the deprecated spellings.
    args$title.position <- NULL
    args$theme <- args$theme + theme(legend.title.position = "top")
  } else {
    args$barwidth  <- unit(width, "pt")
    args$barheight <- unit(height, "pt")
  }
  do.call(guide_colourbar, args)
}

# --- geom defaults -----------------------------------------------------------
# The theme cannot reach inside a geom, so a default-sized geom_point (size 1.5,
# ~4pt on paper) and 0.5-linewidth lines stay chunky no matter how good the theme
# is. These bring the marks themselves down to print scale. Explicit sizes at a
# call site still win, so this only moves the plots nobody tuned by hand.
paper_geom_defaults <- function() {
  set <- function(geom, vals)
    try(update_geom_defaults(geom, vals), silent = TRUE)
  set("point",   list(size = 1.2, stroke = 0.3))
  set("line",    list(linewidth = pt_line(0.75)))
  set("path",    list(linewidth = pt_line(0.75)))
  set("segment", list(linewidth = pt_line(0.75)))
  set("hline",   list(linewidth = pt_line(0.5), colour = "grey40"))
  set("vline",   list(linewidth = pt_line(0.5), colour = "grey40"))
  set("abline",  list(linewidth = pt_line(0.5)))
  set("smooth",  list(linewidth = pt_line(1.0)))
  # Outline-only marks need more stroke than a data line: at 0.5pt a pale fill
  # colour (the light blue in the hot/cold scale, say) washes out to grey while
  # the thicker median segment inside the same box still reads as its true hue.
  set("boxplot", list(linewidth = pt_line(0.9)))
  set("bar",     list(linewidth = pt_line(0.9)))
  set("col",     list(linewidth = pt_line(0.9)))
  set("errorbar",list(linewidth = pt_line(0.75)))
  # In-panel annotation text at ~7pt. Deliberately does NOT set `family`: these
  # defaults feed annotate("text", ...) and geom_text() alike, and a family the pdf
  # device cannot resolve aborts the knit there (see the font section above).
  # Leaving it unset inherits "" = the device default, same as the theme.
  set("text",    list(size = pt_text(7), colour = "black"))
  set("label",   list(size = pt_text(7), colour = "black"))
  invisible(NULL)
}

# --- apply -------------------------------------------------------------------
# Side effects on source(): this is what makes the style automatic rather than
# something each Rmd has to remember. Idempotent, so sourcing twice is harmless.
theme_set(theme_paper())
paper_geom_defaults()

# Any discrete scale the author did NOT specify falls back to these instead of
# ggplot2's evenly-spaced hue rainbow (which is neither colourblind-safe nor
# printable in greyscale). Passed as a LIST of palettes: ggplot2 takes the first
# one long enough for the variable's levels, so <=8 categories get `oi` and the
# 9..16 range gets `oi_ext` rather than silently reverting to hue. (The list
# entries must be character vectors — ggplot2 rejects functions here.) Beyond 16
# levels ggplot2 does fall back to its own scale; that is a signal the plot needs
# faceting or direct labels, not a longer palette.
options(ggplot2.discrete.colour = list(unname(oi), unname(oi_ext)),
        ggplot2.discrete.fill   = list(unname(oi), unname(oi_ext)),
        ggplot2.continuous.colour = function(...) scale_colour_seq(...),
        ggplot2.continuous.fill   = function(...) scale_fill_seq(...))
