# =============================================================================
# Standalone per-plot PDF export for the workflowr site.
#
# workflowr renders each chunk's plots to PNGs under the site figure path
# (docs/figure/<Rmd>/) for the website. When a chunk also sets
# `dev = c("png", "pdf")`, knitr emits a vector PDF of EVERY plot as well — one
# file per plot, so loops and multi-plot chunks each yield their own PDF, and
# base-R plots are captured just like ggplots.
#
# `export_pdf_figures(slug)` is meant to run in a final `include=FALSE` chunk of
# an analysis. It collects those generated PDFs and copies each into
# output/figures/<slug>/ so every figure is available as a single publication
# PDF, without disturbing the PNGs the HTML points at.
#
# It is deliberately fail-safe: any error is caught and downgraded to a message
# so a figure-export hiccup can never abort a knit. If no PDFs are found (e.g.
# `dev` was not set to include "pdf") it says so and does nothing.
#
# Dependencies: base R + knitr + here only (no tidyverse/sf), so it can be
# sourced from any analysis, including ones that do not load validation_helpers.
# =============================================================================

export_pdf_figures <- function(slug, out_root = here::here("output", "figures")) {
  tryCatch({
    fp <- knitr::opts_chunk$get("fig.path")   # e.g. "figure/clinical_data.Rmd/"

    # Candidate directories that may hold this knit's PDF figures. workflowr's
    # runtime fig.path is relative and the exact on-disk location varies by build
    # stage, so search the obvious roots and keep whatever exists.
    fig_roots <- character(0)
    if (!is.null(fp) && nzchar(fp)) {
      fig_roots <- c(fig_roots, fp, file.path(getwd(), fp), here::here(fp))
    }
    for (base in c(here::here("docs", "figure"),
                   here::here("analysis", "figure"),
                   here::here("figure"))) {
      if (dir.exists(base)) {
        subs <- list.dirs(base, recursive = FALSE)
        # keep only sub-directories whose name refers to this analysis
        fig_roots <- c(fig_roots, subs[grepl(slug, basename(subs), fixed = TRUE)])
      }
    }
    fig_roots <- unique(fig_roots[dir.exists(fig_roots)])

    pdfs <- unlist(lapply(fig_roots, function(d)
      list.files(d, pattern = "\\.pdf$", full.names = TRUE)), use.names = FALSE)
    pdfs <- unique(pdfs)

    if (!length(pdfs)) {
      message('export_pdf_figures("', slug, '"): no PDF figures found — ',
              'is `dev = c("png", "pdf")` set in this Rmd\'s setup chunk?')
      return(invisible(character(0)))
    }

    dest <- file.path(out_root, slug)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(pdfs, file.path(dest, basename(pdfs)), overwrite = TRUE)
    message(sprintf('export_pdf_figures("%s"): copied %d/%d PDF(s) to %s',
                    slug, sum(ok), length(pdfs), dest))
    invisible(file.path(dest, basename(pdfs))[ok])
  }, error = function(e) {
    message('export_pdf_figures("', slug, '") skipped: ', conditionMessage(e))
    invisible(character(0))
  })
}
