# Code

Save command-line scripts and shared R code here.

## Figure style

`plot_theme.R` is the single house style for every figure on the site. It is
sourced by `validation_helpers.R` and by `benchmark_plots.R`, so any analysis
that loads either one is styled automatically — there is nothing to call.

It exports:

| | |
|---|---|
| `theme_paper(base_size, grid, axis_lines)` | the theme; applied via `theme_set()` on source |
| `theme_paper_tile()` | heatmap variant (no axis line, no ticks) |
| `oi`, `oi_ext` | Okabe-Ito categorical palette (8) and its 16-colour extension |
| `scale_*_oi()`, `scale_*_ordinal()` | categorical / ordered-discrete scales |
| `scale_*_div()`, `scale_*_seq()` | diverging (blue-white-red) and sequential ramps |
| `guide_cbar()` | short horizontal colourbar |
| `hotcold_cols()`, `hotcold_order()` | immune-phenotype colours and level order |
| `REF_LINE`, `FIT_LINE` | colours for reference lines and fitted trends |
| `pt_line()`, `pt_text()` | pt -> ggplot2 `linewidth` / geom text `size` |

**The one rule:** never write `theme_classic()` / `theme_bw()` / `theme_minimal()`
in an analysis. Those *replace* the active theme, which is how the figures drifted
apart in the first place. A bare `theme(...)` layers on top and is fine — use it
for genuine per-plot deviations (rotated tick labels, horizontal facet strips).

Sourcing `plot_theme.R` also sets ggplot2's `discrete.colour`/`discrete.fill`
options, so a scale you *don't* specify falls back to the house palette instead of
ggplot2's hue rainbow. New plots are publication-ready by default.

**Do not name a font family.** `theme_paper()` uses `""`, which means "the device's
own default" — Helvetica on `pdf()`, a Helvetica-metric sans on cairo/ragg. A named
family has to be registered in the *device's* database (`pdfFonts()`), not merely
installed on the machine, so an OS face like `Arial` or `DejaVu Sans` aborts every
knit that renders PDFs with `invalid font type` — and the PNG pass gives no warning,
because cairo resolves OS fonts happily. If a journal demands a specific face,
register it with `grDevices::pdfFonts()` and set `options(ihc.plot.family = "...")`
before sourcing; `paper_family()` validates it and degrades with a warning.

For a final multi-panel figure, drop the per-panel `title`/`subtitle`
(`labs(title = NULL, subtitle = NULL)`) and let `patchwork::plot_annotation(tag_levels = "a")`
place the panel letters — a tag and a title share the same slot and will overlap.
