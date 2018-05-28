# Version 2.2.0

- Support for sample subgroups in all plot functions, through the config files. PSIs of samples in a subgroup are averaged and displayed as a single data point. Confidence intervals are computed by fitting a joint beta distribution to a population of points sampled from beta distributions coming from the individual samples, and then applying the same method as for single samples.
- Support for plotting multiple events in the same PSI plot, with events identified with different colours, and samples identified by the color background.
- Support for sample columns containing suffixes (e.g. "Sample1-cRPKM", while the config says only "Sample1").
- Support for gene expression tables containing both cRPKM and raw read counts (similar to the PSI and QC columns for events).
- Big changes in the preprocess_sample_colors function and its output, to allow for all this.
- Use of dplyr (and tidyverse in general) instead of plyr.
- README contents moved to vignette

# Version 2.1.2

- Fixed bug that occurs when dealing with columns that are all NA
- Added CITATION
- Converted some internally-used functions to be non-exported
- `plot_multi()` now supports `pheatmap`
- The `test_data` folder is now `inst/extdata`

# Version 2.1.1

Some bug fixes for `plot_expr()`

# Version 2.1.0

This release contains new features and some internal updates. In addition, this release now requires **R version 3.1** or higher and **ggplot2 version 2.0** or higher.

- Replace `show_guides` with `show.legend` to match API of ggplot2 version 2.*.  
- Update README
- `plot_multi` has been further developed: 
  - Use `heatmap.2` from gplots package to produce heatmap. To minimize dependencies, this package is optional and needs to be manually installed (`install.packages("gplots")`). In the absence of gplots, ggplot2 will be used as before.
  - Events and samples can be clustered using options `cluster_rows` and `cluster_cols`, respectively.
  - Add `fill` option for providing custom colours
  - Simply format of row names
  - See `?plot_multi` for full set of options
- Retire warning for deprecated `lines` option
- The option `plot` in `plot_event` and `plot_expr` is deprecated (will issue warning if used). 
- Other documentation updates


# Version 2.0.1

This release contains minor updates.

- `plot_event`: default `ylim` changed to (0,100) from (1,100).

# Version 2.0.0

Version 2.0.0 now uses ggplot2 to generate plots in `plot_event` and
`plot_expr`. In addition to outputting a plot, a ggplot2 object is returned.

- Printing of the plot can be suppressed with `plot = FALSE`.
- The option `lines` in `plot_event` and `plot_expr` is deprecated (will issue
  warning if used).
- Imports plyr for calculating group means.

### NEW
- `plot_multi`: a new experimental function that will generate a heatmap of PSI
  values for two or more events:
```r
plot_multi(psi)
```

# Version 1.1.2

- Change y-axis labels to be horizontal (e.g. `las = 1`)

# Version 1.1.1

- Changed minimum of y-axis to start at 0 for `plot_expr()` (#2)
- Added check if config and input data does has 0 matching samples (#3)
- Other minor updates to documentation

# Version 1.1.0

- `plot_expr`: psiplot can now plot cRPKMs. For example, using the new sample
  dataset `crpkm`:
```r
plot_expr(crpkm[1,], config = config)
``` 
- Renamed sample datasets to: `psi` and `config` (#1)
- Deprecated `xlim` argument in `plot_event` - it shouldn't really be used

# Version 1.0.0

- This is the first release of psiplot in the form of an R package.
