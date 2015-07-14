# Version 2.0.0

Version 2.0.0 now uses ggplot2 to generate plots in `plot_event` and `plot_expr`. In addition to outputting a plot, a ggplot2 object is returned.

- Printing of the plot can be suppressed with `plot = FALSE`.
- The option `lines` in `plot_event` and `plot_expr` is deprecated (will issue warning if used).
- Imports plyr for calculating group means.

### NEW
- `plot_multi`: a new experimental function that will generate a heatmap of PSI values for two or more events:
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

- `plot_expr`: psiplot can now plot cRPKMs. For example, using the new sample dataset `crpkm`:
```r
plot_expr(crpkm[1,], config = config)
``` 
- Renamed sample datasets to: `psi` and `config` (#1)
- Deprecated `xlim` argument in `plot_event` - it shouldn't really be used

# Version 1.0.0

- This is the first release of psiplot in the form of an R package.
