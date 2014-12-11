<!-- README.md is generated from README.Rmd. Please edit that file -->



psiplot
=======

psiplot is an R package for generating plots of percent spliced-in (PSI) values of alternatively-spliced exons. It is based on the `plot` tool provided by [vast-tools](https://github.com/vastgroup/vast-tools), an RNA-Seq pipeline for alternative splicing analysis.

Currently, the plot code for vast-tools and psiplot are maintained separately, even though they are technically identical. Future work is planned to integrate psiplot into the vast-tools `plot` code.

Installation
------------

The most up-to-date development version can be obtained via devtools:

``` {.r}
install.packages("devtools")
devtools::install_github("kcha/psiplot")
```

Usage
-----

psiplot takes as input the PSI results generated by vast-tools (e.g. after running `vast-tools combine` or `vast-tools diff`).

``` {.r}
library(psiplot)

# Plot an event using provided sample dataset
plot_event(psi[1,])
```

### Customizing plots

#### The `.config` file way

In `vast-tools plot`, an optional config file can be used to customize the plots' visual appearance. The same config file can be supplied here as well.

``` {.r}
plot_event(psi[1,], config="/path/to/config")

# config can also be pre-loaded into a data frame
cfg <- read.table("/path/to/config", header=T, sep="\t", stringsAsFactor=FALSE)
plot_event(psi[1,], config=cfg)
```

The color and ordering of samples can be customized by supplying a plot configuration file. This file is tab-delimited and must be manually created. For example:

``` {.r}
config
#>   Order SampleName GroupName RColorCode
#> 1     1    Sample4    Neural    #ff0000
#> 2     2    Sample3    Neural        red
#> 3     3    Sample2    Muscle       blue
#> 4     4    Sample1    Muscle    #0000ff
```

etc..

-   **Order**: The ordering of the samples from left to right.
-   **SampleName**: Name of the sample. MUST match sample name in input table.
-   **GroupName**: Group name. Use for plotting the average PSI of samples belonging to the same group (enable by setting `groupmean=TRUE`)
-   **RColorCode**: An R color specification:
    1.  color name (as specified by `colors()`)
    2.  hex color code (\#rrggbb)

The samples under SampleName MUST MATCH the names in the PSI input table. Only the samples listed in the config file will be represented in the resulting plots. Other samples in the PSI table but not in the config file will be ignored. This may be useful if you want to customize the type of samples in your plots.

#### The R way

The colors and other graphical paramters can also be configured in R. `plot_events()` provides some limited options. See `?plot_events` for more details on the available options.

``` {.r}
plot_event(psi[1,], config = config, pch = 9, ylim = c(20, 80))
```

Related Projects
----------------

-   [vast-tools](https://github.com/vastgroup/vast-tools)
