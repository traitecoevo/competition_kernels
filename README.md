# Competition kernels

Ideas and code for the competition kernels.

Requires:
  - remake
  - callr
  - plant (>= 0.2.2)

## Compilation instructions

First, install remake:

```r
install.packages("drat")
drat:::add("traitecoevo")
install.packages("remake")
```

Then, run install any missing packages (probably the three above, plus dependencies and a few extras: see the [remake file](https://github.com/traitecoevo/competition_kernels/blob/master/remake.yml)) and run the actual analysis.

```r
remake::install_missing_packages()
remake::make()
```

Compilation will take a few hours and will use up to two cores by default.  Using `options(mc.cores=4L)` will use more cores on the slowest bits.

Generating the PDF documents (manuscript, supplementary information and one vignette) will require a reasonably complete LaTeX installation and [pandoc](http://pandoc.org/installing.html).

Installing `plant` is not currently possible on Windows (see the [plant repo](https://github.com/traitecoevo/plant)) but we are hopeful that this will be resolved shortly.
