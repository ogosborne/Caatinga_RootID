# Caatinga_RootID

This repository contains code used to match known trees to unknown mixed root samples in a plot in the Caatinga region of Brazil. The analysis starts with ddRADseq data from leaves of known trees and from mixed roots sampled from across the plot. The data is first processed by the *STACKS* pipeline, and stacks outputs are analysed using the *RootID* package in *R*.

To reproduce the analysis, sequence data should first be downloaded from the European Nucleotide Archive (n.b. accession to be added on manuscript publication), and placed in a directory: *seq_data*. The shell scripts were run in a linux environment (Red Hat Enterprise Linux 8.5 running a bash shell) and the R scripts were run in R version 4.1.2 on a MacBook Pro (MacOS 12.3). To replicate the analysis, [STACKS](https://catchenlab.life.illinois.edu/stacks/) (version used: 2.52) and R should be installed and available in the path, and the following R packages should be installed:

* [RootID](https://github.com/ogosborne/RootID) (version used: 1.0.0)
* [UpSetR](https://cran.r-project.org/package=UpSetR) (version used: 1.4.0)
* [proxy](https://cran.r-project.org/package=proxy) (version used: 0.4-27)
* [ape](https://cran.r-project.org/package=ape) (version used: 5.6-2)
* [ppcor](https://cran.R-project.org/package=ppcor) (version used: 1.1)
* [viridis](https://cran.r-project.org/package=viridis) (version used: 0.6.2)
* [coin](https://cran.r-project.org/package=coin) (version used: 0.4-2)
* [gtools](https://cran.r-project.org/package=gtools) (version used: 3.9.2)
* [RColorBrewer](https://cran.r-project.org/package=RColorBrewer) (version used: 1.1-2)
* [rgl](https://cran.r-project.org/package=rgl) (version used: 0.109.6)

The scripts can then be run in order as follows:

#### Prepare sequence data.
```bash
# Run process_radtags
bash ./shell/process_radtags_4RootID.sh
# Concatenate reads.
bash ./shell/cat_reads.sh
```

#### Run stacks with multiple parameter settings.
```bash
# run ustacks
bash ./shell/ustacks_4RootID.sh
# run cstacks
bash ./shell/cstacks_4RootID.sh
# run sstacks
bash ./shell/sstacks_4RootID.sh
```

#### Run RootID across parameter combinations, compare and choose optimal parameters.
```bash
# Run RootID for all parameter combinations. Record timing for each run.
Rscript ./R/01_RootID_all_params.R
# Choose optimal parameter combinations.
Rscript ./R/02_Find_optimal_params.R
# Compare matches across parameter combinations.
Rscript ./R/03_Compare_matches_across_params.R
# Rerun match.diag with read depth filter
Rscript ./R/04_Match.diag_depth_filt.R
```

#### Run pipeline validation and assessment analyses
```bash
# Run rarefaction analysis.
Rscript ./R/05_Rarefaction_analysis.R
# Various descriptive statistics and correlations re: pipeline performance etc
Rscript ./R/06_Miscellaneous_stats.R
# Comparison of tree and root position
Rscript ./R/07_tree-root_position_comparison.R
# Make dendrogram based on shared markers
Rscript ./R/08_pres-abs_dendrogram.R
# compare tree and root positions??
```

#### Analyse root distribution patterns
```bash
# Aboveground-belowground correlations
Rscript ./R/09_above-belowground_cor.R
# depth niche analysis
Rscript ./R/10_depth_niches.R
```

#### Visualise 3D root distributions.

The final script ```11_3d_plots.R``` should be run line-by-line in the R GUI or Rstudio so the plots can be manually adjusted (turned/zoomed etc) in the RGL device.
