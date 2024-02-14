# IDS

Code and data for the case study in:

Kéry M, Royle JA, Hallman T, Robinson WD, Strebel N, Kellner KF (2024). Integrated distance sampling models for simple point counts. Ecology. 

The version of `unmarked` with the `IDS()` function used in the analyses below should be on CRAN sometime in mid-2024, but in the meantime a copy can be found in this repository, or [here](https://github.com/rbchan/unmarked).

Until it's available on CRAN, this version of `unmarked` can be installed from the source package.
This requires [RTools](https://cran.r-project.org/bin/windows/Rtools/) on Windows, or [Xcode](https://developer.apple.com/xcode/) on Mac.

```r
depends <- tools::package_dependencies('unmarked')[[1]]
depends <- depends[! depends %in% c("graphics", "methods", "parallel", "stats", "utils")]
install.packages(depends)
install.packages("unmarked_1.4.1.9001.tar.gz", repos = NULL, type="source")
```

or

```r
remotes::install_github("rbchan/unmarked")
```

**Contents**

`README.md`: This file

`1_Code_Simulation_1.R`: Code for Simulation 1 (Identifiability of separate observation and process parameters in IDS1 and IDS2)

`2_Code_Simulation_1B.R`: Code for Simulation 1B (Parameter identifiability and estimator performance with IDS1)

`3_Code_Simulation_2.R`: Code for Simulation 2 (Identifiability with distinct covariate effects in the observation model)

`4_Code_Simulation_3.R`: Code for Simulation 3 (How many DS sites are required to obtain adequate estimates of density?)

`5_Code_Simulation_5.R`: Code for Simulation 4 (How well can availability be estimated in an IDS model?)

`IDS1_CaseStudy.R`: Code for American Robin case study (using JAGS and unmarked)

`IDS1_CaseStudy_model_5.txt`: JAGS model for American Robin case study

`data`: Folder containing bird sampling data and spatial data used in the case study (see `IDS1_CaseStudy.R` for descriptions of these files)

`unmarked_1.4.1.9001.tar.gz`: A copy of the R package `unmarked` that contains the `IDS` function used in some simulations and the case study
