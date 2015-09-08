# Response of phylogenetic diversity to human impacts in terrestrial ecosystems

All code is made available for reproducing the analyses as described in
[insert paper URL here]. This is code only, for full completed dataset click
here (not yet available).

## Pipeline

The analysis comes in four stages:

![r_pipeline](https://raw.githubusercontent.com/DomBennett/PREDICTS-PD/master/other/r_pipeline.jpg?token=AC0Jh5LHfPQDyU_Vii3xzDcNnijLRdUTks5V-ByPwA%3D%3D "R analysis pipeline")

Each stage is run once with the exception of `parse.R`, `compare.R` and `metrics.R`
which are each run twice for comparing the results from unconstrained and constrained
pG-lt trees and testing whether there is a difference between total and age
branch length normalisation.

## Running

Make sure all relevant packages are installed:

```{bash}
Rscript install_deps.R
```

Additionally, you will need a compiled version of [pathD8](http://www2.math.su.se/PATHd8/)
in your working directory.

To run entire automated pipeline (UNIX OS):

```{bash}
sh run.sh
```

pG-lt step can only be run separately. It was originally run on Imperial's HPC.
`analysis.R` is an interactive script and it not part of the pipeline.
Note, setup will take a long time to run as they interact
with the [GNR](http://resolver.globalnames.biodinfo.org/).

## Stage details

### Setup

This is run with the `setup.sh` script.

* `setup_presolve.R`: read in published trees and PREDICTS data, search names
against [GNR](http://resolver.globalnames.biodinfo.org/).
* `setup_parse.R`: read in published trees and rate-smooth
* `map.R`: use published phylogenies in `0_data` to generate study-level
phylogenies
* `pGltsetup.R`: identify suitable PREDICTS studies, output into pG-lt friendly
format

### Calculate

This is run with the `calculate.R` script.

* `pGltrun.R`: placeholder script -- doesn't do anything
* `parse.R`: read in pG-lt and mapped phylogenies, check and rate-smooth
* `compare.R`: compare mapped and pG-lt phylogenies
* `metrics.R`: generate PD, PSV and PSE values per site using both pG-lt and
mapped trees.
* `commplots.R`: community plot PREDICTS data onto pG-lt consensus phylogenies

## Directory

```
-- 0_data/
    -- PREDICTS-DATA/
        -- [site-level PREDICTS data]
    -- raw_trees/
        -- [downloaded published trees]
    -- parsed_trees/
        -- [rate-smoothed trees]
-- 1_pGltsetup/
    -- [folders setup for pG-lt]
-- 2_pGltrun/
    -- [pG-lt run results]
-- 3_map/
    -- [study-level trees mapped from published trees]
-- 4_parse/
    -- [rate-smoothed pG-lt trees]
-- 5_compare/
    -- [results from comparing mapped and pG-lt trees]
-- 6_metrics/
    -- [site-level PREDICTS data with phylogenetic metrics]
-- 7_commplots/
    -- [plots of trees with species incidences and abundances per site]
-- 8_analysis/
    -- [results from human disturbance ~ phylogenetic metrics]
-- stages/
    -- [all R stages files]
-- tools/
    -- [all custom R functions]
-- sanity_checks/
    -- [demonstrations of code accuracy]
-- other/
    -- [misc non-essential items]
```

## Data

After cloning this repo to your machine, you can download all missing data files
with this link (not yet available).

Uncompress files:

```{bash}
sh uncompress.sh
```

And compress again:

```{bash}
sh compress.sh
```

## Authors
[Bennett D.J.](https://github.com/DomBennett),
[De Palma A.](https://github.com/adrianadepalma),
[Pearse W.D.](https://github.com/willpearse)
and [Purvis A](https://github.com/AndyPurvis).
