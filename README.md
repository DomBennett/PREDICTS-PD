# Response of phylogenetic diversity to human impacts in terrestrial ecosystems

All code is made available for reproducing the analyses as described in
[insert paper URL here]. This is code only, for full completed dataset click
here (not yet available).

## Pipeline

The analysis comes in four stages:

![r_pipeline](https://raw.githubusercontent.com/DomBennett/PREDICTS-PD/master/other/r_pipeline.jpg?token=AC0Jh5LHfPQDyU_Vii3xzDcNnijLRdUTks5V-ByPwA%3D%3D "R analysis pipeline")

Each stage is run once with the exception of `B1_parse.R`, `B2_compare.R` and `B3_metrics.R`
which are each run twice for comparing the results from unconstrained and constrained
pG-lt trees and testing whether there is a difference between total and age
branch length normalisation.

## Running

Make sure all relevant packages are installed:

```{bash}
Rscript x_install_deps.R
```

Additionally, you will need a compiled version of [pathD8](http://www2.math.su.se/PATHd8/)
in your working directory.

To run entire automated pipeline (UNIX OS):

```{bash}
sh x_run.sh
```

pG-lt step can only be run separately. It was originally run on Imperial's HPC.
Make sure you have `0_pglt/` before running the calculate step.
`C1_analysis.R` is an interactive script and it not part of the pipeline.
Note, setup will take a long time to run as they interact
with the [GNR](http://resolver.globalnames.biodinfo.org/).

## Stage details

### Setup

This is run with the `x_setup.sh` script. You can run pG-lt after it is run.

* `A1_presolve.R`: read in published trees and PREDICTS data, search names
against [GNR](http://resolver.globalnames.biodinfo.org/).
* `A2_parse.R`: read in published trees and rate-smooth
* `A3_pgltsetup.R`: identify suitable PREDICTS studies, output into pG-lt friendly
format. You can run pG-lt with the resulting folder.
* `A4_map.R`: use published phylogenies in `0_data` to generate study-level
phylogenies

### Calculate

This is run with the `x_calculate.R` script.

* `B1_parse.R`: read in pG-lt and mapped phylogenies, check and rate-smooth
* `B2_compare.R`: compare mapped and pG-lt phylogenies
* `B3_metrics.R`: generate PD, PSV and PSE values per site using both pG-lt and
mapped trees.
* `B4_plots.R`: community plot PREDICTS data onto pG-lt consensus phylogenies

## Directory (after complete run)

```
-- 0_data/
    -- PREDICTS-DATA/
        -- [site-level PREDICTS data]
    -- raw_trees/
        -- [downloaded published trees]
-- 0_pglt/
-- A1_preresolve/
-- A2_parse/
-- A3_pgltsetup/
-- A4_map/
-- B1_parse/
-- B2_compare/
-- B3_metrics/
-- B4_plots/
-- C1_analysis/
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
sh x_uncompress.sh
```

And compress again:

```{bash}
sh x_compress.sh
```

## Authors
[Bennett D.J.](https://github.com/DomBennett),
[De Palma A.](https://github.com/adrianadepalma),
[Pearse W.D.](https://github.com/willpearse)
and [Purvis A](https://github.com/AndyPurvis).
