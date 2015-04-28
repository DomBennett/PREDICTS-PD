# Response of phylogenetic diversity to human impacts in terrestrial ecosystems

All code is made available for reproducing the analyses as described in
[insert paper URL here]. Certain datasets have not been made available due to
their large size.

## Running

Make sure all relevant packages are installed:

```{bash}
Rscript install_deps.R
```

To run pipeline (UNIX OS):

```{bash}
sh run.sh
```

This will only run stages 1, 3, 4, 5, 6 and 7. Stage 2 can only be run on a HPC
using [pG-lt](https://github.com/DomBennett/pG-lt). Stage 8 requires human
interaction.

The pipeline requires published trees to have been parsed and placed in
`0_data/parsed_trees/`. If they are not, the setup scripts will need to be run:

```{bash}
sh setup.R
```

Note, stage 3 and the setup script will take a long time to run as they interact
with the [GNR](http://resolver.globalnames.biodinfo.org/).

## Setup-stages

1. `setup_presolve.R`: read in published trees and PREDICTS data, search names
against [GNR](http://resolver.globalnames.biodinfo.org/).
2. `setup_parse.R`: read in published trees and rate-smooth

## Stages

1. `pGltsetup.R`: identify suitable PREDICTS studies, output into pG-lt friendly
format
2. `pGltrun.R`: placeholder script -- doesn't do anything
3. `map.R`: use published phylogenies in `0_data` to generate study-level
phylogenies
4. `parse.R`: read in pG-lt and mapped phylogenies, check and rate-smooth
5. `compare.R`: compare mapped and pG-lt phylogenies
6. `metrics.R`: generate PD, PSV and PSE values per site using both pG-lt and
mapped trees.
7. `commplots.R`: community plot PREDICTS data onto pG-lt consensus phylogenies
8. `analysis.R`: analysis script used to analyse phylogeny ~ human disturbance

## Directory structure after run

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
```

## Authors
[Bennett D.J.](https://github.com/DomBennett),
[De Palma A.](https://github.com/adrianadepalma),
[Pearse W.D.](https://github.com/willpearse)
and [Purvis A](https://github.com/AndyPurvis).
