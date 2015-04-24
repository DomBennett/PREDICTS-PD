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

This will only run stages 1, 3, 4, 5 and 6. Stage 2 can only be run on a HPC
using pG-lt. Stage 7 requires human interaction.

The pipeline requires published trees to have been parsed and placed in
`0_data/parsed_trees/`. If they are not, the setup scripts will need to be run:

```{bash}
sh setup.R
```

Note, stage 3 and the setup script will take a long time to run as they interact
with the [GNR](http://resolver.globalnames.biodinfo.org/).

## Directory structure

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
-- 7_analysis/
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
