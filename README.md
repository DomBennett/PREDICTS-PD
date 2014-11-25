# Response of phylogenetic diversity to human impacts in terrestrial ecosystems

All code is made available for reproducing the analyses as described in [insert paper URL here]. Certain datasets have not been made available due to their large size.

## Pipeline
### run.R
Run all stages

#### 1_pGltsetup.R
Extract names from the PREDICTS database → find parentid → keep studies whose parentid is below the class level → keep studies whose proportion of species:taxon is greater than 0.5, ensure no lists of names are repeated [I’ll need a way to deal with this] → order by ease of running → put in to pG-lt readable format → finish

### 2_pGltrun.R
(outside of R pipeline)

### 3_parse.R
Read in trees from pG-lt runs → rate smooth trees → output

### 4_compare.R
Read in pG-lt and published trees → shrink published trees down to pG-lt trees →  compare using a variety of methods (taxon distance, cophenetic matrix difference) → output

### 5_analysis.R
Read in pG-lt trees, read in PREDICTS data → calculate PD for each site → model PD ~ human impact → output

## Files and folders
### functions/
Folder containing custom functions common to more than one stage

### 0_data
Contains raw data, PREDICTS data [too big for GitHub, see supplementary data] and published trees (bees, birds, mammals ...)

#### 1_pGltsetup
Contains the files and folders exactly as pG-lt accepts them

### 2_pGltrun
Contains the files and folders exactly as pG-lt returns them [too big for GitHub, see supplementary data]

### 3_parse
Ultrametric distributions of trees for each study

### 4_compare
Results from 4_compare.R

### 5_analysis
Results from 5_analysis.R

## Authors
[Bennett D.J.](https://github.com/DomBennett), De Palma A., [Pearse W.D.](https://github.com/willpearse) and Purvis A.
