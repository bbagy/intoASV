<img src="inst/figures/logo.png" align="left" width="300"/>

<br/><br clear="left"/>

---

**intoASV** is an R framework for quantifying *intra-taxonomic microdiversity*
from 16S rRNA amplicon data using ASV-level sequence variation.

Unlike conventional 16S workflows that focus on community composition,
intoASV enables **strain-level inference**, **cluster diagnostics**, and
**QC-oriented visualizations**, with a particular emphasis on
low-biomass and amplicon-based studies.

---

## Installation

```r
devtools::install_github("bbagy/intoASV", force = T) 
```

---

## Quick example

```r
library(intoASV)
library(phyloseq)

# example dataset included with the package
data(example_phyloseq)

res <- Go_intoASV(
  psIN = example_phyloseq,
  project = "example",
  clustering_cutoff = 0.995,
  similarity_cutoff = 0.97,
  method = "nucdiv",
  aligner = "DECIPHER",
  weighting = "entropy"
)
```

---

## Quality control

intoASV provides multiple QC diagnostics, including:

- ASV clustering diagnostics
- internal distance distributions
- cluster size vs internal divergence
- ASV-level distance heatmaps

Static QC figures are included in the package (`inst/figures`),
while interactive heatmaps are provided as HTML via GitHub Pages
for exploratory analysis.

---

## Documentation

- Full workflow and QC vignette  
  `vignettes/intoASV_QC_heatmap_workflow.Rmd`

---

## Citation

If you use **intoASV** in your research, please cite the accompanying
manuscript (in preparation).
