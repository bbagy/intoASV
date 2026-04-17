<img src="inst/figures/logo.png" align="left" width="400"/>

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
devtools::install_github("bbagy/intoASV") 
```

---

## Quick example

**intoASV** takes a `phyloseq` object as input and computes intra-taxonomic
microdiversity metrics from ASV-level sequence variation.

The resulting microdiversity estimates are merged into the `sample_data`
slot of the returned `phyloseq` object, enabling downstream analysis
within standard phyloseq-based workflows.

```r
library(intoASV)
library(phyloseq)

# your phyloseq object containing ASV abundance, taxonomy,
# and DNA sequences in refseq(ps)
ps <- your_phyloseq_object

res <- Go_intoASV(
  psIN = ps,
  project = "myproject",
  level = "Genus",
  taxonomy_cluster_cutoff = 0.995,
  method = "nucdiv",
  aligner = "DECIPHER",
  weighting = "entropy"
)

res_similarity <- Go_intoASV(
  psIN = ps,
  project = "myproject",
  global_similarity_cutoff = 0.97,
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

Typical outputs are written under a date-stamped project directory such as:

```text
myproject_260321/intoASV/pi_tab/
myproject_260321/intoASV/cluster_diagnostics/
```

After similarity-based clustering, cluster diagnostics can be generated with:

```r
Go_clusterDiagnostics(
  project = "myproject",
  cluster_map =
    "myproject_260321/intoASV/pi_tab/cluster_map_similarity_0.970_260321.csv"
)
```

---

## Documentation

- **Full workflow and QC vignette (HTML)**  
  https://bbagy.github.io/intoASV/articles/intoASV_QC_heatmap_workflow.html

- **Interactive ASV distance heatmap (HTML)**  
  https://bbagy.github.io/intoASV/ASV_distance_heatmap_20251226.html

---

## Citation

If you use **intoASV** in your research, please cite the accompanying
manuscript (in preparation).
