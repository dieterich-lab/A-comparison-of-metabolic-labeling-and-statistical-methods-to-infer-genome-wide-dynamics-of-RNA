## Supplementary code and results

Analysis and plotting scripts to reproduce the published results.
Under *tables* are the combined results obtained with `analysis/get_estimates.R`.
See *makefile* for examples.

### Requirements

```{r}
pkgs <- c(
  "tidyverse",
  "cowplot",
  "grid",
  "gridExtra",
  "openxslx",
  "wordcloud"
)

install.packages("devtools")
library(devtools)

install_github("dieterich-lab/pulseR", subdir="pkg")
install_github("hms-dbmi/UpSetR", subdir="pkg")
install_cran(pkgs)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("topGO")
BiocManager::install("org.Hs.eg.db")

```

The directory *supplement* contains extra supplementary material referenced in the 
published manuscript:

- Supplementary File 1: Mapping quality control
- Supplementary File 2: pulseR quality control (see *pulseRTc/workflow/mismatches/* )
- Supplementary File 3: GRAND-SLAM quality control (see *grand-slam/* )
- Supplementary File 4: pulseR results using a different variant caller
- Supplementary File 5: Estimating variance and optimal labeling time using the pulseR workflow


In addition, under *supplement/dge*, we include a global gene expression analysis between no labeling
and the different labeling durations for the nucleotide conversion protocols. 

Under *supplement/deu*, we include a differential exon usage analysis using the nucleotide conversion protocols.

**Note**: Additional R packages are required to run the DGE and the DEU analyses.

