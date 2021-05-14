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
published manuscript.

