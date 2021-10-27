## Supplementary resources 

> Etienne Boileau, Janine Altmüller, Isabel S Naarmann-de Vries, Christoph Dieterich.
> A comparison of metabolic labeling and statistical methods to infer genome-wide dynamics of RNA turnover.
> Briefings in Bioinformatics, 2021, bbab219, https://doi.org/10.1093/bib/bbab219

## How to use this repository

You can clone the repository and prune it locally. You might want to change the structure *e.g* if
you want to perform the analysis of a single SLAM experiment. Some core scripts do not need to be modified such
as the **pulseRTc** package, or the **pulseR** main scripts. However, you may have to modify other scripts,
including some utility R definitions, and makefiles.

First prepare BAM files, `make get-sorted-bams`, and `make get-dedups`. Modify scripts where required.
Prepare the GRAND-SLAM index files, and run GRAND-SLAM (or bcftools) to obtain SNPs. 
Finally, run the pulseRTc workflow `make run-workflow`, `make count-`. Modify makefile where required.
Run pulseR (modify `prep.R`, definition of time sets in `utils.R`), then `make prep`, and `make pulse-`.

We hope to be able to provide in a near future a streamlined workflow and/or more general scripts.

For more information, see below.

### pulseRTc

Workflow and results (pulseR)

### grand-slam

GRAND-SLAM results

### index

Gene annotation

### paper

Analysis and plotting scripts, additional results used in the manuscript. 
In particular, the **Supplementary Files** and **Other analyses** referenced in the manuscript are under [paper/supplement/](paper/supplement/).
