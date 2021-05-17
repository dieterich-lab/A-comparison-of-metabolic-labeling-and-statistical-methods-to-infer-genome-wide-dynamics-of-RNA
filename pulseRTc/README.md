# pulseRTc
Split BAM files (RT or nucleotide conversion, *e.g.* SLAM-, TimeLapse-, TUC-seq ) for using as input to pulseR.

### Prerequisites

Pinned version of selected dependencies are listed in the `requirements.txt` file for reproducible installation. The workflow has been tested with Python 3.6. The current version of `pysam` wraps `htslib-1.7`, `samtools-1.7`, and `bcftools-1.6`. A working installation of `samtools` is required to run the workflow. 
Read counting has been done with featureCounts v1.5.1. The pulseR code and analysis scripts are in R (see below for further information).

### Installation

To install the local VCS project in development mode, use the `--editable` or `-e` option, otherwise
this flag can be ignored. 

First create a virtual environment:
 
```
python3 -m venv /path/to/virtual/environment
```

For information about Python virtual environments, see the [venv](https://docs.python.org/3/library/venv.html) documentation.
To activate the new virtual environment and install `pulseRTc`:

```
# Activate the new virtual environment.
source /path/to/virtual/environment/bin/activate

# If necessary, upgrade pip and wheel or additional packages (such as setuptools if installing in editable mode).
pip install --upgrade pip setuptools wheel

# Clone the git repository...
git clone https://github.com/dieterich-lab/ComparisonOfMetabolicLabeling.git
cd ComparisonOfMetabolicLabeling/pulseRTc

# ... or do a sparse checkout e.g.
mkdir MetabolicLabelingWorkflow && cd MetabolicLabelingWorkflow
git init
git remote add -f origin https://github.com/dieterich-lab/ComparisonOfMetabolicLabeling
git config core.sparseCheckout true
echo 'pulseRTc' >> .git/info/sparse-checkout
git pull origin master
cd pulseRTc

# The period is required, it is the local project path
pip --verbose install -r requirements.txt . 2>&1 | tee install.log

```

### Running

To run the workflow, call the main wrapper `run` with a configuration file. See *makefile* for examples.

Scripts and example configuration files are under *run*. Intermediate results are under *workflow*.
The pulseR scripts are under *pulser*, and the final results under *results*.

The pulseR code and analysis scripts are in R. Minimal requirements are


```{r}
install.packages("devtools")
library(devtools)

install_github("dieterich-lab/pulseR", subdir="pkg")
install_cran("tidyverse")

``` 

**Note:** The directory *workflow/mapping* does not contain the alignment files. Mapping is left to the user.
Some large intermediate files are missing from *workflow/mismatches*.






