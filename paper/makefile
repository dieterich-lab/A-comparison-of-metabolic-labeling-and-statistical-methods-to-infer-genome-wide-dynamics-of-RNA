NOW=`date --utc +%y%m%d.%H%M%SZ`

SLURM=srun

estimates:
	$(SLURM) analysis/get_estimates.R > analysis/$@.$(NOW).out 2>&1 &

# figures
figure3:
	$(SLURM) figures/3/plot.R > figures/3/$@.$(NOW).out 2>&1 &
figure4:
	$(SLURM) figures/4/plot.R > figures/4/$@.$(NOW).out 2>&1 &
# supplementary file "Estimating variance and optimal labeling time using the pulseR workflow"
figureD:
	$(SLURM) figures/D/plot.R > figures/D/$@.$(NOW).out 2>&1 &
	
