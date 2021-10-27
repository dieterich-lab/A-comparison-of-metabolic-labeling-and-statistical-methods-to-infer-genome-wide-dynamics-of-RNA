NOW=`date --utc +%y%m%d.%H%M%SZ`

SLURM=srun
SLURM24=srun -c 24 --mem=50GB
SLURM40=srun -c 40 --mem=50GB

OBJS=slam tuc tls

# Adjust paths as required

# pre-processing ---
# sort/index bams, rename 
SRC=run
get-sorted-bams:
	srun --mem=40GB $(SRC)/sortnrename $(SRC)/SAMPLELIST.TXT > ../logs/$@.$(NOW).log 2>&1 &

# remove duplicates, re-index
get-dedups:
	srun --mem=40GB $(SRC)/markdups $(SRC)/SAMPLELIST.TXT > ../logs/$@.$(NOW).log 2>&1 &
	
# workflow ---
# source /path/to/virtual/environment/bin/activate
CFG=run
LOG=run/logs
run-workflow:
	for obj in $(OBJS); do \
		$(SLURM) run $(CFG)/$$obj.yaml --trim5p 5 --trim3p 5 --use-slurm --num-cpus 40 --mem 100G --logging-level DEBUG --log-file $(LOG)/$@.$(NOW).$$obj.out --mail-user boileau@uni-heidelberg.de & \
	done
	
# read counting
# module load subread, featureCounts v1.5.1
std=workflow/mapping/std/*.bam
slam=workflow/mapping/slam/*.bam
tuc=workflow/mapping/tuc/*.bam
tls=workflow/mapping/tls/*.bam
res=workflow/tables/$*-0h-1h-2h-4h-8h.counts.tsv
count-%:
	srun --mem=50G -c 16 featureCounts -t exon -s 2 -p -B -C -M --fraction --ignoreDup -T 16 --tmpDir /scratch/global_tmp -a ../index/GRCh38_85_plus_Spikes_Isabel.gtf -o $(res) $($*) > $(LOG)/$@.$(NOW).log 2>&1 &

# pulseR ---
PSRC=pulser
PLOG=pulser/logs
# prepare count tables and metadata for processing with pulseR
# *** for ERCC, user must symlink std files  
# *** some hard coded features in prep.R, e.g. path to files starting with "workflow" (see above), etc.
# *** scripts use package here

prep:
	${SLURM} $(PSRC)/prep.R	> $(PLOG)/$@.$(NOW).log 2>&1 &
	
# all time points, standard using ERCC
pulse-ercc:
	$(SLURM24) $(PSRC)/fit_std.R 1 0 > $(PLOG)/$@.$(NOW).out 2>&1 &
# subset of time points, standard using ERCC
pulse-ercc-sets:
	$(SLURM24) $(PSRC)/fit_std.R 1 1 > $(PLOG)/$@.$(NOW).out 2>&1 &

# all time points, standard fitting normFactors
pulse-std:
	$(SLURM24) $(PSRC)/fit_std.R > $(PLOG)/$@.$(NOW).out 2>&1 &	
# subset of time points, standard fitting normFactors	
pulse-std-sets:
	$(SLURM24) $(PSRC)/fit_std.R 0 1 > $(PLOG)/$@.$(NOW).out 2>&1 &
	
# all time points, SLAM/TLS/TUC
pulse-rtc:
	for obj in $(OBJS); do \
		$(SLURM24) $(PSRC)/fit_rtc.R $$obj > $(PLOG)/$@.$(NOW).$$obj.out 2>&1 & \
	done
# subset of time points, SLAM/TLS/TUC
pulse-rtc-sets:
	for obj in $(OBJS); do \
		$(SLURM24) $(PSRC)/fit_rtc.R $$obj 1 > $(PLOG)/$@.$(NOW).$$obj.out 2>&1 & \
	done
