NOW=`date --utc +%y%m%d.%H%M%SZ`

LOG=logs
DEST=tables

# source ~/.virtualenvs/python2/bin/activate
SRC=scripts
prep-gtf:
	srun $(SRC)/dexseq_prepare_annotation2.py -r no -f $(DEST)/GRCh38_85_plus_Spikes_Isabel_flat.gtf ../../../index/GRCh38_85_plus_Spikes_Isabel.gtf $(DEST)/GRCh38_85_plus_Spikes_Isabel_flat.gff > $(LOG)/$@.$(NOW).out 2>&1 &

slam=mapping/slam/*.bam
tuc=mapping/tuc/*.bam
tls=mapping/tls/*.bam

res=tables/$*-0h-1h-2h-4h-8h.counts.tsv
count-%:
	srun --mem=50G -c 16 featureCounts -f -O -s 2 -p -B -C -M --fraction --ignoreDup -T 16 --tmpDir /scratch/global_tmp -a $(DEST)/GRCh38_85_plus_Spikes_Isabel_flat.gtf -o $(res) $($*) > $(LOG)/$@.$(NOW).log 2>&1 &

# DEU per time point, including all protocols (SLAM, TLS, and TUC)
OBJS=1 2 4 8
deu-time:
	for obj in $(OBJS); do \
		srun $(SRC)/run-deu-time.R $$obj > $(LOG)/$@.$(NOW).$$obj.out 2>&1 & \
	done
