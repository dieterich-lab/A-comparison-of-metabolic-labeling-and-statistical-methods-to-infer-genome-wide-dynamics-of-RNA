NOW=`date --utc +%y%m%d.%H%M%SZ`

LOG=logs

# read counting - "un-split" files for DGE
# module load subread, featureCounts v1.5.1
slam=mapping/slam/*.bam
tuc=mapping/tuc/*.bam
tls=mapping/tls/*.bam
res=tables/$*-0h-1h-2h-4h-8h.counts.tsv
count-%:
	srun --mem=50G -c 16 featureCounts -t exon -s 2 -p -B -C -M --fraction --ignoreDup -T 16 --tmpDir /scratch/global_tmp -a ../../../index/GRCh38_85_plus_Spikes_Isabel.gtf -o $(res) $($*) > $(LOG)/$@.$(NOW).log 2>&1 &

SRC=scripts
OBJS=slam tuc tls
run-dge:
	for obj in $(OBJS); do \
		srun $(SRC)/run-deseq.R $$obj > $(LOG)/$@.$(NOW).$$obj.out 2>&1 & \
	done
run-dge-all:
	srun $(SRC)/run-deseq-all.R > $(LOG)/$@.$(NOW).out 2>&1 &
run-cellplot:
	srun $(SRC)/run-cellplot.R > $(LOG)/$@.$(NOW).out 2>&1 &
