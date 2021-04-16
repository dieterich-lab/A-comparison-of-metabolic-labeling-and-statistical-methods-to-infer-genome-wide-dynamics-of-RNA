#! /usr/bin/env python3

"""Extract exon-union coordinates from featureCount output
"""

import argparse
import logging

import numpy as np
import pandas as pd

import pbio.utils.bed_utils as bed_utils
import pbio.misc.logging_utils as logging_utils
import pbio.misc.parallel as parallel

logger = logging.getLogger(__name__)


def merge_gene_group(g):

    interval_starts = np.array(g['Start'].split(";")).astype(int)
    # FeatureCounts takes GTF files as an annotation
    # we must add 1 to the ends because gtf is !closed! at the end, but
    # the merge function expects the end to be open
    interval_ends = np.array(g['End'].split(";")).astype(int) + 1

    res = bed_utils.merge_intervals(interval_starts, interval_ends)
    merged_starts, merged_ends, _ = res
    
    # subtract 1 back off the ends
    merged_ends = merged_ends - 1
    
    # convert to bed12
    start = min(merged_starts)
    
    rel_starts = merged_starts - start
    rel_starts_str = ",".join(str(s) for s in rel_starts)
    
    # we now subtract 1 from the start because BED is base-0
    start -= 1
    # we do not subtract 1 from the end, though, because BED
    # is open on the "end"
    
    lengths = merged_ends - merged_starts + 1
    # check 
    if sum(lengths) != g['Length']:
        msg = "Mismatch in exon-union length!"
        logger.error(msg)
    
    lengths = lengths.astype(str)
    lengths_str = ",".join(lengths)
    
    ret = {    
        'seqname': g['Chr'].split(";")[0],
        'start': start,
        'end': max(merged_ends),
        'id': g['Geneid'],
        'score': 0,
        'strand': g['Strand'].split(";")[0],
        'thick_start': -1,
        'thick_end': -1,
        'color': 0,
        'num_exons': len(lengths),
        'exon_lengths': lengths_str,
        'exon_genomic_relative_starts': rel_starts_str,
        # add length for downstream processing
        'length': g['Length']
    }
    
    return ret


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Convert featureCount output to BED12 with exon-union coordinates at
        meta-feature level.""")

    parser.add_argument('tsv', help="The featureCount tsv file")
    parser.add_argument('out', help="The (output) BED12 file, compressed by default")

    parser.add_argument('-p', '--num-cpus', help="The number of CPUs to use",
        type=int, default=12)

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading featureCount tsv file"
    logger.info(msg)
    
    tsv = pd.read_csv(args.tsv, 
                      usecols=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'], 
                      sep='\t', 
                      comment='#')
    
    msg = "Merging..."
    logger.info(msg)
    merged = parallel.apply_parallel(tsv, args.num_cpus, merge_gene_group)  
    merged = pd.DataFrame(merged)
    
    msg = "Sorting..."
    logger.info(msg)
    # We will break ties among transcripts by the order they appear 
    # in the GTF file. This is the same way star breaks ties.
    merged = bed_utils.sort(merged)

    msg = "Writing BED12 to disk"
    logger.info(msg)
    
    fields = bed_utils.bed12_field_names
    fields.append('length')
    bed_utils.write_bed(merged[fields], args.out)


if __name__ == '__main__':
    main()
    
