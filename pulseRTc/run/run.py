#! /usr/bin/env python3

"""Wrapper function for splbam. 
"""

import os
import argparse
import logging
import shlex
import yaml

import run.utils as utils

logger = logging.getLogger(__name__)


default_num_cpus = 1
default_mem = '80G'


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Wrapper function for splbam.py""")

    parser.add_argument('config', help="The yaml configuration file (full path).")

    parser.add_argument('-q', '--base-qual', help="The minimum base quality for any given mismatch (default: 20).",
                        type=int, default=20)
    
    parser.add_argument('--trim5p', help="The number bases to trim at the 5' ends of reads (default: 0).",
                        type=int, default=0)
    
    parser.add_argument('--trim3p', help="The number bases to trim at the 3' ends of reads (default: 0).",
                        type=int, default=0)
    
    parser.add_argument('--overwrite', help="""If this flag is present, existing files 
                        will be overwritten.""", action='store_true')

    utils.add_sbatch_options(parser, num_cpus=default_num_cpus, mem=default_mem)
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    # check that all of the necessary programs are callable
    programs = ['splbam']
    utils.check_programs_exist(programs)
    
    required_keys = ['samples',
                     'mapping',
                     'mismatch']
    utils.check_keys_exist(config, required_keys)

    # handle all option strings to call splbam
    # use defaults for most options (not passed)
    logging_str = utils.get_logging_options_string(args)
    
    all_opt_str = "-q {} --trim5p {} --trim3p {}".format(args.base_qual, args.trim5p, args.trim3p)
   
    vcf_str = ""
    if 'vcf' in config.keys() and config['vcf'] == True:
        vcf_str = "--vcf"
     
    snpdata_str = ""
    if 'snpdata' in config.keys():
        snpdata_str = "-s {} {}".format(config['snpdata'], vcf_str)
    
    # handle do_not_call so that we do call splbam, but that it does not run anything
    call = not args.do_not_call
    do_not_call_str = ""
    if not call:
        do_not_call_str = "--do-not-call"
    args.do_not_call = False

    overwrite_str = ""
    if args.overwrite:
        overwrite_str = "--overwrite"
        
    mem_str = "--mem {}".format(shlex.quote(args.mem))

    for name, bam in config['samples'].items():
        
        cmd = 'splbam {} {} {} {} {} {} --num-cpus {} {} {} {} {}'.format(bam,
                                                                          config['mapping'],
                                                                          config['mismatch'],
                                                                          name,
                                                                          all_opt_str,
                                                                          snpdata_str,
                                                                          args.num_cpus,
                                                                          mem_str,
                                                                          logging_str,
                                                                          overwrite_str,
                                                                          do_not_call_str)
        utils.check_sbatch(cmd, args=args)
    
    
if __name__ == '__main__':
    main()   
