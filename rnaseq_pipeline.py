import os
import configparser
import argparse
from cmd_generator import *

"""
steps:
raw_fastqc
trimmomatic
clean_fastqc
star_index
star_align
bam_index
scallop
salmon_index
salmon_quant
"""
arg_ini = 'arguments.ini'
arg_pool = configparser.ConfigParser()
arg_pool.read(arg_ini)
commands = dict()
project_dir = os.getcwd()


def cmd_dict(cpu=1, mem=100*1024*1024, retry=1, monitor_resource=True,
             monitor_time_step=2, check_resource_before_run=True, **kwargs):
    args = locals()
    args.pop('kwargs')
    for each in args:


step = 'fastqc'
fastq_list = list()
raw_qc_dir = 'fastqc'
if not os.path.exists(raw_qc_dir):
    os.mkdir(raw_qc_dir)

for each in fastq_list:
    outdir = os.path.join(project_dir, raw_qc_dir, each)
    if os.path.exists(outdir):
        os.mkdir(outdir)
    cmd = fastqc(**arg_pool[step], outdir=outdir, tmpdir=outdir, fastq=each)
    commands[step+'_'+each] = dict(
        cmd=cmd,
        cpu=1,
        mem=100*1024*1024,
        depend='',
        retry=1,
        monitor_resource=True,
        monitor_time_step=2,
        check_resource_before_run=True,
    )