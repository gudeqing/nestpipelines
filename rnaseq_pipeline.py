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
             monitor_time_step=2, check_resource_before_run=True, cmd=None):
    return locals()


class InputParser(object):
    def parse_fastq(self, fastq_file):
        fastq_info = dict()
        with open(fastq_file) as f:
            for line in f:
                if line.startswith('#') or (not line.strip()):
                    pass
                tmp_list = line.strip().split('\t')
                sample, fqs = tmp_list[0], tmp_list[1:]
                fastq_info.setdefault(sample, list())
                read1_list = [x.strip() for x in fqs[0].split(';')]
                fastq_info[sample].append(read1_list)
                if len(fqs) >= 2:
                    read2_list = [x.strip() for x in fqs[1].split(';')]
                    fastq_info[sample].append(read2_list)
        return fastq_info


step = 'fastqc'
fastq_list = [
    [name, 'fq1', 'fq2'],
    ['fq11', 'fq22']
]
raw_qc_dir = 'fastqc'
if not os.path.exists(raw_qc_dir):
    os.mkdir(raw_qc_dir)

for r1, r2 in fastq_list:
    outdir = os.path.join(project_dir, raw_qc_dir, each)
    if os.path.exists(outdir):
        os.mkdir(outdir)
    fastqs =
    cmd = fastqc(**arg_pool[step], outdir=outdir, tmpdir=outdir, fastq=each)
    commands[step+'_'+each] = cmd_dict(cmd=cmd)
