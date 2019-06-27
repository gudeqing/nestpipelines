# coding = utf-8
import os


def fastqc(**kwargs):
    cmd = '{} '.format(kwargs['fastqc'])
    cmd += '--outdir {} '.format(kwargs['outdir'])
    cmd += '--threads {} '.format(kwargs['threads'])
    cmd += '--dir {} '.format(kwargs['tmpdir'])
    if 'adapters' in kwargs and kwargs['adapters']:
        cmd += '--adapters {} '.format(kwargs['adapters'])
    cmd += '--extract '
    cmd += '{fastqs} '.format(fastqs=kwargs['fastqs'])
    return cmd


def cutadapt(**kwargs):
    cmd = '{} '.format(kwargs['cutadapt'])
    cmd += '-a {} '.format(kwargs['a'])
    cmd += '-g {} '.format(kwargs['g'])
    cmd += '-e {} '.format(kwargs['error-rate'])
    if kwargs['no-indels'] == 'yes':
        cmd += '--no-indels '
    cmd += '-n {} '.format(kwargs['times'])
    cmd += '-O {} '.format(kwargs['overlap'])
    cmd += '-j {} '.format(kwargs['cores'])
    cmd += '-u {} '.format(kwargs['cut_before'])
    if kwargs['--trim-n'] == 'yes':
        cmd += '--trim-n '.format()
    if kwargs['--discard-untrimmed'] == 'yes':
        cmd += '--discard-untrimmed '
    cmd += '-m {} '.format(kwargs['minimum_length'])
    cmd += '-M {} '.format(kwargs['maximum_length'])
    cmd += '-o {} '.format(kwargs['output_read1'])
    if 'output_read2' in kwargs:
        cmd += '-p {} '.format(kwargs['output_read2'])
    cmd += '{} '.format(kwargs['input_read1'])
    if 'input_read2' in kwargs:
        cmd += '{} '.format(kwargs['input_read2'])
    return cmd


def mageck_count(**kwargs):
    cmd = '{} '.format(kwargs['mageck'])
    cmd += '-l {} '.format(kwargs['-l'])
    cmd += '--norm-method {} '.format(kwargs['--norm-method'])
    cmd += '--control-sgrna {} '.format(kwargs['--control-sgrna'])
    if kwargs['--keep-tmp'] == 'yes':
        cmd += '--keep-tmp '
    cmd += '-n {} '.format(kwargs['out_prefix'])
    cmd += '--sample-label {} '.format(kwargs['sample_label'])
    cmd += '--fastq {} '.format(kwargs['fastq_read1'])
    if 'fastq_read2' in kwargs:
        cmd += '--fastq-2 {} '.format(kwargs['fastq_read2'])
    return cmd


def bowtie2_index(**kwargs):
    cmd = '{} '.format(kwargs['bowtie2-build'])
    cmd += '{} '.format(kwargs['fasta'])
    cmd += '{} '.format(kwargs['basename'])
    return cmd


def bowtie2_align(**kwargs):
    cmd = '{} '.format(kwargs['bowtie2'])
    cmd += '-x {} '.format(kwargs['index'])
    cmd += '-p {} '.format(kwargs['threads'])
    if kwargs['--norc'] == 'yes':
        cmd += '--norc '
    if 'fastq_read2' in kwargs:
        cmd += '-1 {} '.format(kwargs['fastq_read1'])
        cmd += '-2 {} '.format(kwargs['fastq_read2'])
    else:
        cmd += '-U {} '.format(kwargs['fastq_read1'])
    cmd += '| {} view -bS - > {}'.format(kwargs['samtools'], kwargs['out_bam'])
    return cmd


