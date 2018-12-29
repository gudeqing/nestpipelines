# coding=utf-8
import os


def trimmomatic(**kwargs):
    cmd = '{} '.format(kwargs['java'])
    cmd += '{} '.format(kwargs['trimmomatic'])
    cmd += 'PE '
    cmd += '-threads {}'.format(kwargs['threads'])
    cmd += '{} '.format(kwargs['fq1'])
    cmd += '{} '.format(kwargs['fq2'])
    cmd += '{} '.format(kwargs['trimmed_fq1'])
    cmd += '{} '.format(kwargs['unpaired_fq1'])
    cmd += '{} '.format(kwargs['trimmed_fq2'])
    cmd += '{} '.format(kwargs['unpaired_fq2'])
    cmd += 'ILLUMINACLIP:{}'.format(kwargs['adapter_fasta'])
    cmd += ':{}'.format(kwargs['seed_mismatches'])
    cmd += ':{}'.format(kwargs['palindrome_clip_threshold'])
    cmd += ':{} '.format(kwargs['simple_clip_threshold'])
    cmd += 'LEADING:{} '.format(kwargs['leading'])
    cmd += 'TRAILING:{} '.format(kwargs['trailing'])
    cmd += 'SLIDINGWINDOW:{} '.format(kwargs['sliding_window'])
    cmd += 'MINLEN:{} '.format(kwargs['min_length'])
    cmd += '-trimlog {} '.format(os.path.dirname(kwargs['trimmed_fq1'])+'/trim.log')
    return cmd


def star_index(**kwargs):
    """STAR --runThreadN 30
    --runMode genomeGenerate
    --genomeDir STARINDEX_20180118/
    --genomeFastaFiles WholeGenomeFasta/genome.fa
    --sjdbGTFfile ../Annotation/Genes/genes.gtf
    --sjdbOverhang 134
    """
    cmd = '{} '.format(kwargs['star'])
    cmd += '--runThreadN {} '.format(kwargs['runThreadN'])
    cmd += '--runMode genomeGenerate '
    cmd += '--genomeDir {} '.format(kwargs['genomeDir'])
    cmd += '--genomeFastaFiles {} '.format(kwargs['genomeFastaFiles'])
    cmd += '--sjdbGTFfile {} '.format(kwargs['sjdbGTFfile'])
    cmd += '--sjdbOverhang {} '.format(kwargs['sjdbOverhang'])
    return cmd


def star_align(**kwargs):
    cmd = '{} '.format(kwargs['star'])
    cmd += '--runThreadN {} '.format(kwargs['runThreadN'])
    cmd += '--genomeDir {} '.format(kwargs['genomeDir'])
    cmd += '--readFilesIn {} '.format(kwargs['readFilesIn'])
    cmd += '--outFileNamePrefix  {} '.format(kwargs['outFileNamePrefix '])
    cmd += '--outSAMtype {} '.format(kwargs['outSAMtype'])
    cmd += '--outSAMunmapped {} '.format(kwargs['outSAMunmapped'])
    cmd += '--readFilesCommand {} '.format(kwargs['readFilesCommand'])
    cmd += '--twopassMode {} '.format(kwargs['twopassMode'])
    cmd += '--outFilterMultimapNmax {} '.format(kwargs['outFilterMultimapNmax'])
    cmd += '--alignSJoverhangMin {} '.format(kwargs['alignSJoverhangMin'])
    cmd += '--alignSJDBoverhangMin {} '.format(kwargs['alignSJDBoverhangMin'])
    cmd += '--chimSegmentMin {} '.format(kwargs['chimSegmentMin'])
    cmd += '--outFilterMismatchNoverLmax {} '.format(kwargs['outFilterMismatchNoverLmax'])
    cmd += '--outFilterType {} '.format(kwargs['outFilterType'])
    return cmd


def samtools_index(**kwargs):
    cmd = '{} index '.format(kwargs['samtools'])
    cmd += '{} '.format(kwargs['bam'])
    return cmd


def salmon_index(**kwargs):
    cmd = 'salmon index '.format(kwargs['salmon'])
    cmd += '-t {} '.format(kwargs['transcript_fasta'])
    cmd += '-i {} '.format(kwargs['index_prefix'])
    cmd += '-k {} '.format(kwargs['kmer'])
    return cmd


def salmon_quant(**kwargs):
    cmd = '{} quant '.format(kwargs['salmon'])
    cmd += '-i {} '.format(kwargs['index'])
    cmd += '-l A '
    if kwargs['mode'] == 'pair':
        cmd += '-1 {} -2 {} '.format(kwargs['fq1'], kwargs['fq2'])
    else:
        cmd += '-r {} '.format(kwargs['fq'])
    cmd += '-o {} '.format(kwargs['out_prefix'])
    cmd += '--gcBias '
    cmd += '-p {} '.format(kwargs['threads'])
    if 'transcript2gene' in kwargs:
        if kwargs['transcript2gene'].strip():
            cmd += ' -g {} '.format(kwargs['transcript2gene'])
    return cmd


