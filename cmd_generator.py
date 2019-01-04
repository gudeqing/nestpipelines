# coding=utf-8
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


def trimmomatic(**kwargs):
    cmd = '{} -jar '.format(kwargs['java'])
    cmd += '{} '.format(kwargs['trimmomatic'])
    if kwargs['mode'] == 'PE':
        cmd += 'PE '
    else:
        cmd += 'SE '
    cmd += '-threads {} '.format(kwargs['threads'])
    cmd += '{} '.format(kwargs['fq1'])
    if kwargs['mode'] == 'PE':
        cmd += '{} '.format(kwargs['fq2'])
    cmd += '{} '.format(kwargs['trimmed_fq1'])
    if kwargs['mode'] == 'PE':
        cmd += '{} '.format(kwargs['unpaired_fq1'])
        cmd += '{} '.format(kwargs['trimmed_fq2'])
        cmd += '{} '.format(kwargs['unpaired_fq2'])
    cmd += 'ILLUMINACLIP:{}'.format(kwargs['adapter_fasta'])
    cmd += ':{}'.format(kwargs['seed_mismatches'])
    cmd += ':{}'.format(kwargs['palindrome_clip_threshold'])
    cmd += ':{} '.format(kwargs['simple_clip_threshold'])
    cmd += 'LEADING:{} '.format(kwargs['leading'])
    cmd += 'TRAILING:{} '.format(kwargs['trailing'])
    cmd += 'SLIDINGWINDOW:{}:{} '.format(kwargs['sliding_window_size'], kwargs['sliding_window_quality'])
    cmd += 'MINLEN:{} '.format(kwargs['min_length'])
    cmd += '-trimlog {} '.format(os.path.join(os.path.dirname(kwargs['trimmed_fq1']), 'trim.log'))
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
    cmd += '--outFileNamePrefix  {} '.format(kwargs['outFileNamePrefix'])
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
    cmd += '--quantMode {} '.format(kwargs['quantMode'])
    return cmd


def samtools_index(**kwargs):
    cmd = '{} index '.format(kwargs['samtools'])
    cmd += '{} '.format(kwargs['bam'])
    return cmd


def scallop(**kwargs):
    cmd = '{} '.format(kwargs['scallop'])
    cmd += '-i {} '.format(kwargs['bam'])
    cmd += '-o {} '.format(kwargs['out_gtf'])
    cmd += '--library_type {} '.format(kwargs['library_type'])
    cmd += '--min_transcript_coverage {} '.format(kwargs['min_transcript_coverage'])
    cmd += '--min_single_exon_coverage {} '.format(kwargs['min_single_exon_coverage'])
    cmd += '--min_transcript_length_increase {} '.format(kwargs['min_transcript_length_increase'])
    cmd += '--min_transcript_length_base {} '.format(kwargs['min_transcript_length_base'])
    cmd += '--min_mapping_quality {} '.format(kwargs['min_mapping_quality'])
    cmd += '--max_num_cigar {} '.format(kwargs['max_num_cigar'])
    cmd += '--min_bundle_gap {} '.format(kwargs['min_bundle_gap'])
    cmd += '--min_num_hits_in_bundle {} '.format(kwargs['min_num_hits_in_bundle'])
    cmd += '--min_flank_length {} '.format(kwargs['min_flank_length'])
    cmd += '--min_splice_bundary_hits {} '.format(kwargs['min_splice_bundary_hits'])
    return cmd


def salmon_index(**kwargs):
    cmd = '{} index '.format(kwargs['salmon'])
    cmd += '-t {} '.format(kwargs['transcript_fasta'])
    cmd += '-i {} '.format(kwargs['index_prefix'])
    cmd += '-k {} '.format(kwargs['kmer'])
    return cmd


def salmon_quant(**kwargs):
    cmd = '{} quant '.format(kwargs['salmon'])
    cmd += '-i {} '.format(kwargs['index'])
    cmd += '-l A '
    if kwargs['mode'] == 'PE':
        cmd += '-1 {} -2 {} '.format(kwargs['fq1'], kwargs['fq2'])
    else:
        cmd += '-r {} '.format(kwargs['fq'])
    cmd += '-o {} '.format(kwargs['out_prefix'])
    cmd += '--gcBias '
    cmd += '-p {} '.format(kwargs['threads'])
    if 'transcript2gene' in kwargs:
        if kwargs['transcript2gene'].strip():
            cmd += '-g {} '.format(kwargs['transcript2gene'])
    return cmd


def get_new_transcript(**kwargs):
    """
    gffcompare -o gffall -r reference.gtf scallop.gtf
    gtfcuff puniq gffall.scallop.gtf.tmap scallop.gtf reference.gtf unique.gtf
    gffread unique.gtf -g genome -w unique.fa
    cat unique.fa reference.fa > union.fa
    """
    pass


def abundance_estimates_to_matrix(**kwargs):
    pass


class toolbox():
    pass