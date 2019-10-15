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
    if 'head_crop' in kwargs:
        cmd += 'HEADCROP:{} '.format(kwargs['head_crop'])
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
    cmd += '--outFileNamePrefix {} '.format(kwargs['outFileNamePrefix'])
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
    cmd += '--outSAMstrandField {} '.format(kwargs['outSAMstrandField'])
    cmd += '--quantMode {} '.format(kwargs['quantMode'])
    cmd += '--outSAMattrRGline {} '.format(kwargs['outSAMattrRGline'])
    cmd += '--limitBAMsortRAM {} '.format(kwargs['limitBAMsortRAM'])
    cmd += '--limitIObufferSize {} '.format(kwargs['limitIObufferSize'])
    cmd += '--outFilterMatchNminOverLread {} '.format(kwargs['outFilterMatchNminOverLread'])
    cmd += '--limitOutSAMoneReadBytes {} '.format(kwargs['limitOutSAMoneReadBytes'])
    cmd += '--outSAMattrIHstart {} '.format(kwargs['outSAMattrIHstart'])
    cmd += '--alignMatesGapMax {} '.format(kwargs['alignMatesGapMax'])
    cmd += '--alignIntronMax {} '.format(kwargs['alignIntronMax'])
    cmd += '--alignSJstitchMismatchNmax {} '.format(kwargs['alignSJstitchMismatchNmax'])
    cmd += '--chimJunctionOverhangMin {} '.format(kwargs['chimJunctionOverhangMin'])
    cmd += '--chimMultimapScoreRange {} '.format(kwargs['chimMultimapScoreRange'])
    cmd += '--chimSegmentReadGapMax {} '.format(kwargs['chimSegmentReadGapMax'])
    cmd += '--chimMultimapNmax {} '.format(kwargs['chimMultimapNmax'])
    cmd += '--chimNonchimScoreDropMin {} '.format(kwargs['chimNonchimScoreDropMin'])
    cmd += '--chimOutJunctionFormat {} '.format(kwargs['chimOutJunctionFormat'])
    cmd += '--peOverlapNbasesMin {} '.format(kwargs['peOverlapNbasesMin'])
    cmd += '--peOverlapMMp {} '.format(kwargs['peOverlapMMp'])
    cmd += '--chimOutType {} '.format(kwargs['chimOutType'])
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
    cmd += '--validateMappings '
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


def merge_scallop_transcripts(**kwargs):
    """get all transcripts including new ones"""
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += '-gtf {} '.format(kwargs['gtf'])
    cmd += '-refgtf {} '.format(kwargs['refgtf'])
    cmd += '-refgenome {} '.format(kwargs['refgenome'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-reftrans {} '.format(kwargs['reftrans'])
    cmd += '-gtfmerge {} '.format(kwargs['gtfmerge'])
    cmd += '-gffcompare {} '.format(kwargs['gffcompare'])
    cmd += '-gtfcuff {} '.format(kwargs['gtfcuff'])
    cmd += '-gffread {} '.format(kwargs['gffread'])
    return cmd


def abundance_estimates_to_matrix(**kwargs):
    cmd = '{} '.format(kwargs['perl'])
    cmd += '{} '.format(kwargs['program'])
    cmd += '--est_method {} '.format(kwargs['est_method'])
    cmd += '--cross_sample_norm {} '.format(kwargs['cross_sample_norm'])
    cmd += '--name_sample_by_basedir '
    cmd += '--gene_trans_map none '
    cmd += '--out_prefix {} '.format(kwargs['out_prefix'])
    cmd += '{} '.format(kwargs['quant_result'])
    return cmd


def diff_exp(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += '-count {} '.format(kwargs['count_matrix'])
    cmd += '-exp {} '.format(kwargs['exp_matrix'])
    cmd += '-group {} '.format(kwargs['group_info'])
    cmd += '-cmp {} '.format(kwargs['comparison_info'])
    cmd += '-method {} '.format(kwargs['method'])
    cmd += '-output {} '.format(kwargs['result_dir'])
    cmd += '-pool {} '.format(kwargs['threads'])
    cmd += '-pvalue {} '.format(kwargs['pvalue'])
    cmd += '-fc {} '.format(kwargs['fc'])
    cmd += '-sig_type {} '.format(kwargs['sig_type'])
    # cmd += '-padjust_way {} '.format(kwargs['padjust_way'])
    cmd += '--dispersion {} '.format(kwargs['dispersion'])
    if kwargs['plot'].lower() in ['yes', 'true']:
        cmd += '--plot '
    return cmd


def goatools(**kwargs):
    cmd = '{} {} enrich '.format(kwargs['python'], kwargs['script'])
    cmd += '-alpha {} '.format(kwargs['alpha'])
    cmd += '-correct {} '.format(kwargs['correct'])
    cmd += '-goea_out {} '.format(kwargs['goea_out'])
    cmd += '-dag_out {} '.format(kwargs['dag_out'])
    cmd += '-obo {} '.format(kwargs['obo'])
    cmd += '-study {} '.format(kwargs['study'])
    cmd += '-population {} '.format(kwargs['population'])
    cmd += '-gene2go {} '.format(kwargs['gene2go'])
    cmd += '-geneid2symbol {} '.format(kwargs['geneid2symbol'])
    cmd += '-top {} '.format(kwargs['top'])
    cmd += '-dpi {} '.format(kwargs['dpi'])
    return cmd


def kegg_enrich(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += '-deg {} '.format(kwargs['deg'])
    cmd += '-bgn {} '.format(kwargs['bgn'])
    cmd += '-g2k {} '.format(kwargs['g2k'])
    cmd += '-brite {} '.format(kwargs['brite'])
    cmd += '-g2p {} '.format(kwargs['g2p'])
    # cmd += '-k2p {} '.format(kwargs['k2p'])
    cmd += '-k2e {} '.format(kwargs['k2e'])
    cmd += '-dn {} '.format(kwargs['dn'])
    cmd += '--FDR '
    cmd += '-o {} '.format(kwargs['outdir'])
    if kwargs['only_consider_path_annotated_genes'] == 'yes':
        cmd += '--only_consider_path_annotated_genes '
    return cmd


def star_fusion(**kwargs):
    cmd = '{} '.format(kwargs['star_fusion'])
    cmd += '--genome_lib_dir {} '.format(kwargs['genome_lib_dir'])
    cmd += '--chimeric_junction {} '.format(kwargs['Chimeric_out_junction'])
    cmd += '--CPU {} '.format(kwargs['CPU'])
    cmd += '--FusionInspector {} '.format(kwargs['FusionInspector'])
    cmd += '--left_fq {} '.format(kwargs['left_fq'])
    cmd += '--right_fq {} '.format(kwargs['right_fq'])
    if kwargs['examine_coding_effect'].strip().lower() == 'yes':
        cmd += '--examine_coding_effect '
    if kwargs['extract_fusion_reads'].strip().lower() == 'yes':
        cmd += '--extract_fusion_reads '
    cmd += '--output_dir {} '.format(kwargs['outdir'])
    return cmd


def gene_body_coverage(**kwargs):
    cmd = '{} '.format(kwargs['script'])
    cmd += '-i {} '.format(kwargs['bam'])
    cmd += '-r {} '.format(kwargs['bed'])
    cmd += '-l {} '.format(kwargs['min_mRNA_length'])
    cmd += '-f {} '.format(kwargs['image_format'])
    cmd += '-o {} '.format(kwargs['out_prefix'])
    return cmd


def inner_distance(**kwargs):
    cmd = '{} '.format(kwargs['script'])
    cmd += '-i {} '.format(kwargs['bam'])
    cmd += '-r {} '.format(kwargs['bed'])
    cmd += '-o {} '.format(kwargs['out_prefix'])
    cmd += '-k {} '.format(kwargs['sample_size'])
    cmd += '-q {} '.format(kwargs['min_mapping_quality'])
    return cmd


def read_distribution(**kwargs):
    cmd = '{} '.format(kwargs['script'])
    cmd += '-i {} '.format(kwargs['bam'])
    cmd += '-r {} '.format(kwargs['bed'])
    cmd += '> {} '.format(kwargs['outfile'])
    return cmd


def read_duplication(**kwargs):
    cmd = '{} '.format(kwargs['script'])
    cmd += '-i {} '.format(kwargs['bam'])
    cmd += '-q {} '.format(kwargs['min_mapping_quality'])
    cmd += '-o {} '.format(kwargs['out_prefix'])
    return cmd


def rna_fragment_size(**kwargs):
    cmd = '{} '.format(kwargs['script'])
    cmd += '-i {} '.format(kwargs['bam'])
    cmd += '-r {} '.format(kwargs['bed'])
    cmd += '-q {} '.format(kwargs['min_mapping_quality'])
    cmd += '> {} '.format(kwargs['outfile'])
    return cmd


def rpkm_saturation(**kwargs):
    cmd = '{} '.format(kwargs['script'])
    cmd += '-i {} '.format(kwargs['bam'])
    cmd += '-o {} '.format(kwargs['out_prefix'])
    cmd += '-r {} '.format(kwargs['bed'])
    if kwargs['strand_rule'].lower() != 'none':
        cmd += '-d {} '.format(kwargs['strand_rule'])
    cmd += '-c {} '.format(kwargs['rpkm_cutoff'])
    cmd += '-q {} '.format(kwargs['min_mapping_quality'])
    return cmd


def tpm_saturation(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += '-bam {} '.format(kwargs['bam'])
    cmd += '-annotation {} '.format(kwargs['annotation'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-step {} '.format(kwargs['step'])
    cmd += '-samtools {} '.format(kwargs['samtools'])
    cmd += '-outlier_limit {} '.format(kwargs['outlier_limit'])
    cmd += '-threads {} '.format(kwargs['threads'])
    cmd += '-pool_size {} '.format(kwargs['pool_size'])
    cmd += '-featureCounts {} '.format(kwargs['featureCounts'])
    if kwargs['paired'] == 'yes':
        cmd += '--paired '
    return cmd


def exp_analysis(**kwargs):
    cmd = '{} '.format(kwargs['perl'])
    cmd += "{} ".format(kwargs['trinity_ptr'])
    cmd += "--matrix {} ".format(kwargs['matrix'])
    cmd += "--output {} ".format(kwargs['out_prefix'])
    cmd += " --min_gene_prevalence {} ".format(kwargs['min_gene_prevalence'])
    cmd += "--min_gene_expr_val {} ".format(kwargs['min_gene_expr_val'])
    cmd += "--sample_cor_matrix "
    cmd += "--sample_cor {} ".format(kwargs['sample_cor'])
    cmd += "--sample_cor_scale_limits 0.3,1 "
    cmd += "--boxplot_log2_dist {} ".format(kwargs['boxplot_log2_dist'])
    cmd += "--log2 "
    cmd += "--prin_comp 2"
    return cmd


def get_alignment_summary(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += '-bam {} '.format(kwargs['bam'])
    cmd += '-bed {} '.format(kwargs['bed'])
    cmd += '-rRNA_bed {} '.format(kwargs['rRNA_bed'])
    cmd += '-overlap {} '.format(kwargs['overlap'])
    cmd += '-rRNA_overlap {} '.format(kwargs['rRNA_overlap'])
    cmd += '-bedtools {} '.format(kwargs['bedtools'])
    cmd += '-samtools {} '.format(kwargs['samtools'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-threads {} '.format(kwargs['threads'])
    cmd += '-cov_limit {} '.format(kwargs['cov_limit'])
    cmd += '-step {} '.format(kwargs['step'])
    return cmd


def chromosome_read_distribution(**kwargs):
    cmd = '{} idxstats '.format(kwargs['samtools'])
    cmd += '{} > {}'.format(kwargs['bam'], kwargs['outfile'])
    return cmd


def CollectAlignmentSummaryMetrics(**kwargs):
    "https://broadinstitute.github.io/picard/command-line-overview.html#CollectAlignmentSummaryMetrics"
    cmd = '{} -jar '.format(kwargs['java'])
    cmd += '{} '.format(kwargs['picard'])
    cmd += 'CollectAlignmentSummaryMetrics '
    cmd += 'R={} '.format(kwargs['genome'])
    cmd += 'I={} '.format(kwargs['bam'])
    cmd += 'O={} '.format(kwargs['outfile'])
    return cmd


def CollectInsertSizeMetrics(**kwargs):
    cmd = '{} -jar '.format(kwargs['java'])
    cmd += '{} '.format(kwargs['picard'])
    cmd += 'CollectInsertSizeMetrics '
    cmd += 'I={} '.format(kwargs['bam'])
    cmd += 'O={} '.format(kwargs['outfile'])
    cmd += 'H={} '.format(kwargs['outimage'])
    cmd += 'M={} '.format(kwargs['min_percent'])
    return cmd


def CollectTargetedPcrMetrics(**kwargs):
    cmd = '{} -jar '.format(kwargs['java'])
    cmd += '{} '.format(kwargs['picard'])
    cmd += 'CollectTargetedPcrMetrics  '
    cmd += 'I={} '.format(kwargs['bam'])
    cmd += 'O={} '.format(kwargs['outfile'])
    cmd += 'R={} '.format(kwargs['genome'])
    cmd += 'AMPLICON_INTERVALS={} '.format(kwargs['amplicon_interval_list'])
    cmd += 'TARGET_INTERVALS={} '.format(kwargs['targets_interval_list'])
    cmd += 'PER_TARGET_COVERAGE={} '.format(kwargs['per_target_coverage_outfile'])
    return cmd


def CollectRnaSeqMetrics(**kwargs):
    cmd = '{} -jar '.format(kwargs['java'])
    cmd += '{} '.format(kwargs['picard'])
    cmd += 'CollectRnaSeqMetrics  '
    cmd += 'I={} '.format(kwargs['bam'])
    cmd += 'O={} '.format(kwargs['outfile'])
    cmd += 'CHART={} '.format(kwargs['outimage'])
    cmd += 'REF_FLAT={} '.format(kwargs['ref_flat'])
    cmd += 'STRAND={} '.format(kwargs['strand'])
    cmd += 'RIBOSOMAL_INTERVALS={} '.format(kwargs['ribosomal_interval_list'])
    return cmd


def arriba(**kwargs):
    cmd = "{} ".format(kwargs['arriba'])
    cmd += '-x {} '.format(kwargs['x'])
    cmd += '-g {} '.format(kwargs['g'])
    cmd += '-a {} '.format(kwargs['a'])
    cmd += '-b {} '.format(kwargs['b'])
    cmd += '-o {} '.format(kwargs['o'])
    cmd += '-O {} '.format(kwargs['discarded'])
    cmd += '-s {} '.format(kwargs['strandness'])
    cmd += '-T '.format(kwargs['assemble_transcript'])
    cmd += '-P '.format(kwargs['translate_pipetide'])
    if kwargs['known_fusion']:
        cmd += '-k '.format(kwargs['known_fusion'])
    return cmd
