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


def FastqToSam(**kwargs):
    cmd = '{} -jar '.format(kwargs['java'])
    cmd += '{} '.format(kwargs['picard'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += 'F1={} '.format(kwargs['fq'])
    cmd += 'F2={} '.format(kwargs['fq2'])
    cmd += 'O={} '.format(kwargs['out'])
    cmd += 'SM={} '.format(kwargs['sample_name'])
    cmd += 'PL={} '.format(kwargs['PL'])
    cmd += 'READ_GROUP_NAME={} '.format(kwargs['sample_name'])
    cmd += 'LIBRARY_NAME={} '.format(kwargs['sample_name'])
    return cmd


def MergeBamAlignment(**kwargs):
    cmd = '{} -jar '.format(kwargs['java'])
    cmd += '{} '.format(kwargs['picard'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += 'ALIGNED_BAM={} '.format(kwargs['ALIGNED'])
    cmd += 'UNMAPPED_BAM={} '.format(kwargs['UNMAPPED'])
    cmd += 'O={} '.format(kwargs['out'])
    cmd += 'R={} '.format(kwargs['genome_fasta'])
    cmd += 'VALIDATION_STRINGENCY={} '.format(kwargs['VALIDATION_STRINGENCY'])
    cmd += 'EXPECTED_ORIENTATIONS={} '.format(kwargs['EXPECTED_ORIENTATIONS'])
    cmd += 'ATTRIBUTES_TO_RETAIN={} '.format(kwargs['ATTRIBUTES_TO_RETAIN'])
    cmd += 'SORT_ORDER={} '.format(kwargs['SORT_ORDER'])
    cmd += 'IS_BISULFITE_SEQUENCE={} '.format(kwargs['IS_BISULFITE_SEQUENCE'])
    cmd += 'ALIGNED_READS_ONLY={} '.format(kwargs['ALIGNED_READS_ONLY'])
    cmd += 'MAX_RECORDS_IN_RAM={} '.format(kwargs['MAX_RECORDS_IN_RAM'])
    cmd += 'ADD_MATE_CIGAR={} '.format(kwargs['ADD_MATE_CIGAR'])
    cmd += 'MAX_INSERTIONS_OR_DELETIONS={} '.format(kwargs['MAX_INSERTIONS_OR_DELETIONS'])
    cmd += 'PRIMARY_ALIGNMENT_STRATEGY={} '.format(kwargs['PRIMARY_ALIGNMENT_STRATEGY'])
    cmd += 'UNMAPPED_READ_STRATEGY={} '.format(kwargs['UNMAPPED_READ_STRATEGY'])
    cmd += 'UNMAP_CONTAMINANT_READS={} '.format(kwargs['UNMAP_CONTAMINANT_READS'])
    cmd += 'INCLUDE_SECONDARY_ALIGNMENTS={} '.format(kwargs['INCLUDE_SECONDARY_ALIGNMENTS'])
    return cmd


def MarkDuplicates(**kwargs):
    cmd = '{} -jar '.format(kwargs['java'])
    cmd += '{} '.format(kwargs['picard'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += 'I={} '.format(kwargs['input'])
    cmd += 'O={} '.format(kwargs['output'])
    cmd += 'M={} '.format(kwargs['metrics'])
    cmd += 'VALIDATION_STRINGENCY={} '.format(kwargs['VALIDATION_STRINGENCY'])
    cmd += 'CREATE_INDEX={}  '.format(kwargs['CREATE_INDEX'])
    cmd += 'OPTICAL_DUPLICATE_PIXEL_DISTANCE={}  '.format(kwargs['OPTICAL_DUPLICATE_PIXEL_DISTANCE'])
    cmd += 'ASSUME_SORT_ORDER={}  '.format(kwargs['ASSUME_SORT_ORDER'])
    return cmd


def MarkDuplicatesSpark(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += '-I {} '.format(kwargs['input'])
    cmd += '-O {} '.format(kwargs['output'])
    cmd += '-M {} '.format(kwargs['metrics'])
    cmd += '--create-output-bam-index {}  '.format(kwargs['CREATE_INDEX'])
    cmd += '--optical-duplicate-pixel-distance {}  '.format(kwargs['OPTICAL_DUPLICATE_PIXEL_DISTANCE'])
    for each in kwargs['conf'].split():
        cmd += '--conf {} '.format(each)
    return cmd


def bwa_mem(**kwargs):
    cmd = '{} mem '.format(kwargs['bwa'])
    cmd += '-t {} '.format(kwargs['threads'])
    cmd += '-v {} '.format(kwargs['verbosity'])
    if kwargs['-M'] == 'yes':
        cmd += '-M '
    if kwargs['-Y'] == 'yes':
        cmd += '-Y '
    cmd += '-K 100000000 '
    cmd += "-R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA' ".format(sample=kwargs['sample'])
    cmd += '-o {} '.format(kwargs['out'])
    cmd += '{} '.format(kwargs['genome_fasta'])
    cmd += '{} '.format(kwargs['fq1'])
    cmd += '{} '.format(kwargs['fq2'])
    return cmd


def samtools_index(**kwargs):
    cmd = '{} index '.format(kwargs['samtools'])
    cmd += '-@ {} '.format(kwargs['threads'])
    cmd += '{} '.format(kwargs['bam'])
    return cmd


def sam2bam(**kwargs):
    cmd = '{} view '.format(kwargs['samtools'])
    cmd += '-b '
    cmd += '-@ {} '.format(kwargs['threads'])
    cmd += '-o {} '.format(kwargs['out'])
    cmd += '{} '.format(kwargs['sam'])
    return cmd


def SortAndFixTags(**kwargs):
    cmd = '{} -jar '.format(kwargs['java'])
    cmd += '{} '.format(kwargs['picard'])
    cmd += '{} '.format('SortSam')
    cmd += 'I={} '.format(kwargs['input'])
    cmd += 'O=/dev/stdout '
    cmd += 'SORT_ORDER=coordinate '
    cmd += '| '
    cmd += '{} -jar '.format(kwargs['java'])
    cmd += '{} '.format(kwargs['picard'])
    cmd += '{} '.format('SetNmMdAndUqTags')
    cmd += 'I=/dev/stdin '
    cmd += 'O={} '.format(kwargs['output'])
    cmd += 'R={} '.format(kwargs['genome_fasta'])
    cmd += 'CREATE_INDEX={} '.format(kwargs['CREATE_INDEX'])
    cmd += 'CREATE_MD5_FILE={} '.format(kwargs['CREATE_MD5_FILE'])
    return cmd


def FixTags(**kwargs):
    cmd = '{} -jar '.format(kwargs['java'])
    cmd += '{} '.format(kwargs['picard'])
    cmd += '{} '.format('SetNmMdAndUqTags')
    cmd += 'I={} '.format(kwargs['input'])
    cmd += 'O={} '.format(kwargs['output'])
    cmd += 'R={} '.format(kwargs['genome_fasta'])
    cmd += 'CREATE_INDEX={} '.format(kwargs['CREATE_INDEX'])
    cmd += 'CREATE_MD5_FILE={} '.format(kwargs['CREATE_MD5_FILE'])
    return cmd


def BaseRecalibrator(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += '--use-original-qualities {} '.format(kwargs['use-original-qualities'])
    cmd += '-I {} '.format(kwargs['input'])
    cmd += '-O {} '.format(kwargs['output'])
    cmd += '-R {} '.format(kwargs['genome_fasta'])
    for each in kwargs['intervals']:
        cmd += '-L {} '.format(each)
    for each in kwargs['known_sites'].split():
        cmd += '--known-sites {} '.format(each)
    return cmd


def GatherBQSRReports(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    for each in kwargs['report_list']:
        cmd += '-I {} '.format(each)
    cmd += '-O {} '.format(kwargs['output'])
    return cmd


def ApplyBQSR(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += '--use-original-qualities {} '.format(kwargs['use-original-qualities'])
    cmd += '--create-output-bam-md5 {} '.format(kwargs['create-output-bam-md5'])
    cmd += '--add-output-sam-program-record true '
    cmd += '-I {} '.format(kwargs['input'])
    for each in kwargs['intervals']:
        cmd += '-L {} '.format(each)
    cmd += '--bqsr-recal-file {} '.format(kwargs['bqsr-recal-file'])
    cmd += '-O {} '.format(kwargs['output'])
    cmd += '-R {} '.format(kwargs['genome_fasta'])
    for each in kwargs['static-quantized-quals'].split():
        cmd += '--static-quantized-quals {} '.format(each)
    return cmd


def GatherBamFiles(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += '--CREATE_INDEX {} '.format(kwargs['CREATE_INDEX'])
    cmd += '--CREATE_MD5_FILE {} '.format(kwargs['CREATE_MD5_FILE'])
    cmd += '-R {} '.format(kwargs['genome_fasta'])
    for each in kwargs['bam_list']:
        cmd += '-I {} '.format(each)
    cmd += '-O {} '.format(kwargs['output'])
    return cmd


def MuTect2(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += '--bam-output {} '.format(kwargs['bam-output'])
    cmd += '--bam-writer-type {} '.format(kwargs['bam-writer-type'])
    cmd += '--f1r2-tar-gz {} '.format(kwargs['f1r2-tar-gz'])
    cmd += '--germline-resource {} '.format(kwargs['germline-resource'])
    cmd += '--intervals {} '.format(kwargs['intervals'])
    cmd += '--downsampling-stride {} '.format(kwargs['downsampling-stride'])
    cmd += '--max-reads-per-alignment-start {} '.format(kwargs['max-reads-per-alignment-start'])
    cmd += '--max-suspicious-reads-per-alignment-start {} '.format(kwargs['max-suspicious-reads-per-alignment-start'])
    if kwargs['panel-of-normals'].strip():
        cmd += '--panel-of-normals {} '.format(kwargs['panel-of-normals'])
    if kwargs['mitochondria-mode'] == 'true':
        cmd += '--mitochondria-mode {} '.format(kwargs['mitochondria-mode'])
        cmd += '--median-autosomal-coverage 30 '
    cmd += '-R {} '.format(kwargs['genome_fasta'])
    cmd += '-I {} '.format(kwargs['tumour_bam'])
    # cmd += '--tumour-sample {} '.format(kwargs['tumour-sample'])
    # 如果不输入normal bam参数, 则进入tumour only mode
    if kwargs['normal_bam'].endswith('.bam'):
        cmd += '-I {} '.format(kwargs['normal_bam'])
        cmd += '--normal-sample {} '.format(kwargs['normal-sample'])
    cmd += '-O {} '.format(kwargs['output'])
    return cmd


def GetPileupSummaries(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += '-R {} '.format(kwargs['genome_fasta'])
    cmd += '-I {} '.format(kwargs['bam'])
    cmd += '-O {} '.format(kwargs['output'])
    cmd += '-L {} '.format(kwargs['vcf_contamination'])
    cmd += '-V {} '.format(kwargs['vcf_contamination'])
    return cmd


def SplitIntervals(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += '-R {} '.format(kwargs['genome_fasta'])
    cmd += '--intervals {} '.format(kwargs['intervals'])
    cmd += '--scatter-count {} '.format(kwargs['scatter-count'])
    cmd += '-O {} '.format(kwargs['output'])
    return cmd


def MergeVcfs(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    for each in kwargs['input']:
        cmd += '-I {} '.format(each)
    cmd += '-O {} '.format(kwargs['output'])
    return cmd


def SortSam(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += '--SORT_ORDER {} '.format(kwargs['SORT_ORDER'])
    cmd += '-VALIDATION_STRINGENCY {} '.format(kwargs['VALIDATION_STRINGENCY'])
    cmd += '-I {} '.format(kwargs['input'])
    cmd += '-O {} '.format(kwargs['output'])
    return cmd


def BuildBamIndex(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += '-VALIDATION_STRINGENCY {} '.format(kwargs['VALIDATION_STRINGENCY'])
    cmd += '-I {} '.format(kwargs['input'])
    return cmd


def MergeMutectStats(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    for each in kwargs['stats']:
        cmd += '--stats {} '.format(each)
    cmd += '-O {} '.format(kwargs['output'])
    return cmd


def GatherPileupSummaries(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += '--sequence-dictionary {} '.format(kwargs['sequence-dictionary'])
    for each in kwargs['input']:
        cmd += '-I {} '.format(each)
    cmd += '-O {} '.format(kwargs['output'])
    return cmd


def LearnReadOrientationModel(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    for each in kwargs['input']:
        cmd += '-I {} '.format(each)
    cmd += '-O {} '.format(kwargs['output'])
    return cmd


def CalculateContamination(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += '-I {} '.format(kwargs['tumor_pileups'])
    cmd += '-O {} '.format(kwargs['output'])
    cmd += '--tumor-segmentation {} '.format(kwargs['tumor-segmentation'])
    if kwargs['normal_pileups'] != 'dynamic':
        cmd += '--matched-normal {} '.format(kwargs['normal_pileups'])
    return cmd


def FilterMutectCalls(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += '-R {} '.format(kwargs['genome_fasta'])
    cmd += '-V {} '.format(kwargs['variant'])
    cmd += '-O {} '.format(kwargs['output'])
    if kwargs['contamination-table'] != 'dynamic':
        cmd += '--contamination-table {} '.format(kwargs['contamination-table'])
    if kwargs['tumor-segmentation'] != 'dynamic':
        cmd += '--tumor-segmentation {} '.format(kwargs['tumor-segmentation'])
    cmd += '--ob-priors {} '.format(kwargs['ob-priors'])
    cmd += '--stats {} '.format(kwargs['stats'])
    cmd += '--filtering-stats {} '.format(kwargs['filtering-stats'])
    return cmd


def FilterAlignmentArtifacts(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += '-V {} '.format(kwargs['variant'])
    cmd += '-I {} '.format(kwargs['bam'])
    cmd += '--bwa-mem-index-image {} '.format(kwargs['bwa-mem-index-image'])
    cmd += '-O {} '.format(kwargs['output'])
    return cmd


def HaplotypeCaller(**kwargs):
    cmd = '{} '.format(kwargs['gatk'])
    cmd += '{} '.format(kwargs['tool'])
    cmd += '-R {} '.format(kwargs['genome_fasta'])
    cmd += '-I {} '.format(kwargs['input'])
    cmd += '-O {} '.format(kwargs['output'])
    cmd += '-ERC {} '.format(kwargs['emit-ref-confidence'])
    cmd += '-bamout {} '.format(kwargs['bamout'])
    cmd += '--bam-writer-type {} '.format(kwargs['bam-writer-type'])
    cmd += '-contamination {} '.format(kwargs['contamination'])
    cmd += '--intervals {} '.format(kwargs['intervals'])
    for each in kwargs['annotation-group'].split():
        cmd += '-G {} '.format(each)
    return cmd
