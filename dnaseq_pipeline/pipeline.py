#! /data/users/dqgu/anaconda3/bin/python
# coding=utf-8
import sys
import os


script_path = os.path.abspath(__file__)
if os.path.islink(script_path):
    script_path = os.readlink(script_path)
sys.path.append(os.path.dirname(os.path.dirname(script_path)))

from dnaseq_pipeline.batch_cmd_generator import NestedCmd
from basic.workflow_basic import basic_arg_parser


parser = basic_arg_parser()
# 可以增加新的流程参数
parser.add_argument('-fastq_info', required=False,
                    help="第一列为样本名,分析结果使用名; "
                         "第二列为read1的fastq路径,如有多个,需要分号分隔; "
                         "第三列为可选, read2的fastq路径,单端测序时则无第三列")
parser.add_argument('-tumour_normal', required=False,
                    help="第一列是tumour sample id, "
                         "第二列是normal sample id, 对配对样本call somatic mutation")
parser.add_argument('-only_tumour', required=False,
                    help="每一行是一个tumour sample id, 对样本call somatic mutation")
parser.add_argument('-germline', required=False,
                    help="每一行是一个normal sample id, 对样本call germline mutation")
parser.add_argument('--disable_markdup_spark',   default=False, action='store_true',
                    help='默认使用MarkDuplicatesSpark, 如果设置该参数，则不使用')
# 收集参数和记录命令行信息
args = parser.parse_args()
args.script_path = script_path

# 可以在这里增加对参数的判断
if args.only_show_steps and args.fastq_info is None:
    args.fastq_info = os.path.join(os.path.dirname(args.script_path), 'example.fastq.info.txt')
    args.pair_info = os.path.join(os.path.dirname(args.script_path), 'example.pair.info.txt')


# 从这里开始写pipeline
def pipeline():
    nc = NestedCmd(args)
    _ = nc.fastqc_raw_data_cmds(args.fastq_info)
    trim_cmds = nc.trimmomatic_cmds(args.fastq_info)
    _ = nc.fastqc_trimmed_data_cmds(trim_cmds)
    if list(trim_cmds.keys())[0].split('_', 1)[0] in nc.workflow_arguments.skip:
        fastq2sam_cmds = nc.RawFastqToSam_cmds(args.fastq_info)
        bwa_cmds = nc.bwa_align_rawdata_cmds(args.fastq_info)
    else:
        fastq2sam_cmds = nc.FastqToSam_cmds(trim_cmds)
        bwa_cmds = nc.bwa_align_trimmed_data_cmds(trim_cmds)
    sam2bam_cmds = nc.sam2bam_cmds(bwa_cmds)
    merge_bam_cmds = nc.MergeBamAlignment_cmds(sam2bam_cmds, fastq2sam_cmds)
    mark_dup_cmds = nc.MarkDuplicates_cmds(merge_bam_cmds)
    sort_bam_cmds = nc.SortAndFixTags_cmds(mark_dup_cmds)
    hs_metric_cmds = nc.CollectHsMetrics_cmds(sort_bam_cmds)
    recalibrator_cmds = nc.BaseRecalibrator_cmds(sort_bam_cmds)
    gather_report_cmds = nc.GatherBQSRReports_cmds(recalibrator_cmds)
    apply_bqsr_cmds = nc.ApplyBQSR_cmds(gather_report_cmds)
    gather_bam_cmds = nc.GatherBamFiles_cmds(apply_bqsr_cmds)
    gather_bam_cmds = nc.samtools_index_cmds(gather_bam_cmds)
    get_pileup_cmds = nc.GetPileupSummaries_cmds(gather_bam_cmds)
    split_intervals_cmds = nc.SplitIntervals_cmds()
    if args.tumour_normal:
        mutect_cmds = nc.MuTect2_cmds(args.tumour_normal, gather_bam_cmds, split_intervals_cmds)
        lrom_cmds = nc.LearnReadOrientationModel_cmds(mutect_cmds)
        merge_vcf_cmds = nc.MergeVcfs_cmds(mutect_cmds)
        gather_mutect_bam_cmds = nc.GatherMutect2Bam_cmds(mutect_cmds)
        sort_sam_cmds = nc.SortSam_cmds(gather_mutect_bam_cmds)
        build_index_cmds = nc.BuildBamIndex_cmds(sort_sam_cmds)
        merge_stats_cmds = nc.MergeMutectStats_cmds(mutect_cmds)
        cal_contaminate_cmds = nc.CalculateContamination_cmds(args.tumour_normal, get_pileup_cmds)
        filter_calls_cmds = nc.FilterMutectCalls_cmds(merge_vcf_cmds, merge_stats_cmds, lrom_cmds, cal_contaminate_cmds)
        filter_artifact_cmds = nc.FilterAlignmentArtifacts_cmds(filter_calls_cmds)

    if args.only_tumour:
        mutect_cmds2 = nc.MuTect2_cmds(args.only_tumour, gather_bam_cmds, split_intervals_cmds, step_name='mutectSingle')
        lrom_cmds2 = nc.LearnReadOrientationModel_cmds(mutect_cmds2, step_name='LROMSingle')
        merge_vcf_cmds2 = nc.MergeVcfs_cmds(mutect_cmds2, step_name='MergeVcfsSingle')
        gather_mutect_bam_cmds2 = nc.GatherMutect2Bam_cmds(mutect_cmds2, step_name='MergeMutect2BamSingle')
        sort_sam_cmds2 = nc.SortSam_cmds(gather_mutect_bam_cmds2, step_name='SortBamSingle')
        build_index_cmds2 = nc.BuildBamIndex_cmds(sort_sam_cmds2, step_name='IndexBamSingle')
        merge_stats_cmds2 = nc.MergeMutectStats_cmds(mutect_cmds2, step_name='MergeVcfStatsSingle')
        contaminate_cmds2 = nc.CalculateContamination_cmds(args.only_tumour, get_pileup_cmds, step_name='CalContaminateSingle')
        filter_calls_cmds2 = nc.FilterMutectCalls_cmds(merge_vcf_cmds2, merge_stats_cmds2,
                                                       lrom_cmds2, contaminate_cmds2,
                                                       step_name='FilterVcfSingle')
        filter_artifact_cmds2 = nc.FilterAlignmentArtifacts_cmds(filter_calls_cmds2, step_name='FilterArtifactsSingle')

    if args.germline:
        haplotype_cmds = nc.HaplotypeCaller_cmds(args.germline, gather_bam_cmds, split_intervals_cmds)
        merge_vcf_cmds3 = nc.MergeVcfs_cmds(haplotype_cmds, step_name='MergeVcfsGermline')

    nc.run()


if __name__ == '__main__':
    pipeline()
