#! /data/users/dqgu/anaconda3/bin/python
# coding=utf-8
import sys
import os
import time

script_path = os.path.abspath(__file__)
if os.path.islink(script_path):
    script_path = os.readlink(script_path)
sys.path.append(os.path.dirname(os.path.dirname(script_path)))

from dnaseq_pipeline.batch_cmd_generator import NestedCmd
from basic.workflow_basic import basic_arg_parser


if len(sys.argv) <= 1:
    exit('please provide at least one argument, use -h for help')
parser = basic_arg_parser()
# 可以增加新的流程参数
parser.add_argument('-fastq_info', required=False,
                    help="第一列为样本名,分析结果使用名; "
                         "第二列为read1的fastq路径,如有多个,需要分号分隔; "
                         "第三列为可选, read2的fastq路径,单端测序时则无第三列")
parser.add_argument('-pair_info', required=False,
                    help="第一列是tumour sample id, "
                         "第二列是normal sample id,如果没有则第二列为空. 注意: 第一列不能有重复")
parser.add_argument('-germline_for_unpaired', default='no')
parser.add_argument('--disable_markdup_spark',   default=False, action='store_true',
                    help='默认使用MarkDuplicatesSpark, 如果设置该参数，则不使用')
# 收集参数和记录命令行信息
args = parser.parse_args()
args.script_path = script_path
with open("cmd." + str(time.time()) + ".txt", 'w') as f:
    f.write(' '.join(sys.argv) + '\n')
    print(sys.argv)
    f.write('Detail: \n')
    for k, v in args.__dict__.items():
        f.write('{}: {}\n'.format(k, v))

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
    recalibrator_cmds = nc.BaseRecalibrator_cmds(sort_bam_cmds)
    gather_report_cmds = nc.GatherBQSRReports_cmds(recalibrator_cmds)
    apply_bqsr_cmds = nc.ApplyBQSR_cmds(gather_report_cmds)
    gather_bam_cmds = nc.GatherBamFiles_cmds(apply_bqsr_cmds)
    gather_bam_cmds = nc.samtools_index_cmds(gather_bam_cmds)
    get_pileup_cmds = nc.GetPileupSummaries_cmds(gather_bam_cmds)
    split_intervals_cmds = nc.SplitIntervals_cmds()
    mutect_cmds = nc.MuTect2_cmds(args.pair_info, gather_bam_cmds, split_intervals_cmds)
    lrom_cmds = nc.LearnReadOrientationModel_cmds(mutect_cmds)
    merge_vcf_cmds = nc.MergeVcfs_cmds(mutect_cmds)
    gather_mutect_bam_cmds = nc.GatherMutect2Bam_cmds(mutect_cmds)
    sort_sam_cmds = nc.SortSam_cmds(gather_mutect_bam_cmds)
    build_index_cmds = nc.BuildBamIndex_cmds(sort_sam_cmds)
    merge_stats_cmds = nc.MergeMutectStats_cmds(mutect_cmds)
    cal_contaminate_cmds = nc.CalculateContamination_cmds(args.pair_info, get_pileup_cmds)
    filter_calls_cmds = nc.FilterMutectCalls_cmds(merge_vcf_cmds, merge_stats_cmds, lrom_cmds, cal_contaminate_cmds)
    filter_artifact_cmds = nc.FilterAlignmentArtifacts_cmds(filter_calls_cmds)
    if args.germline_for_unpaired != 'no':
        "加入允许call germline 突变的功能"
        pass
    nc.run()


if __name__ == '__main__':
    pipeline()
