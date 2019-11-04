#! /data/users/dqgu/anaconda3/bin/python
# coding=utf-8
import os
import sys
import time


script_path = os.path.abspath(__file__)
if os.path.islink(script_path):
    script_path = os.readlink(script_path)
sys.path.append(os.path.dirname(os.path.dirname(script_path)))

from rnaseq_pipeline.batch_cmd_generator import NestedCmd
from basic.workflow_basic import basic_arg_parser

parser = basic_arg_parser()
parser.add_argument('-bam', required=False, help="只有一列的文件, 每行记录一个bam的绝对路径, 且bam已经有索引")
parser.add_argument('-fastq_info', required=False,
                    help="如果不提供, 则直接用bam进行call变异, 省略掉把ubam和bam合并的步骤. 第一列为样本名,分析结果使用名;"
                         " 第二列为read1的fastq路径,如有多个,需要分号分隔; "
                         "第三列为可选, read2的fastq路径,单端测序时则无第三列")
args = parser.parse_args()
args.script_path = script_path

if args.only_show_steps:
    if args.fastq_info is None:
        test_data_dir = os.path.join(os.path.dirname(args.script_path), 'testdata')
        if os.path.exists(test_data_dir):
            args.fastq_info = os.path.join(test_data_dir, 'fastq_info.txt')
            args.compare = os.path.join(test_data_dir, 'compare')
            args.group = os.path.join(test_data_dir, 'group')

nc = NestedCmd(args)
if (args.pipeline_cfg is None) and (not args.bam):
    raise Exception('-bam or -pipeline_cfg is needed! Use -h for help')

with open("cmd." + str(time.time()) + ".txt", 'w') as f:
    f.write(' '.join(sys.argv) + '\n')
    f.write('Detail: \n')
    for k, v in args.__dict__.items():
        f.write('{}: {}\n'.format(k, v))


def pipeline():
    """
    注意：
    * 为了能正常跳过一些步骤,步骤名即step_name不能包含'_'
    * step_name不能重复, 保证最后生成的步骤名不能有重复**
    """
    # variant_calling
    faked_index_bam_cmds = dict()
    with open(nc.workflow_arguments.bam) as fr:
        for line in fr:
            bam = line.strip()
            sample_name = os.path.basename(bam.split('.', 1)[0])
            faked_index_bam_cmds[f'step_{sample_name}'] = nc.cmd_dict(
                cmd='faked_cmd',
                sample_name=sample_name,
                sorted_bam=bam,
                out=bam,
            )
    if nc.workflow_arguments.fastq_info:
        # fastq -> ubam
        fastq_info_dict = nc.parse_fastq_info(nc.workflow_arguments.fastq_info)
        fastq2sam_cmds = nc.RawFastqToSam_cmds(fastq_info_dict)
        merge_bam_cmds = nc.MergeBamAlignment_cmds(faked_index_bam_cmds, fastq2sam_cmds)
        # 修改对faked_index_bam_cmd的依赖
        for each_step in nc.workflow:
            if each_step.startswith('MergeBamAlignment_'):
                depend = nc.workflow[each_step]['depend'].split(',')[1]
                nc.workflow[each_step]['depend'] = depend
        markdup_cmds = nc.MarkDuplicates_cmds(merge_bam_cmds)
    else:
        merge_bam_cmds = faked_index_bam_cmds
        markdup_cmds = nc.MarkDuplicates_cmds(merge_bam_cmds)
        # 修改对faked_index_bam_cmd的依赖为空
        for each_step in nc.workflow:
            if each_step.startswith('MarkDuplicates_'):
                nc.workflow[each_step]['depend'] = ''

    split_cigar_cmds = nc.SplitNCigarReads_cmds(markdup_cmds)
    recal_cmds = nc.BaseRecalibrator_cmds(split_cigar_cmds)
    apply_cmds = nc.ApplyBQSR_cmds(split_cigar_cmds, recal_cmds)
    call_var_cmds = nc.HaplotypeCaller_cmds(apply_cmds)
    filter_vcf_cmds = nc.VariantFiltration_cmds(call_var_cmds)
    nc.run()


if __name__ == '__main__':
    pipeline()
