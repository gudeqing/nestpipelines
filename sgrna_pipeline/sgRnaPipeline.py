#! /data/users/dqgu/anaconda3/bin/python
# coding=utf-8
import sys
import os

script_path = os.path.abspath(__file__)
if os.path.islink(script_path):
    script_path = os.readlink(script_path)
sys.path.append(os.path.dirname(os.path.dirname(script_path)))

from sgrna_pipeline.batch_cmd_generator import NestedCmd
from basic.workflow_basic import basic_arg_parser

parser = basic_arg_parser()
parser.add_argument('-fastq_info', required=False,
                    help="第一列为样本名,分析结果使用名; 第二列为read1的fastq路径,如有多个,需要分号分隔; "
                         "第三列为可选, read2的fastq路径,单端测序时则无第三列")

args = parser.parse_args()
args.script_path = os.path.abspath(__file__)

nc = NestedCmd(args)
if (args.pipeline_cfg is None) and (not args.fastq_info):
    raise Exception('-fastq_info or -pipeline_cfg is needed! Use -h for help')


def pipeline():
    cutadapt_cmds = nc.cutadapt_cmds(args.fastq_info)
    fastqc_cmds = nc.fastqc_cmds(cutadapt_cmds)
    mageck_cmds = nc.mageck_count(cutadapt_cmds)
    nc.run()


if __name__ == '__main__':
    pipeline()
