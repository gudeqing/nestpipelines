# coding=utf-8
import sys
import os
sys.path.append('..')

from sentieon_rnaseq.batch_cmd_generator import NestedCmd
from basic.workflow_basic import basic_arg_parser

parser = basic_arg_parser()
parser.add_argument('-b', nargs='+', help="bam文件路径, 空格分隔", required=False)
args = parser.parse_args()
args.script_path = os.path.abspath(__file__)

nc = NestedCmd(args)
if (args.pipeline_cfg is None) and (not args.b) and(not args.e):
    raise Exception('-b or -be or -pipeline_cfg is needed! Use -h for help')


def pipeline():
    score_cmds = nc.score_cmds(args.b)
    dedup_cmds = nc.dedup_cmds(score_cmds)
    split_cmds = nc.split_cmds(dedup_cmds)
    realign_cmds = nc.realign_cmds(split_cmds)
    recalibrate_cmds = nc.recalibrate_cmds(realign_cmds)
    recalibrate2_cmds = nc.recalibrate2_cmds(recalibrate_cmds)
    recalibrate3_cmds = nc.recalibrate3_cmds(recalibrate2_cmds)
    recalibrate4_cmds = nc.recalibrate4_cmds(recalibrate3_cmds)
    calling_cmds = nc.calling_cmds(recalibrate2_cmds)
    nc.run()


if __name__ == '__main__':
    pipeline()
