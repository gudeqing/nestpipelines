#! /data/users/dqgu/anaconda3/bin/python
# coding=utf-8
import sys
import os

script_path = os.path.abspath(__file__)
if os.path.islink(script_path):
    script_path = os.readlink(script_path)
sys.path.append(os.path.dirname(os.path.dirname(script_path)))

from sentieon_rnaseq.batch_cmd_generator import NestedCmd
from basic.workflow_basic import basic_arg_parser

parser = basic_arg_parser()
parser.add_argument('-b', nargs='+', help="bam文件路径, 空格分隔; 也可以是一个存储bam路径的文件", required=False)
args = parser.parse_args()
if args.b:
    if len(args.b) == 1 and os.path.exists(args.b[0]):
        args.b = [x.strip() for x in open(args.b[0])]
args.script_path = os.path.abspath(__file__)

nc = NestedCmd(args)
if (args.pipeline_cfg is None) and (not args.b):
    raise Exception('-b or -pipeline_cfg is needed! Use -h for help')


def pipeline():
    score_cmds = nc.score_cmds(args.b)
    dedup_cmds = nc.dedup_cmds(score_cmds)
    split_cmds = nc.split_cmds(dedup_cmds)
    realign_cmds = nc.realign_cmds(split_cmds)
    recalibrate_cmds = nc.recalibrate_cmds(realign_cmds)
    calling_cmds = nc.calling_cmds(recalibrate_cmds)
    nc.run()


if __name__ == '__main__':
    pipeline()
