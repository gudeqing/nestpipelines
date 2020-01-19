#! /data/users/dqgu/anaconda3/bin/python
# coding=utf-8
import sys
import os
import time

script_path = os.path.abspath(__file__)
if os.path.islink(script_path):
    script_path = os.readlink(script_path)
sys.path.append(os.path.dirname(os.path.dirname(script_path)))

from tcr_pipeline.batch_cmd_generator import NestedCmd
from basic.workflow_basic import basic_arg_parser

parser = basic_arg_parser()
args = parser.parse_args()
args.script_path = os.path.abspath(__file__)


def pipeline():
    nc = NestedCmd(args)
    nc.CalcBasicStats_cmds()
    nc.CalcSegmentUsage_cmds()
    nc.CalcSpectratype_cmds()
    nc.PlotFancySpectratype_cmds()
    cmds = nc.PlotFancyVJUsage_cmds()
    nc.Plot3dVJUsage_cmds(depend_cmds=cmds)
    nc.PlotCDRVenn_cmds()
    nc.PlotSpectratypeV_cmds()
    nc.PlotQuantileStats_cmds()
    nc.RarefactionPlot_cmds()
    nc.CalcDiversityStats_cmds()
    distance_cmd = nc.CalcPairwiseDistances_cmds()
    nc.ClusterSamples_cmds(depend_cmd=distance_cmd)
    nc.run()


if __name__ == '__main__':
    pipeline()
