#! /data/users/dqgu/anaconda3/bin/python
# coding=utf-8
import sys
import os
import time

script_path = os.path.abspath(__file__)
if os.path.islink(script_path):
    script_path = os.readlink(script_path)
sys.path.append(os.path.dirname(os.path.dirname(script_path)))

from xxx_pipeline.batch_cmd_generator import NestedCmd
from basic.workflow_basic import basic_arg_parser

# 初始化参数
if len(sys.argv) <= 1:
    exit('please provide at least one argument, use -h for help')
parser = basic_arg_parser()
# 可以增加新的流程参数
# parser.add_argument('-new_arg', required=False)

# 收集参数和记录命令行信息
args = parser.parse_args()
args.script_path = script_path
with open("cmd." + str(time.time()) + ".txt", 'w') as f:
    f.write(' '.join(sys.argv) + '\n')
    f.write('Detail: \n')
    for k, v in args.__dict__.items():
        f.write('{}: {}\n'.format(k, v))

# 可以在这里增加对参数的判断
if args.pipeline_cfg is not None:
    print('')


# 从这里开始写pipeline, 有时一个步骤有两种选择, 请在这里自行判断

def pipeline():
    """
    注意：
    * 为了能正常跳过一些步骤,步骤名即step_name不能包含'_'
    * step_name不能重复, 保证最后生成的步骤名不能有重复**
    """
    nc = NestedCmd(args)
    nc.which_cmds()
    cmds = nc.which_cmds()
    nc.which2_cmds(depend_cmds=cmds)
    nc.run()


if __name__ == '__main__':
    pipeline()
