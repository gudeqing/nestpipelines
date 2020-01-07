import os
import sys
sys.path.append("..")
from basic.workflow_basic import Basic
import xxx_pipeline.one_cmd_generator as cmdx


class NestedCmd(Basic):
    def __init__(self, workflow_arguments):
        super().__init__(workflow_arguments)
        terminate = self.do_some_pre_judge(workflow_arguments)
        if terminate:
            exit(0)
        # self.arg_pool is from Basic after processing workflow_arguments
        # self.cmd_dict is from Basic after processing workflow_arguments
        # self.project_dir is from Basic after processing workflow_arguments
        # self.workflow is from Basic after processing workflow_arguments

    def show_cmd_example(self, cmd_name):
        """
        :param cmd_name: cmd_generator中的函数名,也是arguments.ini中的section名
        :return: None
        """
        if not self.workflow_arguments.arg_cfg:
            raise Exception("please first input arg_cfg ")
        if cmd_name not in self.arg_pool:
            raise Exception('please provide valid cmd_name, refer --list_cmd_names')
        exec("print(cmdx.{}(**self.arg_pool['{}']))".format(cmd_name, cmd_name))

    # 从下开始修改, 每个函数用来批量生成一组cmds, 并且囊括了其依赖的cmd信息
    def which_cmds(self, step_name='1.BasicStat', main_step_name='1.QC'):
        commands = dict()
        outdir = os.path.join(self.project_dir, main_step_name, step_name)
        os.makedirs(outdir, exist_ok=True)
        # 1. get args from pool
        args = dict(self.arg_pool['which'])
        # 2. update args
        args['prefix'] = os.path.join(outdir, 'all')
        # 3. format cmds
        cmd = cmdx.which(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            cpu=2,
            mem=1024 ** 3 * 1,
            metadata=args['metadata'],
            out_prefix=args['prefix'],
        )
        # 4. update workflow
        self.workflow.update(commands)
        return commands

    def which2_cmds(self, depend_cmds, step_name='3.VJ-pair-bar3d', main_step_name='3.VJ-Usage'):
        commands = dict()
        outdir = os.path.join(self.project_dir, main_step_name, step_name)
        os.makedirs(outdir, exist_ok=True)
        for step, cmd_info in depend_cmds.items():
            sample = cmd_info['sample']
            args = dict(self.arg_pool['Plot3dVJUsage'])
            args['data'] = cmd_info['out_prefix'] + '.fancyvj.wt.txt'
            args['out'] = os.path.join(outdir, sample + '.VJ.3dBar.png')
            cmd = cmdx.which2(**args)
            commands[step_name+'_'+sample] = self.cmd_dict(
                depend=step,
                cmd=cmd,
                cpu=2,
                mem=1024 ** 3 * 1,
                sample=sample,
                out=args['out']
            )
        self.workflow.update(commands)
        return commands
