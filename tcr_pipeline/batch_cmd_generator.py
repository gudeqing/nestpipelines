import os
import sys
sys.path.append("..")

import tcr_pipeline.one_cmd_generator as cmdx
from basic.workflow_basic import Basic


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

    def CalcBasicStats_cmds(self, step_name='1.BasicStats', main_step_name='1.BasicStats'):
        commands = dict()
        outdir = os.path.join(self.project_dir, main_step_name)
        self.mkdir(outdir)
        outdir = os.path.join(outdir, step_name)
        self.mkdir(outdir)
        # get args from pool
        args = dict(self.arg_pool['CalcBasicStats'])
        # update args
        args['prefix'] = os.path.join(outdir, 'all')
        # get cmd
        cmd = cmdx.CalcBasicStats(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            cpu=2,
            mem=1024 ** 3 * 1,
            metadata=args['metadata'],
            out_prefix=args['prefix'],
        )
        self.workflow.update(commands)
        return commands

    def CalcSegmentUsage_cmds(self, step_name='2.SegmentUsage', main_step_name='1.BasicStats'):
        commands = dict()
        outdir = os.path.join(self.project_dir, main_step_name)
        self.mkdir(outdir)
        outdir = os.path.join(outdir, step_name)
        self.mkdir(outdir)
        # get args from pool
        args = dict(self.arg_pool['CalcSegmentUsage'])
        # update args
        args['prefix'] = os.path.join(outdir, 'all')
        # get cmd
        cmd = cmdx.CalcSegmentUsage(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            cpu=2,
            mem=1024 ** 3 * 1,
            metadata=args['metadata'],
            out_prefix=args['prefix'],
        )
        self.workflow.update(commands)
        return commands

    def CalcSpectratype_cmds(self, step_name='3.Spectratype', main_step_name='1.BasicStats'):
        commands = dict()
        outdir = os.path.join(self.project_dir, main_step_name)
        self.mkdir(outdir)
        outdir = os.path.join(outdir, step_name)
        self.mkdir(outdir)
        # get args from pool
        args = dict(self.arg_pool['CalcSpectratype'])
        # update args
        args['prefix'] = os.path.join(outdir, 'all')
        # get cmd
        cmd = cmdx.CalcSpectratype(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            cpu=2,
            mem=1024 ** 3 * 1,
            metadata=args['metadata'],
            out_prefix=args['prefix'],
        )
        self.workflow.update(commands)
        return commands

    def PlotFancySpectratype_cmds(self, step_name='4.PlotSpectratype', main_step_name='1.BasicStats'):
        commands = dict()
        outdir = os.path.join(self.project_dir, main_step_name)
        self.mkdir(outdir)
        outdir = os.path.join(outdir, step_name)
        self.mkdir(outdir)
        import pandas as pd
        input_args = dict(self.arg_pool['vdjtools'])
        meta_df = pd.read_csv(input_args['metadata'], sep='\t')
        for each_file, sample in zip(meta_df['files'], meta_df[input_args['label_field']]):
            # get args from pool
            args = dict(self.arg_pool['PlotFancySpectratype'])
            # update args
            args['prefix'] = os.path.join(outdir, sample)
            args['sample'] = each_file
            # get cmd
            cmd = cmdx.PlotFancySpectratype(**args)
            commands[step_name+'_'+sample] = self.cmd_dict(
                cmd=cmd,
                cpu=2,
                mem=1024 ** 3 * 1,
                out_prefix=args['prefix'],
            )
        self.workflow.update(commands)
        return commands

    def PlotFancyVJUsage_cmds(self, step_name='5.VJUsage', main_step_name='1.BasicStats'):
        commands = dict()
        outdir = os.path.join(self.project_dir, main_step_name)
        self.mkdir(outdir)
        outdir = os.path.join(outdir, step_name)
        self.mkdir(outdir)
        import pandas as pd
        input_args = dict(self.arg_pool['vdjtools'])
        meta_df = pd.read_csv(input_args['metadata'], sep='\t')
        for each_file, sample in zip(meta_df['files'], meta_df[input_args['label_field']]):
            # get args from pool
            args = dict(self.arg_pool['PlotFancyVJUsage'])
            # update args
            args['prefix'] = os.path.join(outdir, sample)
            args['sample'] = each_file
            # get cmd
            cmd = cmdx.PlotFancyVJUsage(**args)
            commands[step_name+'_'+sample] = self.cmd_dict(
                cmd=cmd,
                cpu=2,
                mem=1024 ** 3 * 1,
                out_prefix=args['prefix'],
            )
        self.workflow.update(commands)
        return commands

    def PlotSpectratypeV_cmds(self, step_name='6.PlotSpectratypeV', main_step_name='1.BasicStats'):
        commands = dict()
        outdir = os.path.join(self.project_dir, main_step_name)
        self.mkdir(outdir)
        outdir = os.path.join(outdir, step_name)
        self.mkdir(outdir)
        import pandas as pd
        input_args = dict(self.arg_pool['vdjtools'])
        meta_df = pd.read_csv(input_args['metadata'], sep='\t')
        for each_file, sample in zip(meta_df['files'], meta_df[input_args['label_field']]):
            # get args from pool
            args = dict(self.arg_pool['PlotSpectratypeV'])
            # update args
            args['prefix'] = os.path.join(outdir, sample)
            args['sample'] = each_file
            # get cmd
            cmd = cmdx.PlotSpectratypeV(**args)
            commands[step_name+'_'+sample] = self.cmd_dict(
                cmd=cmd,
                cpu=2,
                mem=1024 ** 3 * 1,
                out_prefix=args['prefix'],
            )
        self.workflow.update(commands)
        return commands

    def PlotQuantileStats_cmds(self, step_name='1.PlotQuantileStats', main_step_name='2.Diversity'):
        commands = dict()
        outdir = os.path.join(self.project_dir, main_step_name)
        self.mkdir(outdir)
        outdir = os.path.join(outdir, step_name)
        self.mkdir(outdir)
        import pandas as pd
        input_args = dict(self.arg_pool['vdjtools'])
        meta_df = pd.read_csv(input_args['metadata'], sep='\t')
        for each_file, sample in zip(meta_df['files'], meta_df[input_args['label_field']]):
            # get args from pool
            args = dict(self.arg_pool['PlotQuantileStats'])
            # update args
            args['prefix'] = os.path.join(outdir, sample)
            args['sample'] = each_file
            # get cmd
            cmd = cmdx.PlotQuantileStats(**args)
            commands[step_name+'_'+sample] = self.cmd_dict(
                cmd=cmd,
                cpu=2,
                mem=1024 ** 3 * 1,
                out_prefix=args['prefix'],
            )
        self.workflow.update(commands)
        return commands

    def RarefactionPlot_cmds(self, step_name='2.RarefactionPlot', main_step_name='2.Diversity'):
        commands = dict()
        outdir = os.path.join(self.project_dir, main_step_name)
        self.mkdir(outdir)
        outdir = os.path.join(outdir, step_name)
        self.mkdir(outdir)
        # get args from pool
        args = dict(self.arg_pool['RarefactionPlot'])
        # update args
        args['prefix'] = os.path.join(outdir, 'all')
        # get cmd
        cmd = cmdx.RarefactionPlot(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            cpu=2,
            mem=1024 ** 3 * 1,
            metadata=args['metadata'],
            out_prefix=args['prefix'],
        )
        self.workflow.update(commands)
        return commands

    def CalcDiversityStats_cmds(self, step_name='3.DiversityStats', main_step_name='2.Diversity'):
        commands = dict()
        outdir = os.path.join(self.project_dir, main_step_name)
        self.mkdir(outdir)
        outdir = os.path.join(outdir, step_name)
        self.mkdir(outdir)
        # get args from pool
        args = dict(self.arg_pool['CalcDiversityStats'])
        # update args
        args['prefix'] = os.path.join(outdir, 'all')
        # get cmd
        cmd = cmdx.CalcDiversityStats(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            cpu=2,
            mem=1024 ** 3 * 1,
            metadata=args['metadata'],
            out_prefix=args['prefix'],
        )
        self.workflow.update(commands)
        return commands

    def CalcPairwiseDistances_cmds(self, step_name='1.PairwiseDistances', main_step_name='3.Cluster'):
        commands = dict()
        outdir = os.path.join(self.project_dir, main_step_name)
        self.mkdir(outdir)
        outdir = os.path.join(outdir, step_name)
        self.mkdir(outdir)
        # get args from pool
        args = dict(self.arg_pool['CalcPairwiseDistances'])
        # update args
        args['prefix'] = os.path.join(outdir, 'all')
        # get cmd
        cmd = cmdx.CalcPairwiseDistances(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            cpu=2,
            mem=1024 ** 3 * 1,
            metadata=args['metadata'],
            out_prefix=args['prefix'],
        )
        self.workflow.update(commands)
        return commands

    def ClusterSamples_cmds(self, depend_cmd, step_name='2.ClusterSamples', main_step_name='3.Cluster'):
        commands = dict()
        outdir = os.path.join(self.project_dir, main_step_name)
        self.mkdir(outdir)
        outdir = os.path.join(outdir, step_name)
        self.mkdir(outdir)
        for step, cmd_info in depend_cmd.items():
            # get args from pool
            args = dict(self.arg_pool['ClusterSamples'])
            # update args
            args['out_prefix'] = os.path.join(outdir, 'all')
            args['input_prefix'] = cmd_info['out_prefix']
            # get cmd
            cmd = cmdx.ClusterSamples(**args)
            commands[step_name] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                cpu=2,
                mem=1024 ** 3 * 1,
                out_prefix=args['out_prefix'],
            )
        self.workflow.update(commands)
        return commands
