import os
import sys
sys.path.append("..")

import sentieon_rnaseq.one_cmd_generator as cmdx
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

    def score_cmds(self, bam_list, step_name='score'):
        if type(bam_list) == str:
            import glob
            bam_list = glob.glob(bam_list)
        if not bam_list:
            raise Exception('No bam found!')
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for bam in bam_list:
            sample = os.path.basename(bam).split('.', 1)[0]
            args = dict(self.arg_pool['score'])
            args['bam'] = bam
            args['score_txt'] = os.path.join(outdir, '{}.score.txt'.format(sample))
            cmd = cmdx.score(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                cpu=3,
                mem=1024**3*10,
                score_txt=args['score_txt'],
                bam=args['bam'],
                sample=sample,
            )
        self.workflow.update(commands)
        return commands

    def dedup_cmds(self, depend_cmds: dict, step_name='dedup'):
        if not depend_cmds:
            return
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for step, cmd_info in depend_cmds.items():
            sample = cmd_info['sample']
            args = dict(self.arg_pool['dedup'])
            args['bam'] = cmd_info['bam']
            args['score_txt'] = cmd_info['score_txt']
            args['dedup_metric_txt'] = os.path.join(outdir, '{}.dedup_metric.txt'.format(sample))
            args['dedup_bam'] = os.path.join(outdir, '{}.dedup.bam'.format(sample))
            cmd = cmdx.dedup(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                mem=1024 ** 3 * 10,
                cpu=5,
                depend=step,
                sample=sample,
                dedup_bam=args['dedup_bam']
            )
        self.workflow.update(commands)
        return commands

    def split_cmds(self, depend_cmds: dict, step_name='split'):
        if not depend_cmds:
            return
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for step, cmd_info in depend_cmds.items():
            sample = cmd_info['sample']
            args = dict(self.arg_pool['split'])
            args['bam'] = cmd_info['dedup_bam']
            args['split_bam'] = os.path.join(outdir, '{}.split.bam'.format(sample))
            cmd = cmdx.split(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                mem=1024 ** 3 * 10,
                cpu=5,
                depend=step,
                sample=sample,
                split_bam=args['split_bam']
            )
        self.workflow.update(commands)
        return commands

    def realign_cmds(self, depend_cmds: dict, step_name='realign'):
        if not depend_cmds:
            return
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for step, cmd_info in depend_cmds.items():
            sample = cmd_info['sample']
            args = dict(self.arg_pool['realign'])
            args['bam'] = cmd_info['split_bam']
            args['realigned_bam'] = os.path.join(outdir, '{}.realigned.bam'.format(sample))
            cmd = cmdx.realign(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                mem=1024 ** 3 * 10,
                cpu=5,
                depend=step,
                sample=sample,
                realigned_bam=args['realigned_bam']
            )
        self.workflow.update(commands)
        return commands

    def recalibrate_cmds(self, depend_cmds: dict, step_name='recalibrate'):
        if not depend_cmds:
            return
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for step, cmd_info in depend_cmds.items():
            sample = cmd_info['sample']
            args = dict(self.arg_pool['recalibrate'])
            args['bam'] = cmd_info['realigned_bam']
            args['recal_data_table'] = os.path.join(outdir, '{}.recal_data.table'.format(sample))
            cmd = cmdx.recalibrate(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                mem=1024 ** 3 * 10,
                cpu=5,
                depend=step,
                sample=sample,
                realigned_bam=cmd_info['realigned_bam'],
                recal_data_table=args['recal_data_table']
            )
        self.workflow.update(commands)
        return commands

    def calling_cmds(self, depend_cmds: dict, step_name='calling'):
        if not depend_cmds:
            return
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for step, cmd_info in depend_cmds.items():
            sample = cmd_info['sample']
            args = dict(self.arg_pool['calling'])
            args['bam'] = cmd_info['realigned_bam']
            args['recal_data_table'] = cmd_info['recal_data_table']
            args['variant_vcf'] = os.path.join(outdir, '{}.variant.vcf'.format(sample))
            cmd = cmdx.calling(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                mem=1024 ** 3 * 12,
                cpu=6,
                depend=step,
                sample=sample,
                variant_vcf=args['variant_vcf']
            )
        self.workflow.update(commands)
        return commands
