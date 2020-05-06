import sys
import os
sys.path.append("..")

import rnaseq_report.one_cmd_generator as cmdx
from basic.workflow_basic import Basic


class NestedCmd(Basic):
    def __init__(self, workflow_arguments):
        super().__init__(workflow_arguments)
        terminate = self.do_some_pre_judge(workflow_arguments)
        if terminate:
            exit(0)
        # ----------------------------------------------------
        self.script_dir = os.path.dirname(self.workflow_arguments.script_path)
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

    def gene_body_coverage_cmd(self, files:list, step_name='1.GeneBodyCoverage', parent_dir='3.QCFigures'):
        commands = dict()
        self.mkdir(os.path.join(self.project_dir, parent_dir))
        outdir = os.path.join(self.project_dir, parent_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['gene_body_coverage'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        if 'gene_body_coverage' in self.arg_pool['image_description']:
            args['desc'] = self.arg_pool['image_description']['gene_body_coverage']
        cmd = cmdx.gene_body_coverage(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def inner_distance_cmd(self, files:list, step_name='2.InnerDistance', parent_dir='3.QCFigures'):
        commands = dict()
        outdir = os.path.join(self.project_dir, parent_dir, step_name)
        os.makedirs(outdir, exist_ok=True)
        args = dict(self.arg_pool['inner_distance'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        if 'inner_distance' in self.arg_pool['image_description']:
            args['desc'] = self.arg_pool['image_description']['inner_distance']
        cmd = cmdx.inner_distance(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def read_distribution_cmd(self, files:list, step_name='3.ReadDistribution', parent_dir='3.QCFigures'):
        commands = dict()
        outdir = os.path.join(self.project_dir, parent_dir, step_name)
        os.makedirs(outdir, exist_ok=True)
        args = dict(self.arg_pool['read_distribution'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        if 'read_distribution' in self.arg_pool['image_description']:
            args['desc'] = self.arg_pool['image_description']['read_distribution']
        cmd = cmdx.read_distribution(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def read_duplication_cmd(self, files:list, step_name='4.ReadDuplication', parent_dir='3.QCFigures'):
        commands = dict()
        outdir = os.path.join(self.project_dir, parent_dir, step_name)
        os.makedirs(outdir, exist_ok=True)
        args = dict(self.arg_pool['read_duplication'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        if 'read_duplication' in self.arg_pool['image_description']:
            args['desc'] = self.arg_pool['image_description']['read_duplication']
        cmd = cmdx.read_duplication(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def chromosome_read_distribution_cmd(self, files:list, step_name='5.ChrReadDistribution', parent_dir='3.QCFigures'):
        commands = dict()
        outdir = os.path.join(self.project_dir, parent_dir, step_name)
        os.makedirs(outdir, exist_ok=True)
        args = dict(self.arg_pool['chromosome_read_distribution'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        if 'chromosome_read_distribution' in self.arg_pool['image_description']:
            args['desc'] = self.arg_pool['image_description']['chromosome_read_distribution']
        cmd = cmdx.chromosome_read_distribution(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def exp_saturation_cmd(self, files:list, step_name='6.ExpSaturation', parent_dir='3.QCFigures'):
        commands = dict()
        outdir = os.path.join(self.project_dir, parent_dir, step_name)
        os.makedirs(outdir, exist_ok=True)
        args = dict(self.arg_pool['exp_saturation'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        if 'exp_saturation' in self.arg_pool['image_description']:
            args['desc'] = self.arg_pool['image_description']['exp_saturation']
        cmd = cmdx.exp_saturation(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def exp_pca_cmd(self, exp_table, step_name='1.ExpPCA', parent_dir='4.ExpAnalysis'):
        commands = dict()
        outdir = os.path.join(self.project_dir, parent_dir, step_name)
        os.makedirs(outdir, exist_ok=True)
        args = dict(self.arg_pool['exp_pca'])
        args['exp_table'] = exp_table
        args['group_dict'] = self.workflow_arguments.group_dict
        args['outdir'] = outdir
        if 'exp_pca' in self.arg_pool['image_description']:
            args['desc'] = self.arg_pool['image_description']['exp_pca']
        cmd = cmdx.exp_pca(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def exp_density_cmd(self, exp_table, step_name='2.ExpDensity', parent_dir='4.ExpAnalysis'):
        commands = dict()
        outdir = os.path.join(self.project_dir, parent_dir, step_name)
        os.makedirs(outdir, exist_ok=True)
        args = dict(self.arg_pool['exp_density'])
        args['exp_table'] = exp_table
        # args['group_dict'] = self.workflow_arguments.group_dict
        args['outdir'] = outdir
        if 'exp_density' in self.arg_pool['image_description']:
            args['desc'] = self.arg_pool['image_description']['exp_density']
        cmd = cmdx.exp_density(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def sample_correlation_cmd(self, exp_table, step_name='3.SampleCorrelation', parent_dir='4.ExpAnalysis'):
        commands = dict()
        outdir = os.path.join(self.project_dir, parent_dir, step_name)
        os.makedirs(outdir, exist_ok=True)
        args = dict(self.arg_pool['sample_correlation'])
        args['data_file'] = exp_table
        args['sample_group'] = self.workflow_arguments.group_dict
        args['out_prefix'] = os.path.join(outdir, 'CorrelationCluster')
        if 'sample_correlation' in self.arg_pool['image_description']:
            args['desc'] = self.arg_pool['image_description']['sample_correlation']
        cmd = cmdx.sample_correlation(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def merge_qc_metrics_cmd(self, project_outdir, step_name='2.QCTable', parent_dir='2.QCSummary'):
        commands = dict()
        outdir = os.path.join(self.project_dir, parent_dir, '1.SampleInfo')
        os.makedirs(outdir, exist_ok=True)
        outdir = os.path.join(self.project_dir, parent_dir, step_name)
        os.makedirs(outdir, exist_ok=True)
        args = dict(self.arg_pool['merge_qc_metrics'])
        args['project_outdir'] = project_outdir
        args['outdir'] = outdir
        cmd = cmdx.merge_qc_metrics(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands
