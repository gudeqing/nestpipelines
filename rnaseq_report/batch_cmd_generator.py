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
        self.slider_dir = os.path.join(self.project_dir, 'Summary')
        self.script_dir = os.path.dirname(self.workflow_arguments.script_path)
        if self.arg_pool:
            self.arg_pool['make_slider']['template'] = os.path.join(self.script_dir, 'templates', 'slide.jinja2')
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

    def gene_body_coverage_cmd(self, files:list, step_name='geneBodyCoverage'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['gene_body_coverage'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        cmd = cmdx.gene_body_coverage(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def fragment_length_cmd(self, files:list, step_name='fragmentLength'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['fragment_length'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        cmd = cmdx.fragment_length(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def inner_distance_cmd(self, files:list, step_name='innerDistance'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['inner_distance'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        cmd = cmdx.inner_distance(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def read_distribution_cmd(self, files:list, step_name='readDistribution'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['read_distribution'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        cmd = cmdx.read_distribution(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def read_duplication_cmd(self, files:list, step_name='readDuplication'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['read_duplication'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        cmd = cmdx.read_duplication(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def chromosome_read_distribution_cmd(self, files:list, step_name='chrReadDistribution'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['chromosome_read_distribution'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        cmd = cmdx.chromosome_read_distribution(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def exp_saturation_cmd(self, files:list, step_name='expSaturation'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['exp_saturation'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        cmd = cmdx.exp_saturation(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def exp_pca_cmd(self, exp_table, step_name='expPCA'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['exp_pca'])
        args['exp_table'] = exp_table
        args['group_dict'] = self.workflow_arguments.group_dict
        args['outdir'] = outdir
        cmd = cmdx.exp_pca(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def exp_density_cmd(self, exp_table, step_name='expDensity'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['exp_density'])
        args['exp_table'] = exp_table
        # args['group_dict'] = self.workflow_arguments.group_dict
        args['outdir'] = outdir
        cmd = cmdx.exp_density(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def sample_correlation_cmd(self, exp_table, step_name='sampleCorrelation'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['sample_correlation'])
        args['data_file'] = exp_table
        args['sample_group'] = self.workflow_arguments.group_dict
        args['out_name'] = os.path.join(outdir, 'CorrelationCluster.html')
        cmd = cmdx.sample_correlation(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def CollectAlignmentSummaryMetrics(self, files:list, step_name='CollectAlignmentSummaryMetrics'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['CollectAlignmentSummaryMetrics'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        cmd = cmdx.CollectAlignmentSummaryMetrics(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def CollectInsertSizeMetrics(self, files:list, step_name='CollectInsertSizeMetrics'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['CollectInsertSizeMetrics'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        cmd = cmdx.CollectInsertSizeMetrics(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def CollectRnaSeqMetrics(self, files:list, step_name='CollectRnaSeqMetrics'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['CollectRnaSeqMetrics'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        cmd = cmdx.CollectRnaSeqMetrics(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def CollectTargetedPcrMetrics(self, files:list, step_name='CollectTargetedPcrMetrics'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['CollectTargetedPcrMetrics'])
        args['files'] = ' '.join(files)
        args['outdir'] = outdir
        cmd = cmdx.CollectTargetedPcrMetrics(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            outdir=outdir,
        )
        self.workflow.update(commands)
        return commands

    def make_gene_body_cov_slider(self, depend_cmd:dict, step_name='GeneBodyCovSlider'):
        commands = dict()
        slider_dir = self.slider_dir
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'GeneBodyCoverage.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = self.arg_pool['image_description']['gene_body_coverage']
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_fragment_size_slider(self, depend_cmd:dict, step_name='FragmentSizeSlider'):
        commands = dict()
        slider_dir = self.slider_dir
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'FragmentSize.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = self.arg_pool['image_description']['fragment_length']
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_inner_distance_slider(self, depend_cmd:dict, step_name='InnerDistanceSlider'):
        commands = dict()
        slider_dir = self.slider_dir
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'InnerDistance.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = self.arg_pool['image_description']['inner_distance']
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_read_distribution_slider(self, depend_cmd:dict, step_name='ReadDistributionSlider'):
        commands = dict()
        slider_dir = self.slider_dir
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'ReadDistribution.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = self.arg_pool['image_description']['read_distribution']
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_read_duplication_slider(self, depend_cmd:dict, step_name='ReadDuplicationSlider'):
        commands = dict()
        slider_dir = self.slider_dir
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'ReadDuplication.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = self.arg_pool['image_description']['read_duplication']
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_exp_saturation_slider(self, depend_cmd:dict, step_name='ExpSaturationSlider'):
        commands = dict()
        slider_dir = self.slider_dir
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'ExpSaturation.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = self.arg_pool['image_description']['exp_saturation']
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_chr_read_distribution_slider(self, depend_cmd:dict, step_name='ChrReadDistributionSlider'):
        commands = dict()
        slider_dir = self.slider_dir
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'ChrReadDistribution.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = self.arg_pool['image_description']['chromosome_read_distribution']
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_exp_distribution_slider(self, depend_cmd:dict, step_name='ExpDistributionSlider'):
        commands = dict()
        slider_dir = self.slider_dir
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'ExpressionDistribution.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = self.arg_pool['image_description']['exp_density']
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_corr_cluster_slider(self, depend_cmd:dict, step_name='CorrelationClusterSlider'):
        commands = dict()
        slider_dir = self.slider_dir
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'CorrelationCluster.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = self.arg_pool['image_description']['sample_correlation']
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_pca_slider(self, depend_cmd:dict, step_name='PCASlider'):
        commands = dict()
        slider_dir = self.slider_dir
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'PCA.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = self.arg_pool['image_description']['exp_pca']
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_AlignmentSummary_slider(self, depend_cmd:dict, step_name='AlignmentSummarySlider'):
        commands = dict()
        slider_dir = self.slider_dir
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'CollectAlignmentSummaryMetrics.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = self.arg_pool['image_description']['CollectAlignmentSummaryMetrics']
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_InsertSizeSlider_slider(self, depend_cmd:dict, step_name='InsertSizeSlider'):
        commands = dict()
        slider_dir = self.slider_dir
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'CollectInsertSizeMetrics.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = self.arg_pool['image_description']['CollectInsertSizeMetrics']
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_RnaSeqMetrics_slider(self, depend_cmd:dict, step_name='RnaSeqMetricsSlider'):
        commands = dict()
        slider_dir = self.slider_dir
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'CollectRnaSeqMetrics.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = self.arg_pool['image_description']['CollectRnaSeqMetrics']
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_TargetedSummary_slider(self, depend_cmd:dict, step_name='TargetedSummarySlider'):
        commands = dict()
        slider_dir = self.slider_dir
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'TargetedSummaryMetrics.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = self.arg_pool['image_description']['CollectTargetedPcrMetrics']
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands


