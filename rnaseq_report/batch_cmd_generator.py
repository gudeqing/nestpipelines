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
        slider_dir = os.path.join(self.project_dir, 'htmls')
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'GeneBodyCoverage.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = "基因覆盖度, 正常为钟罩型, 纵坐标为覆盖度百分比，即横坐标区域对应的覆盖度除以整个基因区域的覆盖度的上四分位值"
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_fragment_size_slider(self, depend_cmd:dict, step_name='FragmentSizeSlider'):
        commands = dict()
        slider_dir = os.path.join(self.project_dir, 'htmls')
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'FragmentSize.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = "Pair-read 长度分布图"
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_inner_size_slider(self, depend_cmd:dict, step_name='InnerSizeSlider'):
        commands = dict()
        slider_dir = os.path.join(self.project_dir, 'htmls')
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'InnerSize.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = "pair-end read 之间的相对距离，负数表明有测序有overlap"
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_read_distribution_slider(self, depend_cmd:dict, step_name='ReadDistributionSlider'):
        commands = dict()
        slider_dir = os.path.join(self.project_dir, 'htmls')
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'ReadDistribution.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = "read比对到各个基因区域的比例分布图"
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_read_duplication_slider(self, depend_cmd:dict, step_name='ReadDuplicationSlider'):
        commands = dict()
        slider_dir = os.path.join(self.project_dir, 'htmls')
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'ReadDuplication.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = "冗余度分析图，正常时左边快速下降，尾部平滑无突起"
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_exp_saturation_slider(self, depend_cmd:dict, step_name='ExpSaturationSlider'):
        commands = dict()
        slider_dir = os.path.join(self.project_dir, 'htmls')
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'ExpSaturation.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = "基因测序饱和度分析，曲线(中位线)越快到达最低点表示测序越饱和",
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_chr_read_distribution_slider(self, depend_cmd:dict, step_name='ChrReadDistributionSlider'):
        commands = dict()
        slider_dir = os.path.join(self.project_dir, 'htmls')
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'ChrReadDistribution.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = "read比对到各个染色体或scaffold的统计分布图"
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_exp_distribution_slider(self, depend_cmd:dict, step_name='ExpDistributionSlider'):
        commands = dict()
        slider_dir = os.path.join(self.project_dir, 'htmls')
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'ExpressionDistribution.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = "表达量分布密度图"
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_corr_cluster_slider(self, depend_cmd:dict, step_name='CorrelationClusterSlider'):
        commands = dict()
        slider_dir = os.path.join(self.project_dir, 'htmls')
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'CorrelationCluster.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = "基于相关性作为距离度量的样本聚类分析，用于查看样本重复性和离群样本",
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_pca_slider(self, depend_cmd:dict, step_name='PCASlider'):
        commands = dict()
        slider_dir = os.path.join(self.project_dir, 'htmls')
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'PCA.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = "主成份分析，可用于查看样本重复性和离群样本"
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_AlignmentSummary_slider(self, depend_cmd:dict, step_name='AlignmentSummarySlider'):
        commands = dict()
        slider_dir = os.path.join(self.project_dir, 'htmls')
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'CollectAlignmentSummaryMetrics.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = "比对结果统计表"
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_InsertSizeSlider_slider(self, depend_cmd:dict, step_name='InsertSizeSlider'):
        commands = dict()
        slider_dir = os.path.join(self.project_dir, 'htmls')
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'CollectInsertSizeMetrics.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = "Insert Size (包含内含子的长度) 统计表"
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_RnaSeqMetrics_slider(self, depend_cmd:dict, step_name='RnaSeqMetricsSlider'):
        commands = dict()
        slider_dir = os.path.join(self.project_dir, 'htmls')
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'CollectRnaSeqMetrics.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = "基于RNA尺度的比对结果统计表"
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands

    def make_TargetedSummary_slider(self, depend_cmd:dict, step_name='TargetedSummarySlider'):
        commands = dict()
        slider_dir = os.path.join(self.project_dir, 'htmls')
        args = dict(self.arg_pool['make_slider'])
        for step, cmd_info in depend_cmd.items():
            args['out'] = os.path.join(slider_dir, 'TargetedSummaryMetrics.html')
            args['images'] = os.path.join(cmd_info['outdir'], '*.html')
            args['image_ids'] = None
            args['image_desc'] = "比对到捕获区域的统计表"
            cmd = cmdx.make_slider(**args)
            commands[step_name] = self.cmd_dict(depend=step, cmd=cmd, outdir=slider_dir)
        self.workflow.update(commands)
        return commands


