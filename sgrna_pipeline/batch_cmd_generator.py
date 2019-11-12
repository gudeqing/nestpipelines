import os
import sys
sys.path.append("..")

import sgrna_pipeline.one_cmd_generator as cmdx
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

    def cutadapt_cmds(self, fastq_info, step_name='cutadpat'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)

        import pandas as pd
        df = pd.read_csv(fastq_info, header=None, index_col=0, sep='\t')
        for sample, cols in df.iterrows():
            cols = list(cols)
            args = dict(self.arg_pool['cutadapt'])
            args['input_read1'] = cols[0]
            args['output_read1'] = os.path.join(outdir,  sample + '.clean.R1.fq.gz')

            if len(cols) >= 2:
                args['input_read2'] = cols[1]
                args['output_read2'] = os.path.join(outdir, sample + '.clean.R2.fq.gz')

            cmd = cmdx.cutadapt(**args)
            commands[step_name+'_'+sample] = self.cmd_dict(
                cmd=cmd,
                cpu=3,
                mem=1024 ** 3 * 1,
                sample=sample,
                output_read1=args['output_read1'],
            )
            if len(cols) > 2:
                commands[step_name+'_'+sample]['output_read2'] = args['output_read2']
        self.workflow.update(commands)
        return commands

    def fastqc_cmds(self, depend_cmds: dict, step_name='fastqc'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)

        for step, cmd_info in depend_cmds.items():
            sample = cmd_info['sample']
            args = dict(self.arg_pool['fastqc'])
            args['outdir'] = outdir
            args['tmpdir'] = outdir
            args['fastqs'] = cmd_info['output_read1']
            if 'output_read2' in cmd_info:
                args['fastqs'] += ' {}'.format(cmd_info['output_read2'])
            cmd = cmdx.fastqc(**args)
            commands[step_name+'_'+sample] = self.cmd_dict(
                cmd=cmd,
                cpu=5,
                mem=1024 ** 3 * 1,
                sample=sample,
                depend=step
            )
        self.workflow.update(commands)
        return commands

    def mageck_count(self, depend_cmds:dict, step_name='MageckCount'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        prefix = os.path.join(outdir, 'mageck')
        samples = list()
        fastq_read1 = list()
        fastq2_read2 = list()
        steps = list()
        for step, cmd_info in depend_cmds.items():
            steps.append(step)
            samples.append(cmd_info['sample'])
            fastq_read1.append(cmd_info['output_read1'])
            if 'output_read2' in cmd_info:
                fastq2_read2.append(cmd_info['output_read2'])
        args = dict(self.arg_pool['mageck_count'])
        args['out_prefix'] = prefix
        args['fastq_read1'] = ' '.join(fastq_read1)
        if fastq2_read2:
            args['fastq_read2'] = ' '.join(fastq2_read2)
        args['sample_label'] = ','.join(samples)
        cmd = cmdx.mageck_count(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            cpu=2,
            mem=1024 ** 3 * 1,
            depend=','.join(steps),
            out_prefix = args['out_prefix']
        )
        self.workflow.update(commands)
        return commands

    def mageck_count_with_rawfastq(self, fastq_info, step_name='MageckRawCount'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        prefix = os.path.join(outdir, 'mageck')
        samples = list()
        fastq_read1 = list()
        fastq_read2 = list()
        import pandas as pd
        df = pd.read_csv(fastq_info, header=None, index_col=0, sep='\t')
        for sample, cols in df.iterrows():
            samples.append(sample)
            cols = list(cols)
            fastq_read1.append(cols[0])
            if len(cols) >= 2:
                fastq_read2.append(cols[1])
        args = dict(self.arg_pool['mageck_count'])
        args['out_prefix'] = prefix
        args['fastq_read1'] = ' '.join(fastq_read1)
        if fastq_read2:
            args['fastq_read2'] = ' '.join(fastq_read2)
        args['sample_label'] = ','.join(samples)
        cmd = cmdx.mageck_count(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            cpu=2,
            mem=1024 ** 3 * 1,
            out_prefix=args['out_prefix']
        )
        self.workflow.update(commands)
        return commands

    def bowtie2_index(self, step_name='bowtie2Index'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args = dict(self.arg_pool['bowtie2_index'])
        fasta = args['fasta']
        args['basename'] = fasta.rsplit('.')[0] + '.index'
        cmd = cmdx.bowtie2_index(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            cpu=5,
            mem=1024 ** 3 * 1,
            index=args['basename']
        )
        self.workflow.update(commands)
        return commands

    def bowtie2_align(self, depend_cmds, depend_index_cmd, step_name='bowtie2Align'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for step, cmd_info in depend_index_cmd.items():
            index_step = step
            index = cmd_info['index']
        for step, cmd_info in depend_cmds.items():
            sample = cmd_info['sample']
            args = dict(self.arg_pool['bowtie2_align'])
            args['index'] = index
            args['out_bam'] = os.path.join(outdir, '{}.bam'.format(sample))
            args['fastq_read1'] = cmd_info['output_read1']
            if 'output_read2' in cmd_info:
                args['fastq_read2'] = ' {}'.format(cmd_info['output_read2'])
            cmd = cmdx.bowtie2_align(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                cpu=5,
                mem=1024 ** 3 * 1,
                sample=sample,
                depend=step + ',' + index_step,
                out_bam=args['out_bam']
            )
        self.workflow.update(commands)
        return commands

    def mageck_count_with_bam(self, depend_cmds:dict, step_name='MageckBamCount'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        prefix = os.path.join(outdir, 'mageck')
        samples = list()
        bams = list()
        steps = list()
        for step, cmd_info in depend_cmds.items():
            steps.append(step)
            samples.append(cmd_info['sample'])
            bams.append(cmd_info['out_bam'])
        args = dict(self.arg_pool['mageck_count'])
        args['out_prefix'] = prefix
        args['sample_label'] = ','.join(samples)
        args['fastq_read1'] = ' '.join(bams)
        cmd = cmdx.mageck_count(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            cpu=2,
            mem=1024 ** 3 * 1,
            depend=','.join(steps),
            out_prefix=args['out_prefix']
        )
        self.workflow.update(commands)
        return commands

