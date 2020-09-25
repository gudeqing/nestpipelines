import os
import sys
sys.path.append('..')

import rnaseq_pipeline.cmd_generator as cmdx
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

    def fastqc_raw_data_cmds(self, fastq_info_dict, step_name='RawDataQC'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for sample, fq_list in fastq_info_dict.items():
            fastqs = ' '.join(x for y in fq_list for x in y)
            args = dict(self.arg_pool['fastqc'])
            args['outdir'] = outdir
            args['tmpdir'] = outdir
            args['fastqs'] = fastqs
            cmd = cmdx.fastqc(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(cmd=cmd)
        self.workflow.update(commands)
        return commands

    def fastqc_trimmed_data_cmds(self, trimming_cmds, step_name='TrimmedDataQC'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        samples = set(trimming_cmds[x]['sample_name'] for x in trimming_cmds)
        for sample in samples:
            depend_steps = list()
            trimmed_fq_list = list()
            for step, cmd_info in trimming_cmds.items():
                if cmd_info['sample_name'] == sample:
                    depend_steps.append(step)
                    trimmed_fq_list.append(cmd_info['trimmed_fq1'])
                    if 'trimmed_fq2' in cmd_info:
                        trimmed_fq_list.append(cmd_info['trimmed_fq2'])
            args = dict(self.arg_pool['fastqc'])
            args['outdir'] = outdir
            args['tmpdir'] = outdir
            args['fastqs'] = ' '.join(trimmed_fq_list)
            cmd = cmdx.fastqc(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd, depend=','.join(depend_steps), result_dir=outdir, sample_name=sample,
            )
        self.workflow.update(commands)
        return commands

    def trimmomatic_cmds(self, fastq_info_dict, step_name='Trim'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for sample, fq_list in fastq_info_dict.items():
            true_sample = sample
            mode = 'PE' if len(fq_list) == 2 else "SE"
            if mode == 'PE':
                fq1_list = fq_list[0]
                fq2_list = fq_list[1]
                if len(fq1_list) != len(fq2_list):
                    raise Exception('fastq文件的read1和read2必须一一对应')
                for ind, (fq1, fq2) in enumerate(zip(fq1_list, fq2_list)):
                    if ind >= 1:
                        sample = sample + '_' + str(ind)
                    if 'MapSplice' in self.workflow_arguments.skip:
                        trimmed_fq1 = os.path.join(outdir, sample + '_clean_R1.fq.gz')
                        trimmed_fq2 = os.path.join(outdir, sample + '_clean_R2.fq.gz')
                    else:
                        # MapSplice需要非压缩的fastq作为输入
                        trimmed_fq1 = os.path.join(outdir, sample + '_clean_R1.fq')
                        trimmed_fq2 = os.path.join(outdir, sample + '_clean_R2.fq')
                    unpaired_fq1 = os.path.join(outdir, sample + '_unpaired_clean_R1.fq.gz')
                    unpaired_fq2 = os.path.join(outdir, sample + '_unpaired_clean_R2.fq.gz')
                    args = dict(self.arg_pool['trimmomatic'])
                    args.update(dict(
                        fq1=fq1, fq2=fq2, mode=mode,
                        trimmed_fq1=trimmed_fq1, trimmed_fq2=trimmed_fq2,
                        unpaired_fq1=unpaired_fq1, unpaired_fq2=unpaired_fq2
                    ))
                    cmd = cmdx.trimmomatic(**args)
                    commands[step_name + '_' + sample] = self.cmd_dict(
                        cmd=cmd, sample_name=true_sample,
                        trimmed_fq1=trimmed_fq1, trimmed_fq2=trimmed_fq2
                    )
            else:
                for ind, fq1 in enumerate(fq_list[0]):
                    if ind >= 1:
                        sample = sample + '_' + str(ind)
                    trimmed_fq1 = os.path.join(outdir, sample + '_clean_R1.fq.gz')
                    args = dict(self.arg_pool['trimmomatic'])
                    args.update(dict(fq1=fq1, trimmed_fq1=trimmed_fq1, mode=mode))
                    cmd = cmdx.trimmomatic(**args)
                    commands[step_name + '_' + sample] = self.cmd_dict(
                        cmd=cmd, sample_name=true_sample, trimmed_fq1=trimmed_fq1
                    )
        self.workflow.update(commands)
        return commands

    def star_index_cmd(self, step_name='AlignIndex'):
        commands = dict()
        if os.path.exists(self.arg_pool['star_index']['genomeDir']):
            if len(os.listdir(self.arg_pool['star_index']['genomeDir'])) <= 1:
                self.mkdir(self.arg_pool['star_index']['genomeDir'])
            else:
                self.logger.warning('STAR index existed, and skip this indexing step!')
                return commands
        else:
            self.mkdir(self.arg_pool['star_index']['genomeDir'])
        cmd = cmdx.star_index(**self.arg_pool['star_index'])
        commands[step_name] = self.cmd_dict(
            cmd=cmd, cpu=1, mem=2 * 1024 ** 3, retry=1,
            monitor_time_step=5)
        self.workflow.update(commands)
        return commands

    def star_align_cmds(self, trimming_cmds, index_cmd, step_name='Align', chimeric_in_bam=False):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        samples = set(trimming_cmds[x]['sample_name'] for x in trimming_cmds)
        for sample in samples:
            depend_steps = list() if not index_cmd else list(index_cmd.keys())
            trimmed_fq1_list = list()
            trimmed_fq2_list = list()
            for step, cmd_info in trimming_cmds.items():
                if cmd_info['sample_name'] == sample:
                    depend_steps.append(step)
                    trimmed_fq1_list.append(cmd_info['trimmed_fq1'])
                    if 'trimmed_fq2' in cmd_info:
                        trimmed_fq2_list.append(cmd_info['trimmed_fq2'])
            args = dict(self.arg_pool['star_align'])
            args['readFilesIn'] = ','.join(trimmed_fq1_list) + ' ' + ','.join(trimmed_fq2_list)
            result_dir = os.path.join(outdir, sample)
            self.mkdir(result_dir)
            prefix = os.path.join(result_dir, sample + '.')
            args['outFileNamePrefix'] = prefix
            args['outSAMattrRGline'] = 'ID:{prefix} SM:{prefix} PL:{platform}'.format(prefix=sample, platform="illumina")
            if chimeric_in_bam:
                args['chimOutType'] = 'WithinBAM'
                args['chimMultimapNmax'] = 0
            cmd = cmdx.star_align(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd, mem=1024 ** 3 * 10, cpu=2,
                monitor_time_step=5, depend=','.join(depend_steps),
                sorted_bam='{}Aligned.sortedByCoord.out.bam'.format(prefix),
                Chimeric_out_junction='{}Chimeric.out.junction'.format(prefix),
                sample_name=sample,
                result_dir=result_dir,
                readFilesIn=args['readFilesIn']
            )
        self.workflow.update(commands)
        return commands

    def star_align_with_rawdata_cmds(self, fastq_info_dict, index_cmd, step_name='Align', chimeric_in_bam=False):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for sample, fq_list in fastq_info_dict.items():
            depend_steps = list() if not index_cmd else list(index_cmd.keys())
            mode = 'PE' if len(fq_list) == 2 else "SE"
            result_dir = os.path.join(outdir, sample)
            self.mkdir(result_dir)
            args = dict(self.arg_pool['star_align'])
            if mode in 'PE':
                args['readFilesIn'] = ','.join(fq_list[0]) + ' ' + ','.join(fq_list[1])
            else:
                args['readFilesIn'] = ','.join(fq_list[0])
            result_dir = os.path.join(outdir, sample)
            self.mkdir(result_dir)
            prefix = os.path.join(result_dir, sample + '.')
            args['outFileNamePrefix'] = prefix
            args['outSAMattrRGline'] = 'ID:{prefix} SM:{prefix} PL:{platform}'.format(prefix=sample, platform="illumina")
            if chimeric_in_bam:
                args['chimOutType'] = 'WithinBAM'
                args['chimMultimapNmax'] = 0
            cmd = cmdx.star_align(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd, mem=1024 ** 3 * 10, cpu=2,
                monitor_time_step=5, depend=','.join(depend_steps),
                sorted_bam='{}Aligned.sortedByCoord.out.bam'.format(prefix),
                Chimeric_out_junction='{}Chimeric.out.junction'.format(prefix),
                sample_name=sample,
                result_dir=result_dir,
                readFilesIn=args['readFilesIn']
            )
        self.workflow.update(commands)
        return commands

    def bam_index_cmds(self, align_cmds, step_name='IndexBam'):
        commands = dict()
        for step, cmd_info in align_cmds.items():
            sample = cmd_info['sample_name']
            args = dict(self.arg_pool['samtools_index'])
            args['bam'] = cmd_info['sorted_bam']
            cmd = cmdx.samtools_index(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd, mem=1024 ** 3 * 2, cpu=2, monitor_time_step=5,
                depend=step,
                sample_name=sample,
                sorted_bam=args['bam']
            )
        self.workflow.update(commands)
        return commands

    def star_fusion_cmds(self, align_cmds, step_name='StarFusion'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for step, cmd_info in align_cmds.items():
            sample = cmd_info['sample_name']
            result_dir = os.path.join(outdir, sample)
            self.mkdir(result_dir)
            args = dict(self.arg_pool['star_fusion'])
            args['outdir'] = result_dir
            args['Chimeric_out_junction'] = cmd_info['Chimeric_out_junction']
            args['left_fq'] = cmd_info['readFilesIn'].strip().split()[0]
            if len(cmd_info['readFilesIn'].strip().split()) >= 2:
                args['right_fq'] = cmd_info['readFilesIn'].strip().split()[1]
            cmd = cmdx.star_fusion(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd, mem=1024 ** 3 * 2, cpu=5, monitor_time_step=5,
                depend=step,
                sample_name=sample,
                outdir=args['outdir']
            )
        self.workflow.update(commands)
        return commands

    def arriba_cmds(self, align_cmds, step_name='ArribaFusion'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for step, cmd_info in align_cmds.items():
            sample = cmd_info['sample_name']
            result_dir = os.path.join(outdir, sample)
            self.mkdir(result_dir)
            args = dict(self.arg_pool['arriba'])
            args['o'] = f'{result_dir}/{sample}.fusions.tsv'
            args['discarded'] = f'{result_dir}/{sample}.discarded.fusions.tsv'
            args['x'] = cmd_info['sorted_bam']
            cmd = cmdx.arriba(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd, mem=1024 ** 3 * 5, cpu=5, monitor_time_step=5,
                depend=step,
                sample_name=sample,
                outdir=outdir,
                fusions=args['o']
            )
        self.workflow.update(commands)
        return commands

    def scallop_cmds(self, align_cmds, step_name='Assembly'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for step, cmd_info in align_cmds.items():
            sample = cmd_info['sample_name']
            args = dict(self.arg_pool['scallop'])
            args['bam'] = cmd_info['sorted_bam']
            result_dir = os.path.join(outdir, sample)
            self.mkdir(result_dir)
            args['out_gtf'] = os.path.join(result_dir, '{}.scallop.gtf'.format(sample))
            cmd = cmdx.scallop(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd, mem=1024 ** 3 * 2, cpu=2, monitor_time_step=5,
                depend=step,
                sample_name=sample,
                out_gtf=args['out_gtf']
            )
        self.workflow.update(commands)
        return commands

    def merge_scallop_transcripts_cmd(self, assemble_cmds, step_name='MergeTranscript'):
        gtfs = list()
        for step, cmd_info in assemble_cmds.items():
            gtfs.append(cmd_info['out_gtf'])
        args = dict(self.arg_pool['merge_scallop_transcripts'])
        args['gtf'] = ' '.join(gtfs)
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        args['outdir'] = outdir
        cmd = cmdx.merge_scallop_transcripts(**args)
        commands = dict()
        commands[step_name] = self.cmd_dict(
            cmd=cmd, cpu=2, depend=','.join(assemble_cmds.keys()),
            result_dir=outdir, all_transcripts=os.path.join(outdir, 'all.transcripts.fa'),
            all_gtf=os.path.join(outdir, 'all.transcripts.gtf'),
            novel_transcripts=os.path.join(outdir, 'novel.transcripts.fa')
        )
        self.workflow.update(commands)
        return commands

    def salmon_index_cmd(self, step_name='QuantIndex', merge_transcript_cmd=None):
        commands = dict()
        if not merge_transcript_cmd:
            if os.path.exists(self.arg_pool['salmon_index']['index_prefix']):
                self.logger.warning('salmon index existed, and skip this indexing step!')
                self.workflow.update(commands)
                return commands
            args = dict(self.arg_pool['salmon_index'])
            cmd = cmdx.salmon_index(**args)
            commands[step_name] = self.cmd_dict(cmd=cmd)
        else:
            outdir = os.path.join(self.project_dir, step_name)
            self.mkdir(outdir)
            args = dict(self.arg_pool['salmon_index'])
            merge_step = list(merge_transcript_cmd.keys())[0]
            args['transcript_fasta'] = merge_transcript_cmd[merge_step]['all_transcripts']
            args['index_prefix'] = outdir
            cmd = cmdx.salmon_index(**args)
            commands[step_name] = self.cmd_dict(
                cmd=cmd, index_prefix=outdir, depend=merge_step,
                all_gtf=merge_transcript_cmd[merge_step]['all_gtf']
            )
        self.workflow.update(commands)
        return commands

    def salmon_quant_with_raw_data_cmds(self, fastq_info_dict, index_cmd, step_name='Quant'):
        commands = dict()
        result_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(result_dir)
        for sample, fq_list in fastq_info_dict.items():
            depend_steps = list() if not index_cmd else list(index_cmd.keys())
            mode = 'PE' if len(fq_list) == 2 else "SE"
            if mode == 'PE':
                fq1_list = fq_list[0]
                fq2_list = fq_list[1]
                if len(fq1_list) != len(fq2_list):
                    raise Exception('fastq文件的read1和read2必须一一对应')
                args = dict(self.arg_pool['salmon_quant'])
                if index_cmd:
                    if 'all_gtf' in list(index_cmd.values())[0]:
                        args['transcript2gene'] = list(index_cmd.values())[0]['all_gtf']
                        args['index'] = list(index_cmd.values())[0]['index_prefix']
                args['fq1'] = ' '.join(fq1_list)
                args['fq2'] = ' '.join(fq2_list)
                args['mode'] = mode
                prefix = os.path.join(result_dir, sample)
                args['out_prefix'] = prefix
                cmd = cmdx.salmon_quant(**args)
                commands[step_name + '_' + sample] = self.cmd_dict(
                    cmd=cmd, mem=1024 ** 3 * 2, cpu=2, monitor_time_step=5,
                    depend=','.join(depend_steps),
                    sample_name=sample,
                    result_dir=result_dir,
                    prefix=prefix
                )
            else:
                args = dict(self.arg_pool['salmon_quant'])
                if index_cmd:
                    if 'all_gtf' in index_cmd.values()[0]:
                        args['transcript2gene'] = list(index_cmd.values())[0]['all_gtf']
                        args['index'] = list(index_cmd.values())[0]['index_prefix']
                args['fq'] = ' '.join(fq_list[0])
                args['mode'] = mode
                prefix = os.path.join(result_dir, sample)
                args['out_prefix'] = prefix
                cmd = cmdx.salmon_quant(**args)
                commands[step_name + '_' + sample] = self.cmd_dict(
                    cmd=cmd, mem=1024 ** 3 * 2, cpu=2, monitor_time_step=5,
                    depend=','.join(depend_steps),
                    sample_name=sample,
                    result_dir=result_dir,
                    prefix=prefix
                )
        self.workflow.update(commands)
        return commands

    def salmon_quant_with_clean_data_cmds(self, trimming_cmds, index_cmd, step_name='Quant'):
        commands = dict()
        result_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(result_dir)
        samples = set(trimming_cmds[x]['sample_name'] for x in trimming_cmds)
        for sample in samples:
            depend_steps = list() if not index_cmd else list(index_cmd.keys())
            trimmed_fq1_list = list()
            trimmed_fq2_list = list()
            for step, cmd_info in trimming_cmds.items():
                if cmd_info['sample_name'] == sample:
                    depend_steps.append(step)
                    trimmed_fq1_list.append(cmd_info['trimmed_fq1'])
                    if 'trimmed_fq2' in cmd_info:
                        trimmed_fq2_list.append(cmd_info['trimmed_fq2'])
            mode = 'PE' if trimmed_fq2_list else 'SE'
            if mode == 'PE':
                if len(trimmed_fq1_list) != len(trimmed_fq2_list):
                    raise Exception('fastq文件的read1和read2必须一一对应')
                args = dict(self.arg_pool['salmon_quant'])
                if index_cmd:
                    if 'all_gtf' in list(index_cmd.values())[0]:
                        args['transcript2gene'] = list(index_cmd.values())[0]['all_gtf']
                        args['index'] = list(index_cmd.values())[0]['index_prefix']
                args['fq1'] = ' '.join(trimmed_fq1_list)
                args['fq2'] = ' '.join(trimmed_fq2_list)
                args['mode'] = mode
                prefix = os.path.join(result_dir, sample)
                args['out_prefix'] = prefix
                cmd = cmdx.salmon_quant(**args)
                commands[step_name + '_' + sample] = self.cmd_dict(
                    cmd=cmd, mem=1024 ** 3 * 2, cpu=2, monitor_time_step=5,
                    depend=','.join(depend_steps),
                    sample_name=sample,
                    result_dir=result_dir,
                    prefix=prefix
                )
            else:
                args = dict(self.arg_pool['salmon_quant'])
                if index_cmd:
                    if 'all_gtf' in list(index_cmd.values())[0]:
                        args['transcript2gene'] = list(index_cmd.values())[0]['all_gtf']
                        args['index'] = list(index_cmd.values())[0]['index_prefix']
                args['fq'] = ' '.join(trimmed_fq1_list)
                args['mode'] = mode
                prefix = os.path.join(result_dir, sample)
                args['out_prefix'] = prefix
                cmd = cmdx.salmon_quant(**args)
                commands[step_name + '_' + sample] = self.cmd_dict(
                    cmd=cmd, mem=1024 ** 3 * 2, cpu=2, monitor_time_step=5,
                    depend=','.join(depend_steps),
                    sample_name=sample,
                    result_dir=result_dir,
                    prefix=prefix
                )
        self.workflow.update(commands)
        return commands

    def merge_quant_cmd(self, quant_cmds, step_name='MergeQuant', quant_method='salmon', level='gene'):
        commands = dict()
        depend = list()
        step_name = step_name + level.capitalize()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        result_dir = ''
        for step, cmd_info in quant_cmds.items():
            result_dir = cmd_info['result_dir']
            depend.append(step)
        args = dict(self.arg_pool['abundance_estimates_to_matrix'])
        args['est_method'] = quant_method
        args['out_prefix'] = os.path.join(out_dir, level)
        if not level == 'gene':
            args['quant_result'] = result_dir + '/*/quant.sf'
            cmd = cmdx.abundance_estimates_to_matrix(**args)
            commands[step_name] = self.cmd_dict(
                cmd=cmd, monitor_resource=False, check_resource_before_run=False,
                depend=','.join(depend), out_prefix=args['out_prefix'],
                transcript_tpm_matrix=args['out_prefix'] + '.isoform.TMM.EXPR.matrix',
                transcript_count_matrix=args['out_prefix'] + '.isoform.counts.matrix',
            )
        else:
            args['quant_result'] = result_dir + '/*/quant.genes.sf'
            cmd = cmdx.abundance_estimates_to_matrix(**args)
            commands[step_name] = self.cmd_dict(
                cmd=cmd, monitor_resource=False, check_resource_before_run=False,
                depend=','.join(depend), out_prefix=args['out_prefix'],
                gene_tpm_matrix=args['out_prefix'] + '.isoform.TMM.EXPR.matrix',
                gene_count_matrix=args['out_prefix'] + '.isoform.counts.matrix',
            )
        self.workflow.update(commands)
        return commands

    def diff_exp_cmd(self, merge_cmd, step_name='Diff', level='gene'):
        commands = dict()
        step_name = step_name + level.capitalize()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        # gene diff exp
        if level == 'gene':
            depend = [x for x in merge_cmd.keys() if x.endswith('Gene')][0]
            depend_info = merge_cmd[depend]
            args = dict(self.arg_pool['diff_exp'])
            args['result_dir'] = out_dir
            args['count_matrix'] = depend_info['gene_count_matrix']
            args['exp_matrix'] = depend_info['gene_tpm_matrix']
        else:
            depend = [x for x in merge_cmd.keys() if x.endswith('Transcript')][0]
            depend_info = merge_cmd[depend]
            args = dict(self.arg_pool['diff_exp'])
            args['result_dir'] = out_dir
            args['count_matrix'] = depend_info['transcript_count_matrix']
            args['exp_matrix'] = depend_info['transcript_tpm_matrix']
        args['group_info'] = self.workflow_arguments.group
        args['comparison_info'] = self.workflow_arguments.compare
        cmd = cmdx.diff_exp(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            cpu=int(args['threads']) + 1,
            depend=depend,
            result_dir=args['result_dir']
        )
        self.workflow.update(commands)
        return commands

    def go_enrich_cmds(self, diffexp_cmd, step_name='GoEnrich', level='gene'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name + level.capitalize())
        self.mkdir(out_dir)
        depend = list(diffexp_cmd.keys())[0]
        depend_info = list(diffexp_cmd.values())[0]
        diff_list_files = list()
        with open(self.workflow_arguments .compare) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                ctrl, test = line.strip().split()
                diff_list_files.append(os.path.join(depend_info['result_dir'], '{}_vs_{}.*.DE.list'.format(ctrl, test)))
        for each in diff_list_files:
            args = dict(self.arg_pool['goatools'])
            args['study'] = each
            cmp_name = os.path.basename(each).split('.', 1)[0]
            args['goea_out'] = os.path.join(out_dir, str(cmp_name) + '.goea.xls')
            args['dag_out'] = os.path.join(out_dir, str(cmp_name) + '.dag.svg')
            if level == 'gene':
                cmd = cmdx.goatools(**args)
            else:
                args['population'] = args['trans_population']
                args['gene2go'] = args['trans_association']
                cmd = cmdx.goatools(**args)
            commands[step_name + level.capitalize() + '_' + str(cmp_name)] = self.cmd_dict(
                cmd=cmd,
                cpu=1,
                depend=depend,
                result_dir=out_dir
            )
        self.workflow.update(commands)
        return commands

    def kegg_enrich_cmds(self, diffexp_cmd, step_name='KeggEnrich', level='gene'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name + level.capitalize())
        self.mkdir(out_dir)
        depend = list(diffexp_cmd.keys())[0]
        depend_info = list(diffexp_cmd.values())[0]
        diff_list_files = list()
        with open(self.workflow_arguments .compare) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                ctrl, test = line.strip().split()
                diff_list_files.append(os.path.join(depend_info['result_dir'], '{}_vs_{}.*.DE.list'.format(ctrl, test)))
        for each in diff_list_files:
            args = dict(self.arg_pool['kegg_enrich'])
            args['deg'] = each
            args['outdir'] = out_dir
            if level == 'gene':
                cmd = cmdx.kegg_enrich(**args)
            else:
                args['g2k'] = args['t2k']
                args['g2p'] = args['t2p']
                cmd = cmdx.kegg_enrich(**args)
            cmp_name = os.path.basename(each).split('.', 1)[0]
            commands[step_name + level.capitalize() + '_' + str(cmp_name)] = self.cmd_dict(
                cmd=cmd,
                cpu=1,
                depend=depend,
                result_dir=out_dir
            )
        self.workflow.update(commands)
        return commands

    def exp_analysis_cmd(self, merge_cmd, step_name='ExpAnalysis', level='gene'):
        commands = dict()
        step_name = step_name + level.capitalize()
        out_dir = os.path.join(self.project_dir, step_name)
        depend = list(merge_cmd.keys())[0]
        depend_info = merge_cmd[depend]
        if ',' not in depend_info['depend']:
            self.logger.warning(f'Skip {level} expression analysis for there is only one sample!')
            return commands
        self.mkdir(out_dir)
        args = dict(self.arg_pool['exp_analysis'])
        if level == 'gene':
            args['out_prefix'] = os.path.join(out_dir, 'gene.tpm')
            args['matrix'] = depend_info['gene_tpm_matrix']
        else:
            args['out_prefix'] = os.path.join(out_dir, 'transcript.tpm')
            args['matrix'] = depend_info['transcript_tpm_matrix']
        cmd = cmdx.exp_analysis(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            cpu=1,
            depend=depend,
            result_dir=out_dir
        )
        self.workflow.update(commands)
        return commands

    def gene_body_coverage_cmds(self, index_bam_cmds, step_name='GeneBodyCoverage'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['gene_body_coverage'])
        for step, cmd_info in index_bam_cmds.items():
            sample = cmd_info['sample_name']
            args['bam'] = cmd_info['sorted_bam']
            args['out_prefix'] = os.path.join(out_dir, sample)
            cmd = cmdx.gene_body_coverage(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                out_prefix=args['out_prefix']
            )
        self.workflow.update(commands)
        return commands

    def inner_distance_cmds(self, index_bam_cmds, step_name='InnerDistance'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['inner_distance'])
        for step, cmd_info in index_bam_cmds.items():
            sample = cmd_info['sample_name']
            args['bam'] = cmd_info['sorted_bam']
            args['out_prefix'] = os.path.join(out_dir, sample)
            cmd = cmdx.inner_distance(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                out_prefix=args['out_prefix']
            )
        self.workflow.update(commands)
        return commands

    def read_distribution_cmds(self, index_bam_cmds, step_name='ReadDistribution'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['read_distribution'])
        for step, cmd_info in index_bam_cmds.items():
            sample = cmd_info['sample_name']
            args['bam'] = cmd_info['sorted_bam']
            args['outfile'] = os.path.join(out_dir, sample + '.read_distribution.txt')
            cmd = cmdx.read_distribution(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                outfile=args['outfile']
            )
        self.workflow.update(commands)
        return commands

    def chromosome_read_distribution_cmds(self, index_bam_cmds, step_name='ChrReadDistribution'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['chromosome_read_distribution'])
        for step, cmd_info in index_bam_cmds.items():
            sample = cmd_info['sample_name']
            args['bam'] = cmd_info['sorted_bam']
            args['outfile'] = os.path.join(out_dir, '{}.chromosome.alignment.stat.txt'.format(sample))
            cmd = cmdx.chromosome_read_distribution(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                outfile=args['outfile']
            )
        self.workflow.update(commands)
        return commands

    def read_duplication_cmds(self, index_bam_cmds, step_name='ReadDuplication'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['read_duplication'])
        for step, cmd_info in index_bam_cmds.items():
            sample = cmd_info['sample_name']
            args['bam'] = cmd_info['sorted_bam']
            args['out_prefix'] = os.path.join(out_dir, sample)
            cmd = cmdx.read_duplication(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                out_prefix=args['out_prefix']
            )
        self.workflow.update(commands)
        return commands

    def rna_fragment_size_cmds(self, index_bam_cmds, step_name='FragmentSize'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['rna_fragment_size'])
        for step, cmd_info in index_bam_cmds.items():
            sample = cmd_info['sample_name']
            args['bam'] = cmd_info['sorted_bam']
            args['outfile'] = os.path.join(out_dir, sample + '.fragment_size.txt')
            cmd = cmdx.rna_fragment_size(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                outfile=args['outfile']
            )
        self.workflow.update(commands)
        return commands

    def rpkm_saturation_cmds(self, index_bam_cmds, step_name='RPKMSaturation'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['rpkm_saturation'])
        for step, cmd_info in index_bam_cmds.items():
            sample = cmd_info['sample_name']
            args['bam'] = cmd_info['sorted_bam']
            args['out_prefix'] = os.path.join(out_dir, sample)
            cmd = cmdx.rpkm_saturation(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                out_prefix=args['out_prefix']
            )
        self.workflow.update(commands)
        return commands

    def tpm_saturation_cmds(self, index_bam_cmds, step_name='TPMSaturation'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['tpm_saturation'])
        for step, cmd_info in index_bam_cmds.items():
            sample = cmd_info['sample_name']
            args['bam'] = cmd_info['sorted_bam']
            args['outdir'] = out_dir
            cmd = cmdx.tpm_saturation(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                exp_matrix='{}.tpm.xls'.format(sample)
            )
        self.workflow.update(commands)
        return commands

    def get_alignment_summary_cmds(self, index_bam_cmds, step_name='AlignmentSummary'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['get_alignment_summary'])
        for step, cmd_info in index_bam_cmds.items():
            sample = cmd_info['sample_name']
            args['bam'] = cmd_info['sorted_bam']
            args['outdir'] = out_dir
            cmd = cmdx.get_alignment_summary(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                alignment_summary='{}.alignment_summary.json'.format(sample)

            )
        self.workflow.update(commands)
        return commands

    def CollectAlignmentSummaryMetrics_cmds(self, index_bam_cmds, step_name='CollectAlignmentSummaryMetrics'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['CollectAlignmentSummaryMetrics'])
        for step, cmd_info in index_bam_cmds.items():
            sample = cmd_info['sample_name']
            args['bam'] = cmd_info['sorted_bam']
            args['outfile'] = os.path.join(out_dir, '{}.{}.xls'.format(sample, step_name))
            cmd = cmdx.CollectAlignmentSummaryMetrics(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                metrics=args['outfile']

            )
        self.workflow.update(commands)
        return commands

    def CollectInsertSizeMetrics_cmds(self, index_bam_cmds, step_name='CollectInsertSizeMetrics'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['CollectInsertSizeMetrics'])
        for step, cmd_info in index_bam_cmds.items():
            sample = cmd_info['sample_name']
            args['bam'] = cmd_info['sorted_bam']
            args['outfile'] = os.path.join(out_dir, '{}.{}.xls'.format(sample, step_name))
            args['outimage'] = os.path.join(out_dir, '{}.{}.pdf'.format(sample, step_name))
            cmd = cmdx.CollectInsertSizeMetrics(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                metrics=args['outfile']

            )
        self.workflow.update(commands)
        return commands

    def CollectTargetedPcrMetrics_cmds(self, index_bam_cmds, step_name='CollectTargetedPcrMetrics'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['CollectTargetedPcrMetrics'])
        for step, cmd_info in index_bam_cmds.items():
            sample = cmd_info['sample_name']
            args['bam'] = cmd_info['sorted_bam']
            args['outfile'] = os.path.join(out_dir, '{}.{}.xls'.format(sample, step_name))
            args['per_target_coverage_outfile'] = os.path.join(out_dir, '{}.PerTargetCoverage.txt'.format(sample))
            cmd = cmdx.CollectTargetedPcrMetrics(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                metrics=args['outfile']

            )
        self.workflow.update(commands)
        return commands

    def CollectRnaSeqMetrics_cmds(self, index_bam_cmds, step_name='CollectRnaSeqMetrics'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['CollectRnaSeqMetrics'])
        for step, cmd_info in index_bam_cmds.items():
            sample = cmd_info['sample_name']
            args['bam'] = cmd_info['sorted_bam']
            args['outfile'] = os.path.join(out_dir, '{}.{}.xls'.format(sample, step_name))
            args['outimage'] = os.path.join(out_dir, '{}.{}.pdf'.format(sample, step_name))
            cmd = cmdx.CollectRnaSeqMetrics(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                metrics=args['outfile']

            )
        self.workflow.update(commands)
        return commands

    def FastqToSam_cmds(self, trimming_cmds, step_name='FastqToSam'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        samples = set(trimming_cmds[x]['sample_name'] for x in trimming_cmds)
        for sample in samples:
            trimmed_fq1_list = list()
            trimmed_fq2_list = list()
            depends = list()
            for step, cmd_info in trimming_cmds.items():
                if cmd_info['sample_name'] == sample:
                    depends.append(step)
                    trimmed_fq1_list.append(cmd_info['trimmed_fq1'])
                    if 'trimmed_fq2' in cmd_info:
                        trimmed_fq2_list.append(cmd_info['trimmed_fq2'])
            if len(trimmed_fq1_list) > 1 or len(trimmed_fq2_list) > 1:
                raise Exception('目前不支持一个样本的原始数据对应多个fastq这种情况，你可以分析前合并fastq')
            args = dict(self.arg_pool['FastqToSam'])
            args['fq'] = trimmed_fq1_list[0]
            args['fq2'] = trimmed_fq2_list[0]
            args['sample_name'] = sample
            args['out'] = os.path.join(outdir, f'{sample}.ubam')
            cmd = cmdx.FastqToSam(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=depends[0],
                out=args['out'],
                sample_name=sample,
                mem=1024 ** 3 * 2,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def RawFastqToSam_cmds(self, fastq_info_dict, step_name='FastqToSam'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for sample, fq_list in fastq_info_dict.items():
            fq1_list = fq_list[0]
            fq2_list = fq_list[1]
            if len(fq1_list) > 1 or len(fq2_list) > 1:
                raise Exception('目前不支持一个样本的原始数据对应多个fastq这种情况，你可以分析前合并fastq')
            args = dict(self.arg_pool['FastqToSam'])
            args['fq'] = fq1_list[0]
            args['fq2'] = fq2_list[0]
            args['sample_name'] = sample
            args['out'] = os.path.join(outdir, f'{sample}.ubam')
            cmd = cmdx.FastqToSam(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                out=args['out'],
                sample_name=sample,
                mem=1024 ** 3 * 2,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def MergeBamAlignment_cmds(self, index_star_bam_cmds, fastq2bam_cmds, step_name='MergeBamAlignment'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['MergeBamAlignment'])
        for step, cmd_info in index_star_bam_cmds.items():
            sample = cmd_info['sample_name']
            args['ALIGNED'] = cmd_info['sorted_bam']
            another_step = [x for x in fastq2bam_cmds if x.endswith(sample)][0]
            args['UNMAPPED'] = fastq2bam_cmds[another_step]['out']
            args['out'] = os.path.join(out_dir, sample+'.merged.bam')
            cmd = cmdx.MergeBamAlignment(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=','.join([step, another_step]),
                out=args['out'],
                sample_name=sample,
                mem=1024 ** 3 * 2,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def MarkDuplicates_cmds(self, merge_bam_cmds, step_name='MarkDuplicates'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['MarkDuplicates'])
        for step, cmd_info in merge_bam_cmds.items():
            sample = cmd_info['sample_name']
            args['input'] = cmd_info['out']
            args['output'] = os.path.join(out_dir, sample+'.markdup.bam')
            args['metrics'] = os.path.join(out_dir, sample+'.mark_duplicated.txt')
            cmd = cmdx.MarkDuplicates(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                output=args['output'],
                sample_name=sample,
                mem=1024 ** 3 * 3,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def SplitNCigarReads_cmds(self, mark_dup_cmds, step_name='SplitNCigarReads'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['SplitNCigarReads'])
        for step, cmd_info in mark_dup_cmds.items():
            sample = cmd_info['sample_name']
            args['input'] = cmd_info['output']
            args['output'] = os.path.join(out_dir, sample+'.SplitNCigar.bam')
            cmd = cmdx.SplitNCigarReads(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                output=args['output'],
                sample_name=sample,
                mem=1024 ** 3 * 3,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def BaseRecalibrator_cmds(self, split_cmds, step_name='BaseRecalibrator'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['BaseRecalibrator'])
        for step, cmd_info in split_cmds.items():
            sample = cmd_info['sample_name']
            args['input'] = cmd_info['output']
            args['output'] = os.path.join(out_dir, sample+'.recal.table')
            cmd = cmdx.BaseRecalibrator(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                output=args['output'],
                sample_name=sample,
                mem=1024 ** 3 * 3,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def ApplyBQSR_cmds(self, split_cmds, recal_cmds, step_name='ApplyBQSR'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['ApplyBQSR'])
        for step, cmd_info in split_cmds.items():
            sample = cmd_info['sample_name']
            args['input'] = cmd_info['output']
            another_step = [x for x in recal_cmds if x.endswith(sample)][0]
            args['bqsr-recal-file'] = recal_cmds[another_step]['output']
            args['output'] = os.path.join(out_dir, sample+'.ready.bam')
            cmd = cmdx.ApplyBQSR(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=','.join([step, another_step]),
                output=args['output'],
                sample_name=sample,
                mem=1024 ** 3 * 4,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def HaplotypeCaller_cmds(self, bqsr_cmds, step_name='HaplotypeCaller'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['HaplotypeCaller'])
        for step, cmd_info in bqsr_cmds.items():
            sample = cmd_info['sample_name']
            args['input'] = cmd_info['output']
            args['output'] = os.path.join(out_dir, sample+'.raw.vcf.gz')
            cmd = cmdx.HaplotypeCaller(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                output=args['output'],
                sample_name=sample,
                mem=1024 ** 3 * 6,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def VariantFiltration_cmds(self, calling_cmds, step_name='VariantFiltration'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['VariantFiltration'])
        for step, cmd_info in calling_cmds.items():
            sample = cmd_info['sample_name']
            args['vcf'] = cmd_info['output']
            args['out'] = os.path.join(out_dir, sample+'.filtered.vcf.gz')
            cmd = cmdx.VariantFiltration(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                output=args['out'],
                sample_name=sample,
                mem=1024 ** 3 * 3,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def MapSplice_cmds(self, trimming_cmds, step_name='MapSplice'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        samples = set(trimming_cmds[x]['sample_name'] for x in trimming_cmds)
        for sample in samples:
            trimmed_fq1_list = list()
            trimmed_fq2_list = list()
            depends = list()
            for step, cmd_info in trimming_cmds.items():
                if cmd_info['sample_name'] == sample:
                    depends.append(step)
                    trimmed_fq1_list.append(cmd_info['trimmed_fq1'])
                    if 'trimmed_fq2' in cmd_info:
                        trimmed_fq2_list.append(cmd_info['trimmed_fq2'])
            if len(trimmed_fq1_list) > 1 or len(trimmed_fq2_list) > 1:
                raise Exception('目前不支持一个样本的原始数据对应多个fastq这种情况，你可以分析前合并fastq')
            args = dict(self.arg_pool['MapSplice'])
            args['fq'] = trimmed_fq1_list[0]
            args['fq2'] = trimmed_fq2_list[0]
            args['sample_name'] = sample
            args['outdir'] = os.path.join(outdir, sample)
            cmd = cmdx.MapSplice(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=depends[0],
                outdir=args['outdir'],
                sample_name=sample,
                mem=1024 ** 3 * 5,
                cpu=4,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def RawDataMapSplice_cmds(self, fastq_info_dict, step_name='MapSplice'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for sample, fq_list in fastq_info_dict.items():
            fq1_list = fq_list[0]
            fq2_list = fq_list[1]
            if len(fq1_list) > 1 or len(fq2_list) > 1:
                raise Exception('目前不支持一个样本的原始数据对应多个fastq这种情况，你可以分析前合并fastq')
            args = dict(self.arg_pool['MapSplice'])
            args['fq'] = fq1_list[0]
            args['fq2'] = fq2_list[0]
            args['sample_name'] = sample
            args['outdir'] = os.path.join(outdir, sample)
            cmd = cmdx.MapSplice(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                out=args['outdir'],
                sample_name=sample,
                mem=1024 ** 3 * 5,
                cpu=4,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def annovar_cmds(self, depend_cmds, step_name='AnnotateVariation'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['annovar'])
        for step, cmd_info in depend_cmds.items():
            sample = cmd_info['sample_name']
            args['input'] = cmd_info['output']
            args['outfile'] = os.path.join(out_dir, sample)
            cmd = cmdx.annovar(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                output=args['outfile'],
                sample_name=sample,
                mem=1024 ** 3 * 4,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def NormVCF_cmds(self, depend_cmds, step_name='NormVCF'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['NormVCF'])
        for step, cmd_info in depend_cmds.items():
            sample = cmd_info['sample_name']
            args['input'] = cmd_info['output']
            args['output'] = os.path.join(out_dir, sample+'.norm.vcf')
            cmd = cmdx.NormVCF(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                output=args['output'],
                sample_name=sample,
                mem=1024 ** 3 * 3,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def NGSCheckMate_cmds(self, fastq_info, step_name='NGSCheckMate'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        reformated_fastq_list = os.path.join(self.project_dir, step_name, 'fastq.list')
        with open(fastq_info) as fr, open(reformated_fastq_list, 'w') as fw:
            for line in fr:
                lst = line.strip().split()
                new_lst = lst[1:] + [lst[0]]
                fw.write('\t'.join(new_lst)+'\n')
        args = dict(self.arg_pool['NGSCheckMate'])
        args['fastq_info'] = reformated_fastq_list
        args['out_name'] = 'check_pair'
        args['outdir'] = out_dir
        cmd = cmdx.NGSCheckMate(**args)
        commands[step_name] = self.cmd_dict(
            cmd=cmd,
            output=args['outdir'],
            mem=1024 ** 3 * 4,
            cpu=2,
            monitor_time_step=5
        )
        self.workflow.update(commands)
        return commands
