import os
import sys
sys.path.append("..")
from basic.workflow_basic import Basic
import dnaseq_pipeline.one_cmd_generator as cmdx


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

    def fastqc_raw_data_cmds(self, fastq_info, step_name='RawDataQC'):
        fastq_info_dict = self.parse_fastq_info(fastq_info)
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for sample, fq_list in fastq_info_dict.items():
            warn = 'As bwa does not support a sample with multiple fastqs, you should cat them first!'
            if len(fq_list) == 2:
                if len(fq_list[0]) >= 2 or len(fq_list[1]) >=2:
                    raise Exception(warn)
            else:
                if len(fq_list[0]) >= 2:
                    raise Exception(warn)
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

    def trimmomatic_cmds(self, fastq_info, step_name='Trim'):
        fastq_info_dict = self.parse_fastq_info(fastq_info)
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

    def bwa_align_rawdata_cmds(self, fastq_info, step_name='BWAmem'):
        fastq_info_dict = self.parse_fastq_info(fastq_info)
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for sample, fq_list in fastq_info_dict.items():
            args = dict(self.arg_pool['bwa_mem'])
            args['out'] = os.path.join(outdir, '{}.bwa.sam'.format(sample))
            args['fq1'] = fq_list[0][0]
            args['sample'] = sample
            if len(fq_list) >= 2:
                args['fq2'] = fq_list[1][0]
            cmd = cmdx.bwa_mem(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                mem=1024 ** 3 * 10,
                cpu=10,
                monitor_time_step=5,
                sample_name=sample,
                outdir=outdir,
                out=args['out']
            )
        self.workflow.update(commands)
        return commands

    def bwa_align_trimmed_data_cmds(self, depend_cmds, step_name='BWAmem'):
        commands = dict()
        outdir = os.path.join(self.project_dir, step_name)
        self.mkdir(outdir)
        for step, cmd_info in depend_cmds.items():
            sample = cmd_info['sample_name']
            args = dict(self.arg_pool['bwa_mem'])
            args['out'] = os.path.join(outdir, '{}.bwa.sam'.format(sample))
            args['fq1'] = cmd_info['trimmed_fq1']
            args['sample'] = sample
            if 'trimmed_fq2' in cmd_info:
                args['fq2'] = cmd_info['trimmed_fq2']
            cmd = cmdx.bwa_mem(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                mem=1024 ** 3 * 10,
                cpu=10,
                monitor_time_step=5,
                sample_name=sample,
                outdir=outdir,
                out=args['out']
            )
        self.workflow.update(commands)
        return commands

    def sam2bam_cmds(self, bwa_cmds, step_name='Sam2Bam'):
        commands = dict()
        for step, cmd_info in bwa_cmds.items():
            sample = cmd_info['sample_name']
            args = dict(self.arg_pool['sam2bam'])
            args['sam'] = cmd_info['out']
            args['out'] = cmd_info['out'][:-3] + 'bam'
            cmd = cmdx.sam2bam(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                mem=1024 ** 3 * 5, cpu=2, monitor_time_step=5,
                depend=step,
                sample_name=sample,
                out=args['out']
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
            args['out'] = os.path.join(outdir, f'{sample}.unmapped.bam')
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

    def RawFastqToSam_cmds(self, fastq_info, step_name='FastqToSam'):
        fastq_info_dict = self.parse_fastq_info(fastq_info)
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
            args['out'] = os.path.join(outdir, f'{sample}.unmapped.bam')
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

    def MergeBamAlignment_cmds(self, sam2bam_cmds, fastq2bam_cmds, step_name='MergeBamAlignment'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['MergeBamAlignment'])
        for step, cmd_info in sam2bam_cmds.items():
            sample = cmd_info['sample_name']
            args['ALIGNED'] = cmd_info['out']
            another_step = [x for x in fastq2bam_cmds if x.endswith(sample)][0]
            args['UNMAPPED'] = fastq2bam_cmds[another_step]['out']
            args['out'] = os.path.join(out_dir, sample+'.merged.bam')
            cmd = cmdx.MergeBamAlignment(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=','.join([step, another_step]),
                out=args['out'],
                sample_name=sample,
                mem=1024 ** 3 * 3,
                cpu=3,
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
            if self.workflow_arguments.disable_markdup_spark:
                cmd = cmdx.MarkDuplicates(**args)
            else:
                cmd = cmdx.MarkDuplicatesSpark(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                output=args['output'],
                sample_name=sample,
                mem=1024 ** 3 * 5,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def SortAndFixTags_cmds(self, mark_dup_cmds, step_name='SortAndFixTags'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['SortAndFixTags'])
        for step, cmd_info in mark_dup_cmds.items():
            sample = cmd_info['sample_name']
            args['input'] = cmd_info['output']
            args['output'] = os.path.join(out_dir, sample+'.markdup.sortedByCoord.bam')
            if self.workflow_arguments.disable_markdup_spark:
                cmd = cmdx.SortAndFixTags(**args)
            else:
                cmd = cmdx.FixTags(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                output=args['output'],
                sample_name=sample,
                mem=1024 ** 3 * 5,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def BaseRecalibrator_cmds(self, sort_cmds, step_name='BaseRecalibrator'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['BaseRecalibrator'])
        for step, cmd_info in sort_cmds.items():
            sample = cmd_info['sample_name']
            args['input'] = cmd_info['output']
            seq_group_lst = [x.strip() for x in open(args['seq_grouping_file'])]
            for ind, group in enumerate(seq_group_lst):
                args['output'] = os.path.join(out_dir, sample+'.{}.recal.table'.format(ind))
                args['intervals'] = group.split()
                cmd = cmdx.BaseRecalibrator(**args)
                commands[step_name + '_' + sample + '_{}'.format(ind)] = self.cmd_dict(
                    cmd=cmd,
                    depend=step,
                    output=args['output'],
                    sample_name=sample,
                    mem=1024 ** 3 * 3,
                    cpu=2,
                    monitor_time_step=5,
                    sorted_bam=args['input']
                )
        self.workflow.update(commands)
        return commands

    def GatherBQSRReports_cmds(self, recal_cmds, step_name='GatherBQSRReports'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['GatherBQSRReports'])
        samples = []
        for step, cmd_info in recal_cmds.items():
            sample = cmd_info['sample_name']
            if sample not in samples:
                samples.append(sample)
        for sample in samples:
            depend_steps = []
            report_lst = []
            sorted_bam = ''
            for step, cmd_info in recal_cmds.items():
                if sample == cmd_info['sample_name']:
                    depend_steps.append(step)
                    report_lst.append(cmd_info['output'])
                    sorted_bam = cmd_info['sorted_bam']
            args['report_list'] = report_lst
            args['output'] = os.path.join(out_dir, sample + '.all.recal.table')
            cmd = cmdx.GatherBQSRReports(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=','.join(depend_steps),
                output=args['output'],
                sample_name=sample,
                sorted_bam=sorted_bam,
                mem=1024 ** 3 * 5,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def ApplyBQSR_cmds(self, gather_report_cmds, step_name='ApplyBQSR'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['ApplyBQSR'])
        for step, cmd_info in gather_report_cmds.items():
            sample = cmd_info['sample_name']
            args['input'] = cmd_info['sorted_bam']
            args['bqsr-recal-file'] = cmd_info['output']
            seq_group_lst = [x.strip() for x in open(args['seq_grouping_file'])]
            for ind, group in enumerate(seq_group_lst):
                args['output'] = os.path.join(out_dir, sample + '.{}.bqsr.bam'.format(ind))
                args['intervals'] = group.split()
                cmd = cmdx.ApplyBQSR(**args)
                commands[step_name + '_' + sample + '_{}'.format(ind)] = self.cmd_dict(
                    cmd=cmd,
                    depend=step,
                    output=args['output'],
                    sample_name=sample,
                    mem=1024 ** 3 * 4,
                    cpu=2,
                    monitor_time_step=5
                )
        self.workflow.update(commands)
        return commands

    def GatherBamFiles_cmds(self, apply_bqsr_cmds, step_name='GatherBamFiles'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['GatherBamFiles'])
        samples = []
        for step, cmd_info in apply_bqsr_cmds.items():
            sample = cmd_info['sample_name']
            if sample not in samples:
                samples.append(sample)
        for sample in samples:
            depend_steps = []
            bam_lst = []
            for step, cmd_info in apply_bqsr_cmds.items():
                if sample == cmd_info['sample_name']:
                    depend_steps.append(step)
                    bam_lst.append(cmd_info['output'])
            args['bam_list'] = bam_lst
            args['output'] = os.path.join(out_dir, sample + '.final.bqsr.bam')
            cmd = cmdx.GatherBamFiles(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=','.join(depend_steps),
                output=args['output'],
                sample_name=sample,
                mem=1024 ** 3 * 8,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def samtools_index_cmds(self, GatherBamFiles_cmds, step_name='IndexGatheredBam'):
        commands = dict()
        args = dict(self.arg_pool['samtools_index'])
        for step, cmd_info in GatherBamFiles_cmds.items():
            sample = cmd_info['sample_name']
            args['bam'] = cmd_info['output']
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmdx.samtools_index(**args),
                mem=1024 ** 3 * 5, cpu=2, monitor_time_step=5,
                depend=step,
                sample_name=sample,
                output=args['bam']
            )
        self.workflow.update(commands)
        return commands

    def SplitIntervals_cmds(self, step_name='SplitIntervals'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['SplitIntervals'])
        args['output'] = out_dir
        out_intervals = []
        for each in range(int(args['scatter-count'])):
            out_intervals.append(os.path.join(out_dir, '{:0>4}-scattered.interval_list'.format(each)))
        commands[step_name] = self.cmd_dict(
            cmd=cmdx.SplitIntervals(**args),
            output=args['output'],
            mem=1024 ** 3 * 2,
            cpu=2,
            monitor_time_step=5,
            intervals = out_intervals
        )
        self.workflow.update(commands)
        return commands

    def MuTect2_cmds(self, pair_info, GatherBamFiles_cmds, split_intervals_cmds, step_name='MuTect2'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['MuTect2'])
        # parse pair info
        pair_info = [x.strip().split() for x in open(pair_info)]
        used_sample = []
        for pair in pair_info:
            depends = list(split_intervals_cmds.keys())
            step1 = [x for x in GatherBamFiles_cmds.keys() if x.endswith(pair[0])][0]
            depends.append(step1)
            args['tumour_bam'] = GatherBamFiles_cmds[step1]['output']
            sample = GatherBamFiles_cmds[step1]['sample_name']
            if sample in used_sample:
                raise Exception(f'Tumour sample id {sample} duplicated!')
            else:
                used_sample.append(sample)
            args['tumour-sample'] = sample
            call_mode = 'Single'
            if len(pair) >= 2:
                call_mode = 'Paired'
                step2 = [x for x in GatherBamFiles_cmds.keys() if x.endswith(pair[1])][0]
                depends.append(step2)
                args['normal_bam'] = GatherBamFiles_cmds[step2]['output']
                args['normal-sample'] = GatherBamFiles_cmds[step2]['sample_name']
            for interval in split_intervals_cmds[depends[0]]['intervals']:
                ind = os.path.basename(interval).split('-')[0]
                args['intervals'] = interval
                args['output'] = os.path.join(out_dir, sample+'.{}.vcf'.format(ind))
                args['bam-output'] = os.path.join(out_dir, sample + '.mutect2.{}.bam'.format(ind))
                args['f1r2-tar-gz'] = os.path.join(out_dir, sample + '.f1r2.{}.tar.gz'.format(ind))
                cmd = cmdx.MuTect2(**args)
                commands[step_name + call_mode + '_' + sample + '_{}'.format(ind)] = self.cmd_dict(
                    cmd=cmd,
                    depend=','.join(depends),
                    output=args['output'],
                    sample_name=sample,
                    fake_sample_name=sample+call_mode,
                    interval_ind=ind,
                    interval=interval,
                    tumour_bam=args['tumour_bam'],
                    normal_bam=args['normal_bam'] if len(pair) >= 2 else '',
                    mem=1024 ** 3 * 10,
                    cpu=2,
                    monitor_time_step=5,
                    outbam=args['bam-output'],
                    f1r2=args['f1r2-tar-gz'],
                    caller_type='somatic',
                    call_mode=call_mode
                )
        self.workflow.update(commands)
        return commands

    def HaplotypeCaller_cmds(self, pair_info, GatherBamFiles_cmds, split_intervals_cmds, step_name='Haplotype'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['HaplotypeCaller'])
        # parse pair info
        pair_info = [x.strip().split() for x in open(pair_info)]
        used_sample = []
        for pair in pair_info:
            if len(pair) != 1:
                continue
            depends = list(split_intervals_cmds.keys())
            step1 = [x for x in GatherBamFiles_cmds.keys() if x.endswith(pair[0])][0]
            depends.append(step1)
            args['input'] = GatherBamFiles_cmds[step1]['output']
            sample = GatherBamFiles_cmds[step1]['sample_name']
            if sample in used_sample:
                raise Exception(f'sample id {sample} duplicated!')
            else:
                used_sample.append(sample)

            for interval in split_intervals_cmds[depends[0]]['intervals']:
                ind = os.path.basename(interval).split('-')[0]
                args['intervals'] = interval
                args['output'] = os.path.join(out_dir, sample + '.{}.vcf'.format(ind))
                args['bamout'] = os.path.join(out_dir, sample + '.haplotype.{}.bam'.format(ind))
                cmd = cmdx.HaplotypeCaller(**args)
                commands[step_name + '_' + sample + '_{}'.format(ind)] = self.cmd_dict(
                    cmd=cmd,
                    depend=','.join(depends),
                    output=args['output'],
                    sample_name=sample,
                    fake_sample_name=sample,
                    interval_ind=ind,
                    interval=interval,
                    mem=1024 ** 3 * 10,
                    cpu=2,
                    monitor_time_step=5,
                    outbam=args['bamout'],
                    caller_type='germline'
                )
        self.workflow.update(commands)
        return commands

    def GetPileupSummaries_cmds(self, GatherBamFiles_cmds, step_name='GetPileup'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['GetPileupSummaries'])
        for step, cmd_info in GatherBamFiles_cmds.items():
            sample = cmd_info['sample_name']
            args['bam'] = cmd_info['output']
            args['output'] = os.path.join(out_dir, sample+'.pileups.table')
            cmd = cmdx.GetPileupSummaries(**args)
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmd,
                depend=step,
                output=args['output'],
                sample_name=sample,
                mem=1024 ** 3 * 5,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def LearnReadOrientationModel_cmds(self, mutect2_cmds, step_name='LROM'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['LearnReadOrientationModel'])
        samples = []
        for _, cmd_info in mutect2_cmds.items():
            sample = cmd_info['sample_name']
            if sample not in samples:
                samples.append(sample)
        for sample in samples:
            f1r2_lst = []
            depends = []
            for step, cmd_info in mutect2_cmds.items():
                if cmd_info['sample_name'] == sample:
                    f1r2_lst.append(cmd_info['f1r2'])
                    depends.append(step)
                    tumour_bam = cmd_info['tumour_bam']
            args['input'] = f1r2_lst
            args['output'] = os.path.join(out_dir, sample+'.artifact-priors.tar.gz')
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmdx.LearnReadOrientationModel(**args),
                depend=','.join(depends),
                output=args['output'],
                tumour_bam=tumour_bam,
                sample_name=sample,
                mem=1024 ** 3 * 8,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def MergeVcfs_cmds(self, mutect2_cmds, step_name='MergeVcfs'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['MergeVcfs'])
        samples = []
        for _, cmd_info in mutect2_cmds.items():
            sample = cmd_info['sample_name']
            if sample not in samples:
                samples.append(sample)
        for sample in samples:
            vcf_lst = []
            depends = []
            for step, cmd_info in mutect2_cmds.items():
                if cmd_info['sample_name'] == sample:
                    vcf_lst.append(cmd_info['output'])
                    depends.append(step)
            args['input'] = vcf_lst
            args['output'] = os.path.join(out_dir, sample + '.raw.vcf.gz')
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmdx.MergeVcfs(**args),
                depend=','.join(depends),
                output=args['output'],
                sample_name=sample,
                mem=1024 ** 3 * 3,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def GatherMutect2Bam_cmds(self, mutect2_cmds, step_name='MergeMutect2Bam'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['GatherBamFiles'])
        samples = []
        for _, cmd_info in mutect2_cmds.items():
            sample = cmd_info['sample_name']
            if sample not in samples:
                samples.append(sample)

        for sample in samples:
            bam_lst = []
            depends = []
            for step, cmd_info in mutect2_cmds.items():
                if cmd_info['sample_name'] == sample:
                    bam_lst.append(cmd_info['outbam'])
                    depends.append(step)
            args['bam_list'] = bam_lst
            args['CREATE_INDEX'] = 'false'
            args['CREATE_MD5_FILE'] = 'false'
            args['output'] = os.path.join(out_dir, sample + '.mutect2.unsorted.bam')
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmdx.GatherBamFiles(**args),
                depend=','.join(depends),
                output=args['output'],
                sample_name=sample,
                outdir=out_dir,
                mem=1024 ** 3 * 6,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def SortSam_cmds(self, merge_bam_cmds, step_name='SortBam'):
        commands = dict()
        args = dict(self.arg_pool['SortSam'])
        for step, cmd_info in merge_bam_cmds.items():
            outdir = cmd_info['outdir']
            sample = cmd_info['sample_name']
            args['input'] = cmd_info['output']
            args['output'] = os.path.join(outdir, sample+'.sorted.bam')
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmdx.SortSam(**args),
                mem=1024 ** 3 * 5, cpu=2, monitor_time_step=5,
                depend=step,
                sample_name=sample,
                outdir=outdir,
                output=args['output']
            )
        self.workflow.update(commands)
        return commands

    def BuildBamIndex_cmds(self, sort_bam_cmds, step_name='IndexBam'):
        commands = dict()
        args = dict(self.arg_pool['BuildBamIndex'])
        for step, cmd_info in sort_bam_cmds.items():
            outdir = cmd_info['outdir']
            sample = cmd_info['sample_name']
            args['input'] = cmd_info['output']
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmdx.BuildBamIndex(**args),
                mem=1024 ** 3 * 5, cpu=2, monitor_time_step=5,
                depend=step,
                sample_name=sample,
                outdir=outdir,
                sorted_bam=args['input']
            )
        self.workflow.update(commands)
        return commands

    def MergeMutectStats_cmds(self, mutect2_cmds, step_name='MergeVcfStats'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['MergeMutectStats'])
        samples = []
        for _, cmd_info in mutect2_cmds.items():
            sample = cmd_info['sample_name']
            if sample not in samples:
                samples.append(sample)
        for sample in samples:
            vcf_lst = []
            depends = []
            for step, cmd_info in mutect2_cmds.items():
                if cmd_info['sample_name'] == sample:
                    vcf_lst.append(cmd_info['output']+'.stats')
                    depends.append(step)
            args['stats'] = vcf_lst
            args['output'] = os.path.join(out_dir, sample + '.vcf.stats')
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmdx.MergeMutectStats(**args),
                depend=','.join(depends),
                output=args['output'],
                sample_name=sample,
                mem=1024 ** 3 * 3,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def CalculateContamination_cmds(self, pair_info, get_pileup_cmds, step_name='CalcContamination'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['CalculateContamination'])
        # parse pair info
        pair_info = [x.strip().split() for x in open(pair_info)]
        for pair in pair_info:
            depends = []
            step1 = [x for x in get_pileup_cmds.keys() if x.endswith(pair[0])][0]
            depends.append(step1)
            args['tumor_pileups'] = get_pileup_cmds[step1]['output']
            if len(pair) >= 2:
                step2 = [x for x in get_pileup_cmds.keys() if x.endswith(pair[1])][0]
                depends.append(step2)
                args['normal_pileups'] = get_pileup_cmds[step2]['output']
            sample = get_pileup_cmds[step1]['sample_name']
            args['tumor-segmentation'] = os.path.join(out_dir, sample+'.segment.table')
            args['output'] = os.path.join(out_dir, sample+'.contamination.table')
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmdx.CalculateContamination(**args),
                depend=','.join(depends),
                output=args['output'],
                sample_name=sample,
                segment=args['tumor-segmentation'],
                mem=1024 ** 3 * 5,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def FilterMutectCalls_cmds(self, merge_vcf_cmds, merge_stats_cmds, lrom_cmds,
                               cal_contamination_cmds, step_name='FilterVcf'):
        commands = dict()
        out_dir = os.path.join(self.project_dir, step_name)
        self.mkdir(out_dir)
        args = dict(self.arg_pool['FilterMutectCalls'])
        samples = []
        for _, cmd_info in merge_vcf_cmds.items():
            sample = cmd_info['sample_name']
            if sample not in samples:
                samples.append(sample)
        for sample in samples:
            step1 = [x for x in merge_vcf_cmds if merge_vcf_cmds[x]['sample_name'] == sample][0]
            step2 = [x for x in merge_stats_cmds if merge_stats_cmds[x]['sample_name'] == sample][0]
            step3 = [x for x in lrom_cmds if lrom_cmds[x]['sample_name'] == sample][0]
            step4 = [x for x in cal_contamination_cmds if cal_contamination_cmds[x]['sample_name']==sample][0]
            args['variant'] = merge_vcf_cmds[step1]['output']
            args['stats'] = merge_stats_cmds[step2]['output']
            args['ob-priors'] = lrom_cmds[step3]['output']
            args['contamination-table'] = cal_contamination_cmds[step4]['output']
            args['tumor-segmentation'] = cal_contamination_cmds[step4]['segment']
            args['output'] = os.path.join(out_dir, sample+'.filtered.vcf.gz')
            args['filtering-stats'] = os.path.join(out_dir, sample+'.filtering.stats')
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmdx.FilterMutectCalls(**args),
                depend=','.join([step1, step2, step3, step4]),
                output=args['output'],
                outdir=out_dir,
                stats=args['filtering-stats'],
                sample_name=sample,
                tumour_bam=lrom_cmds[step3]['tumour_bam'],
                mem=1024 ** 3 * 5,
                cpu=2,
                monitor_time_step=5
            )
        self.workflow.update(commands)
        return commands

    def FilterAlignmentArtifacts_cmds(self, filter_calls_cmds, step_name='FilterArtifacts'):
        commands = dict()
        args = dict(self.arg_pool['FilterAlignmentArtifacts'])
        for step, cmd_info in filter_calls_cmds.items():
            outdir = cmd_info['outdir']
            sample = cmd_info['sample_name']
            args['variant'] = cmd_info['output']
            args['bam'] = cmd_info['tumour_bam']
            args['output'] = os.path.join(outdir, sample+'.filtered2.vcf.gz')
            commands[step_name + '_' + sample] = self.cmd_dict(
                cmd=cmdx.FilterAlignmentArtifacts(**args),
                mem=1024 ** 3 * 5, cpu=2, monitor_time_step=5,
                depend=step,
                sample_name=sample,
                outdir=outdir,
                sorted_bam=args['output']
            )
        self.workflow.update(commands)
        return commands

