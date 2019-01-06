import os
import configparser
import argparse
from cmd_generator import *
from nestcmd import RunCommands


arg_ini = 'arguments.ini'
# arg_pool = configparser.ConfigParser()
arg_pool = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
arg_pool.optionxform = str
arg_pool.read(arg_ini, encoding='utf-8')
project_dir = os.path.join(os.getcwd(), 'result')
if not os.path.exists(project_dir):
    os.mkdir(project_dir)
fastq_info_file = 'fastq_info.txt'


def cmd_dict(cmd, cpu=1, mem=200*1024*1024, retry=1, monitor_resource=True,
             monitor_time_step=2, check_resource_before_run=True, **kwargs):
    args = locals()
    args.pop('kwargs')
    if kwargs:
        args.update(kwargs)
    return args


def mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)


def parse_fastq_info(fastq_info_file) -> dict:
    """
    解析fastq输入信息
    :param fastq_info_file:
    format example, at least two column, the third column is optional
    '''
    sample_name | Read1_1_abspath;Read1_2_abspath | Read2_1_abspath;Read2_2_abspath
    sample_name | Read1_abspath | Read2_abspath
    '''
    :return: dict
    """
    fastq_info = dict()
    with open(fastq_info_file) as f:
        for line in f:
            if line.startswith('#') or (not line.strip()):
                pass
            tmp_list = line.strip().split()
            sample, fqs = tmp_list[0], tmp_list[1:]
            fastq_info.setdefault(sample, list())
            read1_list = [x.strip() for x in fqs[0].split(';')]
            fastq_info[sample].append(read1_list)
            if len(fqs) >= 2:
                read2_list = [x.strip() for x in fqs[1].split(';')]
                fastq_info[sample].append(read2_list)
    return fastq_info


def fastqc_raw_data_cmds(fastq_info_dict, step_name='RawDataQC'):
    commands = dict()
    outdir = os.path.join(project_dir, step_name)
    mkdir(outdir)
    for sample, fq_list in fastq_info_dict.items():
        fastqs = ' '.join(x for y in fq_list for x in y)
        args = dict(arg_pool['fastqc'])
        args['outdir'] = outdir
        args['tmpdir'] = outdir
        args['fastqs'] = fastqs
        cmd = fastqc(**args)
        commands[step_name+'_'+sample] = cmd_dict(cmd=cmd)
    return commands


def fastqc_trimmed_data_cmds(trimming_cmds, step_name='TrimmedDataQC'):
    commands = dict()
    outdir = os.path.join(project_dir, step_name)
    mkdir(outdir)
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
        args = dict(arg_pool['fastqc'])
        args['outdir'] = outdir
        args['tmpdir'] = outdir
        args['fastqs'] = ' '.join(trimmed_fq_list)
        cmd = fastqc(**args)
        commands[step_name + '_' + sample] = cmd_dict(
            cmd=cmd, depend=','.join(depend_steps), result_dir=outdir, sample_name=sample,
        )
    return commands


def trimmomatic_cmds(fastq_info_dict, step_name='Trim'):
    commands = dict()
    outdir = os.path.join(project_dir, step_name)
    mkdir(outdir)
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
                trimmed_fq1 = os.path.join(outdir, sample + '_clean_R1.fq.gz')
                trimmed_fq2 = os.path.join(outdir, sample + '_clean_R2.fq.gz')
                unpaired_fq1 = os.path.join(outdir, sample + '_unpaired_clean_R1.fq.gz')
                unpaired_fq2 = os.path.join(outdir, sample + '_unpaired_clean_R2.fq.gz')
                args = dict(arg_pool['trimmomatic'])
                args.update(dict(
                    fq1=fq1, fq2=fq2, mode=mode,
                    trimmed_fq1=trimmed_fq1, trimmed_fq2=trimmed_fq2,
                    unpaired_fq1=unpaired_fq1, unpaired_fq2=unpaired_fq2
                ))
                cmd = trimmomatic(**args)
                commands[step_name+'_'+sample] = cmd_dict(
                    cmd=cmd, sample_name=true_sample,
                    trimmed_fq1=trimmed_fq1, trimmed_fq2=trimmed_fq2
                )
        else:
            for ind, fq1 in enumerate(fq_list[0]):
                if ind >= 1:
                    sample = sample + '_' + str(ind)
                trimmed_fq1 = os.path.join(outdir, sample + '_clean_R1.fq.gz')
                args = dict(arg_pool['trimmomatic'])
                args.update(dict(fq1=fq1, trimmed_fq1=trimmed_fq1, mode=mode))
                cmd = trimmomatic(**args)
                commands[step_name+'_'+sample] = cmd_dict(
                    cmd=cmd, sample_name=true_sample, trimmed_fq1=trimmed_fq1
                )
    return commands


def star_index_cmd(step_name='AlignIndex'):
    commands = dict()
    if os.path.exists(arg_pool['star_index']['genomeDir']):
        print('STAR index existed, and skip this indexing step!')
        return commands
    cmd = star_index(**arg_pool['star_index'])
    commands[step_name] = cmd_dict(cmd=cmd)
    return commands


def star_align_cmds(trimming_cmds, index_cmd, step_name='Align'):
    commands = dict()
    outdir = os.path.join(project_dir, step_name)
    mkdir(outdir)
    samples = set(trimming_cmds[x]['sample_name'] for x in trimming_cmds)
    for sample in samples:
        depend_steps = list() if not index_cmd else list(index_cmd.keys())
        trimmed_fq_list = list()
        for step, cmd_info in trimming_cmds.items():
            if cmd_info['sample_name'] == sample:
                depend_steps.append(step)
                trimmed_fq_list.append(cmd_info['trimmed_fq1'])
                if 'trimmed_fq2' in cmd_info:
                    trimmed_fq_list.append(cmd_info['trimmed_fq2'])
        args = dict(arg_pool['star_align'])
        args['readFilesIn'] = ' '.join(trimmed_fq_list)
        result_dir = os.path.join(outdir, sample)
        mkdir(result_dir)
        prefix = os.path.join(result_dir, sample + '.')
        args['outFileNamePrefix'] = prefix
        cmd = star_align(**args)
        commands[step_name+'_'+sample] = cmd_dict(
            cmd=cmd, mem=1024**3*2, cpu=2,
            monitor_time_step=5, depend=','.join(depend_steps),
            sorted_bam='{}Aligned.sortedByCoord.out.bam'.format(prefix),
            sample_name=sample,
            result_dir=result_dir
        )
    return commands


def star_align_with_rawdata_cmds(fastq_info_dict, index_cmd, step_name='Align'):
    commands = dict()
    outdir = os.path.join(project_dir, step_name)
    mkdir(outdir)
    for sample, fq_list in fastq_info_dict.items():
        depend_steps = list() if not index_cmd else list(index_cmd.keys())
        mode = 'PE' if len(fq_list) == 2 else "SE"
        result_dir = os.path.join(outdir, sample)
        mkdir(result_dir)
        if mode == 'PE':
            fq_list = [x for y in zip(fq_list[0], fq_list[1]) for x in y]
        else:
            fq_list = fq_list[0]
        args = dict(arg_pool['star_align'])
        args['readFilesIn'] = ' '.join(fq_list)
        result_dir = os.path.join(outdir, sample)
        mkdir(result_dir)
        prefix = os.path.join(result_dir, sample + '_')
        args['outFileNamePrefix'] = prefix
        cmd = star_align(**args)
        commands[step_name + '_' + sample] = cmd_dict(
            cmd=cmd, mem=1024 ** 3 * 2, cpu=2,
            monitor_time_step=5, depend=','.join(depend_steps),
            sorted_bam='{}Aligned.sortedByCoord.out.bam'.format(prefix),
            sample_name=sample,
            result_dir=result_dir
        )
    return commands


def scallop_cmds(align_cmds, step_name='Assembly'):
    commands = dict()
    outdir = os.path.join(project_dir, step_name)
    mkdir(outdir)
    for step, cmd_info in align_cmds.items():
        if 'sample_name' not in cmd_info:
            continue
        sample = cmd_info['sample_name']
        args = dict(arg_pool['scallop'])
        args['bam'] = cmd_info['sorted_bam']
        result_dir = os.path.join(outdir, sample)
        mkdir(result_dir)
        args['out_gtf'] = os.path.join(result_dir, '{}.scallop.gtf'.format(sample))
        cmd = scallop(**args)
        commands[step_name+'_'+sample] = cmd_dict(
            cmd=cmd, mem=1024**3*2, cpu=2, monitor_time_step=5,
            depend=step,
            sample_name=sample,
            out_gtf=args['out_gtf']
        )
    return commands


def merge_scallop_transcripts_cmd(assemble_cmds, step_name='MergeTranscript'):
    gtfs = list()
    for step, cmd_info in assemble_cmds.items():
        gtfs.append(cmd_info['out_gtf'])
    args = dict(arg_pool['merge_scallop_transcripts'])
    args['gtf'] = ' '.join(gtfs)
    outdir = os.path.join(project_dir, step_name)
    mkdir(outdir)
    args['outdir'] = outdir
    cmd = merge_scallop_transcripts(**args)
    commands = dict()
    commands[step_name] = cmd_dict(
        cmd=cmd, cpu=2, depend=','.join(assemble_cmds.keys()),
        result_dir=outdir, all_transcripts=os.path.join(outdir, 'all.transcripts.fa'),
        all_gtf=os.path.join(outdir, 'all.transcripts.gtf')
    )
    return commands


def salmon_index_cmd(step_name='QuantIndex', merge_transcript_cmd=None):
    commands = dict()
    if not merge_transcript_cmd:
        if os.path.exists(arg_pool['salmon_index']['index_prefix']):
            print('salmon index existed, and skip this indexing step!')
            return commands
        args = dict(arg_pool['salmon_index'])
        cmd = salmon_index(**args)
        commands[step_name] = cmd_dict(cmd=cmd)
    else:
        outdir = os.path.join(project_dir, step_name)
        mkdir(outdir)
        args = dict(arg_pool['salmon_index'])
        merge_step = list(merge_transcript_cmd.keys())[0]
        args['transcript_fasta'] = merge_transcript_cmd[merge_step]['all_transcripts']
        args['index_prefix'] = outdir
        cmd = salmon_index(**args)
        commands[step_name] = cmd_dict(
            cmd=cmd, index_prefix=outdir, depend=merge_step,
            all_gtf=merge_transcript_cmd[merge_step]['all_gtf']
        )
    return commands


def salmon_quant_with_raw_data_cmds(fastq_info_dict, index_cmd, step_name='Quant'):
    commands = dict()
    result_dir = os.path.join(project_dir, step_name)
    mkdir(result_dir)
    for sample, fq_list in fastq_info_dict.items():
        depend_steps = list() if not index_cmd else list(index_cmd.keys())
        mode = 'PE' if len(fq_list) == 2 else "SE"
        if mode == 'PE':
            fq1_list = fq_list[0]
            fq2_list = fq_list[1]
            if len(fq1_list) != len(fq2_list):
                raise Exception('fastq文件的read1和read2必须一一对应')
            args = dict(arg_pool['salmon_quant'])
            if index_cmd:
                if 'all_gtf' in list(index_cmd.values())[0]:
                    args['transcript2gene'] = list(index_cmd.values())[0]['all_gtf']
            args['fq1'] = ' '.join(fq1_list)
            args['fq2'] = ' '.join(fq2_list)
            args['mode'] = mode
            prefix = os.path.join(result_dir, sample)
            args['out_prefix'] = prefix
            cmd = salmon_quant(**args)
            commands[step_name + '_' + sample] = cmd_dict(
                cmd=cmd, mem=1024 ** 3 * 2, cpu=2, monitor_time_step=5,
                depend=','.join(depend_steps),
                sample_name=sample,
                result_dir=result_dir,
                prefix=prefix
            )
        else:
            args = dict(arg_pool['salmon_quant'])
            if index_cmd:
                if 'all_gtf' in index_cmd.values()[0]:
                    args['transcript2gene'] = list(index_cmd.values())[0]['all_gtf']
            args['fq'] = ' '.join(fq_list[0])
            args['mode'] = mode
            prefix = os.path.join(result_dir, sample)
            args['out_prefix'] = prefix
            cmd = salmon_quant(**args)
            commands[step_name + '_' + sample] = cmd_dict(
                cmd=cmd, mem=1024 ** 3 * 2, cpu=2, monitor_time_step=5,
                depend=','.join(depend_steps),
                sample_name=sample,
                result_dir=result_dir,
                prefix=prefix
            )
    return commands


def salmon_quant_with_clean_data_cmds(trimming_cmds, index_cmd, step_name='Quant'):
    commands = dict()
    result_dir = os.path.join(project_dir, step_name)
    mkdir(result_dir)
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
            args = dict(arg_pool['salmon_quant'])
            if index_cmd:
                if 'all_gtf' in list(index_cmd.values())[0]:
                    args['transcript2gene'] = list(index_cmd.values())[0]['all_gtf']
            args['fq1'] = ' '.join(trimmed_fq1_list)
            args['fq2'] = ' '.join(trimmed_fq2_list)
            args['mode'] = mode
            prefix = os.path.join(result_dir, sample)
            args['out_prefix'] = prefix
            cmd = salmon_quant(**args)
            commands[step_name + '_' + sample] = cmd_dict(
                cmd=cmd, mem=1024 ** 3 * 2, cpu=2, monitor_time_step=5,
                depend=','.join(depend_steps),
                sample_name=sample,
                result_dir=result_dir,
                prefix=prefix
            )
        else:
            args = dict(arg_pool['salmon_quant'])
            if index_cmd:
                if 'all_gtf' in list(index_cmd.values())[0]:
                    args['transcript2gene'] = list(index_cmd.values())[0]['all_gtf']
            args['fq'] = ' '.join(trimmed_fq1_list)
            args['mode'] = mode
            prefix = os.path.join(result_dir, sample)
            args['out_prefix'] = prefix
            cmd = salmon_quant(**args)
            commands[step_name + '_' + sample] = cmd_dict(
                cmd=cmd, mem=1024 ** 3 * 2, cpu=2, monitor_time_step=5,
                depend=','.join(depend_steps),
                sample_name=sample,
                result_dir=result_dir,
                prefix=prefix
            )
    return commands


def merge_quant_cmds(quant_cmds, step_name='MergeQuant', quant_method='salmon'):
    commands = dict()
    depend = list()
    out_dir = os.path.join(project_dir, step_name)
    mkdir(out_dir)
    result_dir = ''
    for step, cmd_info in quant_cmds.items():
        if 'sample_name' not in cmd_info:
            continue
        result_dir = cmd_info['result_dir']
        depend.append(step)
    args = dict(arg_pool['abundance_estimates_to_matrix'])
    args['est_method'] = quant_method
    args['quant_result'] = result_dir + '/*/quant.sf'
    args['out_prefix'] = os.path.join(out_dir, 'all')
    cmd = abundance_estimates_to_matrix(**args)
    commands[step_name + 'Transcript'] = cmd_dict(
        cmd=cmd, monitor_resource=False, check_resource_before_run=False,
        depend=','.join(depend), out_prefix=args['out_prefix'],
        transcript_tpm_matrix=args['out_prefix'] + '.isoform.TMM.EXPR.matrix',
        transcript_count_matrix=args['out_prefix'] + '.isoform.counts.matrix',
    )
    args['quant_result'] = result_dir + '/*/quant.genes.sf'
    cmd = abundance_estimates_to_matrix(**args)
    commands[step_name + 'Gene'] = cmd_dict(
        cmd=cmd, monitor_resource=False, check_resource_before_run=False,
        depend=','.join(depend), out_prefix=args['out_prefix'],
        gene_tpm_matrix=args['out_prefix'] + '.gene.TMM.EXPR.matrix',
        gene_count_matrix=args['out_prefix'] + '.gene.counts.matrix',
    )
    return commands


def pipeline():
    commands = configparser.ConfigParser()
    commands.optionxform = str
    commands['mode'] = dict(
        threads=3,
        retry=1,
        monitor_resource=True,
        monitor_time_step=2,
        check_resource_before_run=True,
    )
    fastq_info_dict = parse_fastq_info(fastq_info_file)
    commands.update(fastqc_raw_data_cmds(fastq_info_dict))
    trim_cmds = trimmomatic_cmds(fastq_info_dict)
    commands.update(trim_cmds)
    fastqc_cmds = fastqc_trimmed_data_cmds(trimming_cmds=trim_cmds)
    commands.update(fastqc_cmds)
    star_indexing = star_index_cmd()
    commands.update(star_indexing)
    align_cmds = star_align_cmds(trimming_cmds=trim_cmds, index_cmd=star_indexing)
    commands.update(align_cmds)
    assembly_cmds = scallop_cmds(align_cmds=align_cmds)
    commands.update(assembly_cmds)
    merge_trans_cmd = merge_scallop_transcripts_cmd(assembly_cmds)
    commands.update(merge_trans_cmd)
    salmon_indexing = salmon_index_cmd(merge_transcript_cmd=merge_trans_cmd)
    commands.update(salmon_indexing)
    quant_cmds = salmon_quant_with_clean_data_cmds(trimming_cmds=trim_cmds, index_cmd=salmon_indexing)
    commands.update(quant_cmds)
    merge_cmds = merge_quant_cmds(quant_cmds=quant_cmds)
    commands.update(merge_cmds)
    with open('pipeline_cmds.ini', 'w') as configfile:
        commands.write(configfile)
    workflow = RunCommands('pipeline_cmds.ini')
    workflow.parallel_run()


if __name__ == '__main__':
    pipeline()
