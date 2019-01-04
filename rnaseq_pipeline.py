import os
import configparser
import argparse
from cmd_generator import *


arg_ini = 'arguments.ini'
# arg_pool = configparser.ConfigParser()
arg_pool = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
arg_pool.optionxform = str
arg_pool.read(arg_ini, encoding='utf-8')
project_dir = os.getcwd()
fastq_info_file = 'fastq_info.txt'



def cmd_dict(cmd, cpu=1, mem=100*1024*1024, retry=1, monitor_resource=True,
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


def fastqc_cmds(fastq_info_dict, step_name='RawDataQC'):
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
    cmd = star_index(**arg_pool['star_index'])
    commands[step_name] = cmd_dict(cmd=cmd)
    return commands


def star_align_cmds(trimming_cmds, step_name='Align', index_step_name='AlignIndex'):
    commands = dict()
    index_exist = True
    if os.path.exists(arg_pool['star_index']['genomeDir']):
        print('STAR index existed, and skip this indexing step!')
    else:
        index_exist = False
        # print('building index for star...')
        commands = star_index_cmd(step_name=index_step_name)
    outdir = os.path.join(project_dir, step_name)
    mkdir(outdir)
    samples = set(trimming_cmds[x]['sample_name'] for x in trimming_cmds)
    for sample in samples:
        depend_steps = list() if index_exist else [index_step_name]
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


def star_align_with_rawdata_cmds(fastq_info_dict, step_name='Align', index_step_name='AlignIndex'):
    commands = dict()
    index_exist = True
    if os.path.exists(arg_pool['star_index']['genomeDir']):
        print('STAR index existed, and skip this indexing step!')
    else:
        index_exist = False
        # print('building index for star...')
        commands = star_index_cmd(step_name=index_step_name)
    outdir = os.path.join(project_dir, step_name)
    mkdir(outdir)
    for sample, fq_list in fastq_info_dict.items():
        depend_steps = list() if index_exist else [index_step_name]
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


def salmon_index_cmd(step_name='QuantIndex'):
    commands = dict()
    cmd = salmon_index(**arg_pool['salmon_index'])
    commands[step_name] = cmd_dict(cmd=cmd)
    return commands


def salmon_quant_with_rawdata_cmds(fastq_info_dict, step_name='Quant', index_step_name='QuantIndex'):
    commands = dict()
    index_exist = True
    if os.path.exists(arg_pool['salmon_index']['index']):
        print('salmon index existed, and skip this indexing step!')
    else:
        index_exist = False
        # print('building index for salmon...')
        commands = salmon_index_cmd(step_name=index_step_name)
    outdir = os.path.join(project_dir, step_name)
    mkdir(outdir)
    for sample, fq_list in fastq_info_dict.items():
        depend_steps = list() if index_exist else [index_step_name]
        mode = 'PE' if len(fq_list) == 2 else "SE"
        result_dir = os.path.join(outdir, sample)
        mkdir(result_dir)
        if mode == 'PE':
            fq1_list = fq_list[0]
            fq2_list = fq_list[1]
            if len(fq1_list) != len(fq2_list):
                raise Exception('fastq文件的read1和read2必须一一对应')
            fq1 = ' '.join(fq1_list)
            fq2 = ' '.join(fq2_list)
            args = dict(arg_pool['salmon_quant'])
            args['fq1'] = fq1
            args['fq2'] = fq2
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
            fqs = ' '.join(fq_list[0])
            args = dict(arg_pool['salmon_quant'])
            args['fq'] = fqs
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


def salmon_quant_cmds(trimming_cmds, step_name='Quant', index_step_name='QuantIndex'):
    commands = dict()
    index_exist = True
    if os.path.exists(arg_pool['salmon_index']['index_prefix']):
        print('salmon index existed, and skip this indexing step!')
    else:
        index_exist = False
        # print('building index for salmon...')
        commands = salmon_index_cmd(step_name=index_step_name)
    outdir = os.path.join(project_dir, step_name)
    mkdir(outdir)
    samples = set(trimming_cmds[x]['sample_name'] for x in trimming_cmds)
    for sample in samples:
        depend_steps = list() if index_exist else [index_step_name]
        trimmed_fq1_list = list()
        trimmed_fq2_list = list()
        for step, cmd_info in trimming_cmds.items():
            if cmd_info['sample_name'] == sample:
                depend_steps.append(step)
                trimmed_fq1_list.append(cmd_info['trimmed_fq1'])
                if 'trimmed_fq2' in cmd_info:
                    trimmed_fq2_list.append(cmd_info['trimmed_fq2'])
        mode = 'PE' if trimmed_fq2_list else 'SE'
        result_dir = os.path.join(outdir, sample)
        mkdir(result_dir)
        if mode == 'PE':
            if len(trimmed_fq1_list) != len(trimmed_fq2_list):
                raise Exception('fastq文件的read1和read2必须一一对应')
            fq1 = ' '.join(trimmed_fq1_list)
            fq2 = ' '.join(trimmed_fq2_list)
            args = dict(arg_pool['salmon_quant'])
            args['fq1'] = fq1
            args['fq2'] = fq2
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
            fqs = ' '.join(trimmed_fq1_list)
            args = dict(arg_pool['salmon_quant'])
            args['fq'] = fqs
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


def pipeline():
    commands = configparser.ConfigParser()
    commands.optionxform = str
    commands['mode'] = dict(
        threads = 3,
        retry = 1,
        monitor_resource = True,
        monitor_time_step = 2,
        check_resource_before_run = True,
    )
    fastq_info_dict = parse_fastq_info(fastq_info_file)
    commands.update(fastqc_cmds(fastq_info_dict))
    trim_cmds = trimmomatic_cmds(fastq_info_dict)
    commands.update(trim_cmds)
    align_cmds = star_align_cmds(trimming_cmds=trim_cmds)
    commands.update(align_cmds)
    align_cmds.pop('AlignIndex')
    assembly_cmds = scallop_cmds(align_cmds=align_cmds)
    commands.update(assembly_cmds)
    quant_cmds = salmon_quant_cmds(trimming_cmds=trim_cmds)
    commands.update(quant_cmds)
    with open('pipeline_cmds.ini', 'w') as configfile:
        commands.write(configfile)



pipeline()