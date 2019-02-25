# coding=utf-8
import os
import configparser
import argparse
import shutil
from pprint import pprint

from cmd_generator import *
from nestcmd import RunCommands, set_logger

parser = argparse.ArgumentParser()
parser.add_argument('-arg_cfg', required=False, help="config file containing all parameters info")
parser.add_argument('-fastq_info', required=False,
                    help="第一列为样本名,分析结果使用名; 第二列为read1的fastq路径,如有多个,需要分号分隔; "
                         "第三列为可选, read2的fastq路径,单端测序时则无第三列")
parser.add_argument('-group', help="样本分组信息文件,至少两列,第一列样本名,第二列为分组名,其他列也是分组名")
parser.add_argument('-compare', help="比较信息文件,两列,第1列是对照组名,第2列是实验组名")
parser.add_argument('-o', default=os.path.join(os.getcwd(), 'Result'), help='分析目录或结果目录')
parser.add_argument('-skip', default=list(), nargs='+',
                    help='指定要跳过的步骤名, 空格分隔,程序会自动跳过依赖他们的步骤, --only_show_steps可查看步骤名; '
                         '注意: 如果跳过trim步骤, 则用原始数据; 如果跳过assembl或mergeTranscript, 则仅对参考基因定量')
parser.add_argument('--only_show_steps', default=False, action="store_true",
                    help="仅仅显示当前流程包含的主步骤, 且已经排除指定跳过的步骤; "
                         "你还可以通过--list_cmd_names查看当前流程包含哪些命令行")
parser.add_argument('--only_show_detail_steps', default=False, action="store_true",
                    help="仅仅显示当前流程包含的详细步骤, 且已经排除指定跳过的步骤")
parser.add_argument('--only_write_pipeline', default=False, action='store_true',
                    help="仅仅生成流程pipeline.ini")
parser.add_argument('-threads', default=5, type=int, help="允许并行的步骤数")
parser.add_argument('-retry', default=1, type=int,
                    help='某步骤运行失败后再尝试运行的次数, 默认1次. 如需对某一步设置不同的值, 可在运行流程前修改pipeline.ini')
parser.add_argument('--continue_run', default=False, action='store_true',
                    help='流程运行结束后, 从失败的步骤续跑, 记得要用-o指定之前的结果目录, 用-pipeline_cfg指定pipeline.ini; '
                         '如果顺便还想重新跑已经成功运行的步骤, 可通过-rerun_steps指定, 或者在状态表cmd_stat.txt中将其修改为failed即可')
parser.add_argument('-rerun_steps', default=list(), nargs='+',
                    help="使用--continue_run有效, 通过该参数指定重跑已经成功的步骤, 空格分隔, 这样做的可能原因可以是: 你重新设置了参数")
parser.add_argument('-pipeline_cfg', default=None,
                    help="已有的pipeline.ini, 续跑时必须需提供; 如提供该参数, 则此时无视 arg_cfg,fastq_info,group,cmp,skip 等参数")
parser.add_argument('--list_cmd_names', default=False, action='store_true', help="仅输出参数配置文件里的包含的cmd名称")
parser.add_argument('-show_cmd_example', help="提供一个cmd名称,输出该cmd的样例")
parser.add_argument('--no_monitor_resource', default=False, action='store_true',
                    help='是否监控每一步运行时的资源消耗, 如需对某一步设置不同的值, 可在运行流程前修改pipeline.ini')
parser.add_argument('--monitor_time_step', default=3, type=int,
                    help='监控资源时的时间间隔, 默认3秒, 如需对某一步设置不同的值, 可在运行流程前修改pipeline.ini')
parser.add_argument('-wait_resource_time', default=60, type=int,
                    help="等待资源的时间上限, 默认60秒, 等待时间超过这个时间时,资源不足时判定任务失败")
parser.add_argument('--no_check_resource_before_run', default=False, action='store_true',
                    help="指示运行某步骤前检测指定的资源是否足够, 如不足, 则该步骤失败; 如果设置该参数, 则运行前不检查资源. "
                         "如需对某一步设置不同的值,可运行前修改pipeline.ini. "
                         "如需更改指定的资源, 可在运行流程前修改pipeline.ini")
# get script path
script_path = os.path.abspath(__file__)
if os.path.islink(script_path):
    script_path = os.readlink(script_path)

arguments = parser.parse_args()
skip_steps = arguments.skip
arg_pool = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
arg_pool.optionxform = str

project_dir = arguments.o
if arguments.only_show_steps or arguments.only_show_detail_steps:
    if arguments.pipeline_cfg is None and (not arguments.continue_run):
        project_dir = project_dir + '.tmp'
if not os.path.exists(project_dir):
    os.mkdir(project_dir)
project_dir = os.path.abspath(project_dir)
# fastq_info_file = arguments.fastq_info
logger = set_logger(os.path.join(project_dir, 'workflow.log'))

if arguments.arg_cfg:
    arg_pool.read(arguments.arg_cfg, encoding='utf-8')
else:
    if not arguments.pipeline_cfg:
        arg_file = os.path.join(os.path.dirname(script_path), 'arguments.ini')
        logger.warning("You are using unchanged configuration: {}".format(arg_file))
        arguments.arg_cfg = arg_file
        arg_pool.read(arguments.arg_cfg, encoding='utf-8')


def cmd_dict(cmd, cpu=1, mem=200*1024*1024, retry=arguments.retry,
             monitor_resource=not arguments.no_monitor_resource,
             monitor_time_step=arguments.monitor_time_step,
             check_resource_before_run=not arguments.no_check_resource_before_run,
             **kwargs):
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
        if len(os.listdir(arg_pool['star_index']['genomeDir'])) <= 1:
            mkdir(arg_pool['star_index']['genomeDir'])
        else:
            logger.warning('STAR index existed, and skip this indexing step!')
            return commands
    else:
        mkdir(arg_pool['star_index']['genomeDir'])
    cmd = star_index(**arg_pool['star_index'])
    commands[step_name] = cmd_dict(
        cmd=cmd, cpu=1, mem=2*1024**3, retry=1,
        monitor_time_step=5)
    return commands


def star_align_cmds(trimming_cmds, index_cmd, step_name='Align'):
    commands = dict()
    outdir = os.path.join(project_dir, step_name)
    mkdir(outdir)
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
        args = dict(arg_pool['star_align'])
        args['readFilesIn'] = ','.join(trimmed_fq1_list) + ' ' + ','.join(trimmed_fq2_list)
        result_dir = os.path.join(outdir, sample)
        mkdir(result_dir)
        prefix = os.path.join(result_dir, sample + '.')
        args['outFileNamePrefix'] = prefix
        cmd = star_align(**args)
        commands[step_name+'_'+sample] = cmd_dict(
            cmd=cmd, mem=1024**3*1, cpu=2,
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
        args = dict(arg_pool['star_align'])
        if mode in 'PE':
            args['readFilesIn'] = ','.join(fq_list[0]) + ' ' + ','.join(fq_list[1])
        else:
            args['readFilesIn'] = ','.join(fq_list[0])
        result_dir = os.path.join(outdir, sample)
        mkdir(result_dir)
        prefix = os.path.join(result_dir, sample + '.')
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


def bam_index_cmds(align_cmds, step_name='IndexBam'):
    commands = dict()
    for step, cmd_info in align_cmds.items():
        sample = cmd_info['sample_name']
        args = dict(arg_pool['samtools_index'])
        args['bam'] = cmd_info['sorted_bam']
        cmd = samtools_index(**args)
        commands[step_name + '_' + sample] = cmd_dict(
            cmd=cmd, mem=1024 ** 3 * 2, cpu=2, monitor_time_step=5,
            depend=step,
            sample_name=sample,
            sorted_bam=args['bam']
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
        all_gtf=os.path.join(outdir, 'all.transcripts.gtf'),
        novel_transcripts=os.path.join(outdir, 'novel.transcripts.fa')
    )
    return commands


def salmon_index_cmd(step_name='QuantIndex', merge_transcript_cmd=None):
    commands = dict()
    if not merge_transcript_cmd:
        if os.path.exists(arg_pool['salmon_index']['index_prefix']):
            logger.warning('salmon index existed, and skip this indexing step!')
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
                    args['index'] = list(index_cmd.values())[0]['index_prefix']
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
                    args['index'] = list(index_cmd.values())[0]['index_prefix']
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
                    args['index'] = list(index_cmd.values())[0]['index_prefix']
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
                    args['index'] = list(index_cmd.values())[0]['index_prefix']
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


def merge_quant_cmd(quant_cmds, step_name='MergeQuant', quant_method='salmon', level='gene'):
    commands = dict()
    depend = list()
    step_name = step_name + level.capitalize()
    out_dir = os.path.join(project_dir, step_name)
    mkdir(out_dir)
    result_dir = ''
    for step, cmd_info in quant_cmds.items():
        result_dir = cmd_info['result_dir']
        depend.append(step)
    args = dict(arg_pool['abundance_estimates_to_matrix'])
    args['est_method'] = quant_method
    args['out_prefix'] = os.path.join(out_dir, level)
    if not level == 'gene':
        args['quant_result'] = result_dir + '/*/quant.sf'
        cmd = abundance_estimates_to_matrix(**args)
        commands[step_name] = cmd_dict(
            cmd=cmd, monitor_resource=False, check_resource_before_run=False,
            depend=','.join(depend), out_prefix=args['out_prefix'],
            transcript_tpm_matrix=args['out_prefix'] + '.isoform.TMM.EXPR.matrix',
            transcript_count_matrix=args['out_prefix'] + '.isoform.counts.matrix',
        )
    else:
        args['quant_result'] = result_dir + '/*/quant.genes.sf'
        cmd = abundance_estimates_to_matrix(**args)
        commands[step_name] = cmd_dict(
            cmd=cmd, monitor_resource=False, check_resource_before_run=False,
            depend=','.join(depend), out_prefix=args['out_prefix'],
            gene_tpm_matrix=args['out_prefix'] + '.isoform.TMM.EXPR.matrix',
            gene_count_matrix=args['out_prefix'] + '.isoform.counts.matrix',
        )
    return commands


def diff_exp_cmd(merge_cmd, step_name='Diff', level='gene'):
    commands = dict()
    step_name = step_name + level.capitalize()
    out_dir = os.path.join(project_dir, step_name)
    mkdir(out_dir)
    # gene diff exp
    if level == 'gene':
        depend = [x for x in merge_cmd.keys() if x.endswith('Gene')][0]
        depend_info = merge_cmd[depend]
        args = dict(arg_pool['diff_exp'])
        args['result_dir'] = out_dir
        args['count_matrix'] = depend_info['gene_count_matrix']
        args['exp_matrix'] = depend_info['gene_tpm_matrix']
    else:
        depend = [x for x in merge_cmd.keys() if x.endswith('Transcript')][0]
        depend_info = merge_cmd[depend]
        args = dict(arg_pool['diff_exp'])
        args['result_dir'] = out_dir
        args['count_matrix'] = depend_info['transcript_count_matrix']
        args['exp_matrix'] = depend_info['transcript_tpm_matrix']
    args['group_info'] = arguments.group
    args['comparison_info'] = arguments.compare
    cmd = diff_exp(**args)
    commands[step_name] = cmd_dict(
        cmd=cmd,
        cpu=int(args['threads']) + 1,
        depend=depend,
        result_dir=args['result_dir']
    )
    return commands


def go_enrich_cmds(diffexp_cmd, step_name='GoEnrich', level='gene'):
    commands = dict()
    out_dir = os.path.join(project_dir, step_name+level.capitalize())
    mkdir(out_dir)
    depend = list(diffexp_cmd.keys())[0]
    depend_info = list(diffexp_cmd.values())[0]
    diff_list_files = list()
    with open(arguments.compare) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            ctrl, test = line.strip().split()
            diff_list_files.append(os.path.join(depend_info['result_dir'], '{}_vs_{}.*.DE.list'.format(ctrl, test)))
    for each in diff_list_files:
        args = dict(arg_pool['goatools'])
        args['study'] = each
        cmp_name = os.path.basename(each).split('.', 1)[0]
        args['outfile'] = os.path.join(out_dir, str(cmp_name)+'.goea.xls')
        if level == 'gene':
            cmd = goatools(**args)
        else:
            args['population'] = args['trans_population']
            args['association'] = args['trans_association']
            cmd = goatools(**args)
        commands[step_name+level.capitalize()+'_'+str(cmp_name)] = cmd_dict(
            cmd=cmd,
            cpu=1,
            depend=depend,
            result_dir=out_dir
        )
    return commands


def kegg_enrich_cmds(diffexp_cmd, step_name='KeggEnrich', level='gene'):
    commands = dict()
    out_dir = os.path.join(project_dir, step_name+level.capitalize())
    mkdir(out_dir)
    depend = list(diffexp_cmd.keys())[0]
    depend_info = list(diffexp_cmd.values())[0]
    diff_list_files = list()
    with open(arguments.compare) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            ctrl, test = line.strip().split()
            diff_list_files.append(os.path.join(depend_info['result_dir'], '{}_vs_{}.*.DE.list'.format(ctrl, test)))
    for each in diff_list_files:
        args = dict(arg_pool['kegg_enrich'])
        args['deg'] = each
        args['outdir'] = out_dir
        if level == 'gene':
            cmd = kegg_enrich(**args)
        else:
            args['g2k'] = args['t2k']
            args['g2p'] = args['t2p']
            cmd = kegg_enrich(**args)
        cmp_name = os.path.basename(each).split('.', 1)[0]
        commands[step_name+level.capitalize()+'_'+str(cmp_name)] = cmd_dict(
            cmd=cmd,
            cpu=1,
            depend=depend,
            result_dir=out_dir
        )
    return commands


def exp_analysis_cmd(merge_cmd, step_name='ExpAnalysis', level='gene'):
    commands = dict()
    step_name = step_name + level.capitalize()
    out_dir = os.path.join(project_dir, step_name)
    mkdir(out_dir)
    args = dict(arg_pool['exp_analysis'])
    depend = list(merge_cmd.keys())[0]
    depend_info = list(merge_cmd.values())[0]
    if level == 'gene':
        args['out_prefix'] = os.path.join(out_dir, 'gene.tpm')
        args['matrix'] = depend_info['gene_tpm_matrix']
    else:
        args['out_prefix'] = os.path.join(out_dir, 'transcript.tpm')
        args['matrix'] = depend_info['transcript_tpm_matrix']
    cmd = exp_analysis(**args)
    commands[step_name] = cmd_dict(
        cmd=cmd,
        cpu=1,
        depend=depend,
        result_dir=out_dir
    )
    return commands


def gene_body_coverage_cmds(index_bam_cmds, step_name='GeneBodyCoverage'):
    commands = dict()
    out_dir = os.path.join(project_dir, step_name)
    mkdir(out_dir)
    args = dict(arg_pool['gene_body_coverage'])
    for step, cmd_info in index_bam_cmds.items():
        sample = cmd_info['sample_name']
        args['bam'] = cmd_info['sorted_bam']
        args['out_prefix'] = os.path.join(out_dir, sample)
        cmd = gene_body_coverage(**args)
        commands[step_name + '_' + sample] = cmd_dict(
            cmd=cmd,
            depend=step,
            out_prefix=args['out_prefix']
        )
    return commands


def inner_distance_cmds(index_bam_cmds, step_name='InnerDistance'):
    commands = dict()
    out_dir = os.path.join(project_dir, step_name)
    mkdir(out_dir)
    args = dict(arg_pool['inner_distance'])
    for step, cmd_info in index_bam_cmds.items():
        sample = cmd_info['sample_name']
        args['bam'] = cmd_info['sorted_bam']
        args['out_prefix'] = os.path.join(out_dir, sample)
        cmd = inner_distance(**args)
        commands[step_name + '_' + sample] = cmd_dict(
            cmd=cmd,
            depend=step,
            out_prefix=args['out_prefix']
        )
    return commands


def read_distribution_cmds(index_bam_cmds, step_name='ReadDistribution'):
    commands = dict()
    out_dir = os.path.join(project_dir, step_name)
    mkdir(out_dir)
    args = dict(arg_pool['read_distribution'])
    for step, cmd_info in index_bam_cmds.items():
        sample = cmd_info['sample_name']
        args['bam'] = cmd_info['sorted_bam']
        args['outfile'] = os.path.join(out_dir, sample+'.read_distribution.txt')
        cmd = read_distribution(**args)
        commands[step_name + '_' + sample] = cmd_dict(
            cmd=cmd,
            depend=step,
            outfile=args['outfile']
        )
    return commands


def read_duplication_cmds(index_bam_cmds, step_name='ReadDuplication'):
    commands = dict()
    out_dir = os.path.join(project_dir, step_name)
    mkdir(out_dir)
    args = dict(arg_pool['read_duplication'])
    for step, cmd_info in index_bam_cmds.items():
        sample = cmd_info['sample_name']
        args['bam'] = cmd_info['sorted_bam']
        args['out_prefix'] = os.path.join(out_dir, sample)
        cmd = read_duplication(**args)
        commands[step_name + '_' + sample] = cmd_dict(
            cmd=cmd,
            depend=step,
            out_prefix=args['out_prefix']
        )
    return commands


def rna_fragment_size_cmds(index_bam_cmds, step_name='FragmentSize'):
    commands = dict()
    out_dir = os.path.join(project_dir, step_name)
    mkdir(out_dir)
    args = dict(arg_pool['rna_fragment_size'])
    for step, cmd_info in index_bam_cmds.items():
        sample = cmd_info['sample_name']
        args['bam'] = cmd_info['sorted_bam']
        args['outfile'] = os.path.join(out_dir, sample+'.fragment_size.txt')
        cmd = rna_fragment_size(**args)
        commands[step_name + '_' + sample] = cmd_dict(
            cmd=cmd,
            depend=step,
            outfile=args['outfile']
        )
    return commands


# def rpkm_saturation_cmds(index_bam_cmds, step_name='RPKMSaturation'):
#     commands = dict()
#     out_dir = os.path.join(project_dir, step_name)
#     mkdir(out_dir)
#     args = dict(arg_pool['rpkm_saturation'])
#     for step, cmd_info in index_bam_cmds.items():
#         sample = cmd_info['sample_name']
#         args['bam'] = cmd_info['sorted_bam']
#         args['out_prefix'] = os.path.join(out_dir, sample)
#         cmd = rpkm_saturation(**args)
#         commands[step_name + '_' + sample] = cmd_dict(
#             cmd=cmd,
#             depend=step,
#             out_prefix=args['out_prefix']
#         )
#     return commands


def tpm_saturation_cmds(index_bam_cmds, step_name='TPMSaturation'):
    commands = dict()
    out_dir = os.path.join(project_dir, step_name)
    mkdir(out_dir)
    args = dict(arg_pool['tpm_saturation'])
    for step, cmd_info in index_bam_cmds.items():
        sample = cmd_info['sample_name']
        args['bam'] = cmd_info['sorted_bam']
        args['outdir'] = out_dir
        cmd = tpm_saturation(**args)
        commands[step_name + '_' + sample] = cmd_dict(
            cmd=cmd,
            depend=step,
            exp_matrix='{}.tpm.xls'.format(sample)
        )
    return commands


def get_alignment_summary_cmds(index_bam_cmds, step_name='AlignmentSummary'):
    commands = dict()
    out_dir = os.path.join(project_dir, step_name)
    mkdir(out_dir)
    args = dict(arg_pool['get_alignment_summary'])
    for step, cmd_info in index_bam_cmds.items():
        sample = cmd_info['sample_name']
        args['bam'] = cmd_info['sorted_bam']
        args['outdir'] = out_dir
        cmd = get_alignment_summary(**args)
        commands[step_name + '_' + sample] = cmd_dict(
            cmd=cmd,
            depend=step,
            alignment_summary='{}.alignment_summary.json'.format(sample)

        )
    return commands


def CollectAlignmentSummaryMetrics_cmds(index_bam_cmds, step_name='CollectAlignmentSummaryMetrics'):
    commands = dict()
    out_dir = os.path.join(project_dir, step_name)
    mkdir(out_dir)
    args = dict(arg_pool['CollectAlignmentSummaryMetrics'])
    for step, cmd_info in index_bam_cmds.items():
        sample = cmd_info['sample_name']
        args['bam'] = cmd_info['sorted_bam']
        args['outfile'] = os.path.join(out_dir, '{}.{}.xls'.format(sample, step_name))
        cmd = CollectAlignmentSummaryMetrics(**args)
        commands[step_name + '_' + sample] = cmd_dict(
            cmd=cmd,
            depend=step,
            metrics=args['outfile']

        )
    return commands


def CollectInsertSizeMetrics_cmds(index_bam_cmds, step_name='CollectInsertSizeMetrics'):
    commands = dict()
    out_dir = os.path.join(project_dir, step_name)
    mkdir(out_dir)
    args = dict(arg_pool['CollectInsertSizeMetrics'])
    for step, cmd_info in index_bam_cmds.items():
        sample = cmd_info['sample_name']
        args['bam'] = cmd_info['sorted_bam']
        args['outfile'] = os.path.join(out_dir, '{}.{}.xls'.format(sample, step_name))
        args['outimage'] = os.path.join(out_dir, '{}.{}.pdf'.format(sample, step_name))
        cmd = CollectInsertSizeMetrics(**args)
        commands[step_name + '_' + sample] = cmd_dict(
            cmd=cmd,
            depend=step,
            metrics=args['outfile']

        )
    return commands


def CollectTargetedPcrMetrics_cmds(index_bam_cmds, step_name='CollectTargetedPcrMetrics'):
    commands = dict()
    out_dir = os.path.join(project_dir, step_name)
    mkdir(out_dir)
    args = dict(arg_pool['CollectTargetedPcrMetrics'])
    for step, cmd_info in index_bam_cmds.items():
        sample = cmd_info['sample_name']
        args['bam'] = cmd_info['sorted_bam']
        args['outfile'] = os.path.join(out_dir, '{}.{}.xls'.format(sample, step_name))
        cmd = CollectTargetedPcrMetrics(**args)
        commands[step_name + '_' + sample] = cmd_dict(
            cmd=cmd,
            depend=step,
            metrics=args['outfile']

        )
    return commands


def CollectRnaSeqMetrics_cmds(index_bam_cmds, step_name='CollectRnaSeqMetrics'):
    commands = dict()
    out_dir = os.path.join(project_dir, step_name)
    mkdir(out_dir)
    args = dict(arg_pool['CollectRnaSeqMetrics'])
    for step, cmd_info in index_bam_cmds.items():
        sample = cmd_info['sample_name']
        args['bam'] = cmd_info['sorted_bam']
        args['outfile'] = os.path.join(out_dir, '{}.{}.xls'.format(sample, step_name))
        args['outimage'] = os.path.join(out_dir, '{}.{}.pdf'.format(sample, step_name))
        cmd = CollectRnaSeqMetrics(**args)
        commands[step_name + '_' + sample] = cmd_dict(
            cmd=cmd,
            depend=step,
            metrics=args['outfile']

        )
    return commands


def run_existed_pipeline(steps=''):
    if arguments.pipeline_cfg is None or not os.path.exists(arguments.pipeline_cfg):
        raise Exception('Please provide valid pipeline.ini file')
    workflow = RunCommands(arguments.pipeline_cfg, outdir=project_dir, logger=logger)

    if arguments.only_show_steps:
        pprint('----Pipeline has the following main steps----')
        tmp_list = [x.split('_', 1)[0] for x in workflow.names()]
        main_steps = list()
        _ = [main_steps.append(x) for x in tmp_list if x not in main_steps]
        pprint(main_steps)
        return
    elif arguments.only_show_detail_steps:
        pprint('----Pipeline has the following steps----')
        pprint(workflow.names())
        return
    if arguments.continue_run:
        workflow.continue_run(steps=steps)
    else:
        workflow.parallel_run()


def show_cmd_example(cmd_name):
    """
    :param cmd_name: cmd_generator中的函数名,也是arguments.ini中的section名
    :return: None
    """
    if not arguments.arg_cfg:
        raise Exception("please first input arg_cfg ")
    if cmd_name not in arg_pool:
        raise Exception('please provide valid cmd_name, refer --list_cmd_names')
    exec("pprint({}(**arg_pool['{}']))".format(cmd_name, cmd_name))


def list_cmd_names():
    if not arguments.arg_cfg:
        raise Exception("please first input arg_cfg ")
    pprint(list(arg_pool.keys())[1:])


def pipeline():
    """
    为了能正常跳过一些步骤,步骤名即step_name不能包含'_',且保证最后生成的步骤名不能有重复。
    """
    if arguments.only_show_steps:
        if arguments.fastq_info is None:
            test_data_dir = os.path.join(os.path.dirname(script_path), 'testdata')
            if os.path.exists(test_data_dir):
                arguments.fastq_info = os.path.join(test_data_dir, 'fastq_info.txt')
                arguments.compare = os.path.join(test_data_dir, 'compare')
                arguments.group = os.path.join(test_data_dir, 'group')

    if arguments.show_cmd_example:
        show_cmd_example(arguments.show_cmd_example)
        return
    if arguments.list_cmd_names:
        list_cmd_names()
        return

    if arguments.continue_run:
        if not arguments.pipeline_cfg:
            raise Exception('Existed pipeline_cfg must be provided !')
        run_existed_pipeline(steps=arguments.rerun_steps)
        return
    else:
        if arguments.pipeline_cfg:
            logger.warning('You are re-running the whole existed pipeline')
            run_existed_pipeline()
            return

    if not arguments.continue_run or arguments.pipeline_cfg is None:
        if not arguments.arg_cfg:
            raise Exception('-arg_cfg is needed!')
        if not os.path.exists(arguments.arg_cfg):
            raise Exception('arg_cfg file not exist')
        if not arguments.fastq_info:
            if arguments.only_show_detail_steps:
                shutil.rmtree(project_dir)
            raise Exception('-fastq_info is needed!')

    commands = configparser.ConfigParser()
    commands.optionxform = str
    commands['mode'] = dict(
        threads=arguments.threads,
        retry=arguments.retry,
        monitor_resource=not arguments.no_monitor_resource,
        monitor_time_step=arguments.monitor_time_step,
        check_resource_before_run=not arguments.no_check_resource_before_run,
    )

    # fastqc and trimmomatic
    fastq_info_dict = parse_fastq_info(arguments.fastq_info)
    commands.update(fastqc_raw_data_cmds(fastq_info_dict, step_name='RawDataQC'))
    trim_cmds = trimmomatic_cmds(fastq_info_dict, step_name='Trim')
    commands.update(trim_cmds)
    fastqc_cmds = fastqc_trimmed_data_cmds(trimming_cmds=trim_cmds, step_name='TrimmedDataQC')
    commands.update(fastqc_cmds)

    # align and bam index
    star_indexing = star_index_cmd(step_name='AlignIndex')
    commands.update(star_indexing)
    if list(trim_cmds.keys())[0].split('_', 1)[0] in skip_steps:
        align_cmds = star_align_with_rawdata_cmds(fastq_info_dict, index_cmd=star_indexing, step_name='Align')
    else:
        align_cmds = star_align_cmds(trimming_cmds=trim_cmds, index_cmd=star_indexing, step_name='Align')
    commands.update(align_cmds)
    bam_indexing_cmds = bam_index_cmds(align_cmds, step_name='IndexBam')
    commands.update(bam_indexing_cmds)

    # run some RseQC cmds
    gbc_cmds = gene_body_coverage_cmds(bam_indexing_cmds, step_name='GeneBodyCoverage')
    commands.update(gbc_cmds)
    inner_dist_cmds = inner_distance_cmds(bam_indexing_cmds, step_name='InnerDistance')
    commands.update(inner_dist_cmds)
    rd_cmds = read_distribution_cmds(bam_indexing_cmds, step_name='ReadDistribution')
    commands.update(rd_cmds)
    rdup_cmds = read_duplication_cmds(bam_indexing_cmds, step_name='ReadDuplication')
    commands.update(rdup_cmds)
    frag_size_cmds = rna_fragment_size_cmds(bam_indexing_cmds, step_name='FragmentSize')
    commands.update(frag_size_cmds)

    # run TPM saturation
    saturation_cmds = tpm_saturation_cmds(bam_indexing_cmds, step_name='TPMSaturation')
    commands.update(saturation_cmds)

    # 根据star比对结果log文件和使用bedtools intersect 结合 samtools flagstat 统计比对结果
    alignment_summary_cmds = get_alignment_summary_cmds(bam_indexing_cmds, step_name='AlignmentSummary')
    commands.update(alignment_summary_cmds)

    # run some picard tools 获得大量的比对情况的统计结果，包括insert_size 和 gene_body_coverage
    collect_alignment_summary_cmds = CollectAlignmentSummaryMetrics_cmds(bam_indexing_cmds)
    commands.update(collect_alignment_summary_cmds)
    collect_insert_size_cmds = CollectInsertSizeMetrics_cmds(bam_indexing_cmds)
    commands.update(collect_insert_size_cmds)
    collect_target_info_cmds = CollectTargetedPcrMetrics_cmds(bam_indexing_cmds)
    commands.update(collect_target_info_cmds)
    collect_rna_metric_cmds = CollectRnaSeqMetrics_cmds(bam_indexing_cmds)
    commands.update(collect_rna_metric_cmds)

    # assemble and merge
    assembly_cmds = scallop_cmds(align_cmds=align_cmds, step_name='Assembly')
    commands.update(assembly_cmds)
    merge_trans_cmd = merge_scallop_transcripts_cmd(assembly_cmds, step_name='MergeTranscript')
    commands.update(merge_trans_cmd)

    # quant with salmon and merge result
    if list(merge_trans_cmd.keys())[0] in skip_steps or list(assembly_cmds.keys())[0].split('_', 1)[0] in skip_steps:
        salmon_indexing = salmon_index_cmd(merge_transcript_cmd=None, step_name='QuantIndex')
    else:
        salmon_indexing = salmon_index_cmd(merge_transcript_cmd=merge_trans_cmd, step_name='QuantIndex')
    commands.update(salmon_indexing)
    if list(trim_cmds.keys())[0].split('_', 1)[0] not in skip_steps:
        quant_cmds = salmon_quant_with_clean_data_cmds(trim_cmds, index_cmd=salmon_indexing, step_name='Quant')
    else:
        quant_cmds = salmon_quant_with_raw_data_cmds(fastq_info_dict, index_cmd=salmon_indexing, step_name='Quant')
    commands.update(quant_cmds)
    merge_gene_exp_cmd = merge_quant_cmd(quant_cmds=quant_cmds, step_name="MergeQuant")
    commands.update(merge_gene_exp_cmd)

    # expression analysis
    gene_exp_cluster_cmd = exp_analysis_cmd(merge_gene_exp_cmd, step_name='ExpAnalysis', level='gene')
    commands.update(gene_exp_cluster_cmd)
    merge_trans_exp_cmd = merge_quant_cmd(quant_cmds=quant_cmds, level='transcript', step_name='MergeQuant')
    commands.update(merge_trans_exp_cmd)
    trans_exp_cluster_cmd = exp_analysis_cmd(merge_trans_exp_cmd, step_name='ExpAnalysis', level='transcript')
    commands.update(trans_exp_cluster_cmd)

    # diff and enrich
    if arguments.group and arguments.compare:
        diff_gene_cmd = diff_exp_cmd(merge_gene_exp_cmd, level='gene', step_name='Diff')
        commands.update(diff_gene_cmd)
        gene_go_enrich_cmds = go_enrich_cmds(diff_gene_cmd, level='gene', step_name='GoEnrich')
        commands.update(gene_go_enrich_cmds)
        gene_kegg_enrich_cmds = kegg_enrich_cmds(diff_gene_cmd, level='gene', step_name='KeggEnrich')
        commands.update(gene_kegg_enrich_cmds)
        diff_trans_cmd = diff_exp_cmd(merge_trans_exp_cmd, level='transcript', step_name='Diff')
        commands.update(diff_trans_cmd)
        trans_go_enrich_cmds = go_enrich_cmds(diff_trans_cmd, level='transcript', step_name='GoEnrich')
        commands.update(trans_go_enrich_cmds)
        trans_kegg_enrich_cmds = kegg_enrich_cmds(diff_trans_cmd, level='transcript', step_name='KeggEnrich')
        commands.update(trans_kegg_enrich_cmds)

    # -----------skip some steps--------------
    if skip_steps:
        with open(os.path.join(project_dir, 'pipeline.ini'), 'w') as configfile:
            commands.write(configfile)
        workflow = RunCommands(os.path.join(project_dir, 'pipeline_raw.ini'), outdir=project_dir)
        for each in skip_steps:
            skips = [x for x in commands if x == each or x.startswith(each+'_')]
            if not skips:
                exit('Step {} was Not found, please refer --only_show_steps or --only_show_detail_steps'.format(each))
            else:
                _ = [commands.pop(x) for x in skips]
        # skip the step whose dependencies are not all included in the commands
        total_deduced_skips = list()
        while True:
            to_be_skip = list()
            for step in commands.keys():
                depends = workflow.get_dependency(step)
                if set(depends) - set(commands.keys()):
                    # print('Skip {} for at least one of its dependencies were skipped'.format(step))
                    to_be_skip.append(step)
                    total_deduced_skips.append(step)
            _ = [commands.pop(x) for x in to_be_skip]
            if all(len(set(workflow.get_dependency(x)) - set(commands.keys())) == 0 for x in commands):
                break
        if total_deduced_skips:
            logger.warning("Warning: the following steps are also skipped for depending relationship")
            logger.warning(set(x.split('_', 1)[0] for x in total_deduced_skips))

    # ----only show steps----
    tmp_list = [x.split('_', 1)[0] for x in commands.keys()]
    main_steps = list()
    _ = [main_steps.append(x) for x in tmp_list if x not in main_steps]
    if arguments.only_show_steps:
        pprint('----Pipeline has the following main steps----')
        pprint(main_steps[2:])
        shutil.rmtree(project_dir)
        return
    elif arguments.only_show_detail_steps:
        pprint('----Pipeline has the following steps----')
        pprint(list(commands.keys())[2:])
        shutil.rmtree(project_dir)
        return

    # ---------write pipeline cmds--------------------
    with open(os.path.join(project_dir, 'pipeline.ini'), 'w') as configfile:
        commands.write(configfile)
    if arguments.only_write_pipeline:
        return
    # ----------------run-----------------
    workflow = RunCommands(os.path.join(project_dir, 'pipeline.ini'), logger=logger,
                           outdir=project_dir, timeout=arguments.wait_resource_time)
    workflow.parallel_run()


if __name__ == '__main__':
    pipeline()
