import os
import argparse
from pprint import pprint
import configparser
import shutil

from nestcmd import RunCommands, set_logger


class Basic(object):
    def __init__(self, workflow_arguments):
        self.logger = None
        workflow_arguments = self.arg_preprocess(workflow_arguments)
        # self.do_some_pre_judge(workflow_arguments)
        if workflow_arguments.arg_cfg:
            self.arg_pool = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
            self.arg_pool.optionxform = str
            self.arg_pool.read(workflow_arguments.arg_cfg, encoding='utf-8')
        self.project_dir = workflow_arguments.o
        self.workflow_arguments = workflow_arguments
        self.workflow = self.init_workflow_dict()

    def arg_preprocess(self, arguments):
        # get script path
        script_path = os.path.abspath(__file__)
        if os.path.islink(script_path):
            script_path = os.readlink(script_path)

        project_dir = arguments.o
        if arguments.only_show_steps or arguments.only_show_detail_steps:
            if arguments.pipeline_cfg is None and (not arguments.continue_run):
                # the project dir will be deleted after showing help
                project_dir = project_dir + '.tmp'

        if not os.path.exists(project_dir):
            os.mkdir(project_dir)
        project_dir = os.path.abspath(project_dir)
        arguments.o = project_dir
        log_file = os.path.join(project_dir, 'pipeline.log')
        self.logger = set_logger(log_file, 'pipeline')

        if not arguments.arg_cfg:
            if not arguments.pipeline_cfg:
                arg_file = os.path.join(os.path.dirname(script_path), 'arguments.ini')
                self.logger.warning("You are using unchanged configuration: {}".format(arg_file))
                arguments.arg_cfg = arg_file

        return arguments

    def do_some_pre_judge(self, arguments):
        if arguments.show_cmd_example:
            self.show_cmd_example(arguments.show_cmd_example)
            return True

        if arguments.list_cmd_names:
            self.list_cmd_names()
            return True

        if arguments.continue_run:
            if not arguments.pipeline_cfg:
                raise Exception('Existed pipeline_cfg must be provided !')
            self.run_existed_pipeline(steps=arguments.rerun_steps)
            return True
        else:
            if arguments.pipeline_cfg:
                self.logger.warning('You are re-running the whole existed pipeline')
                self.run_existed_pipeline()
                return True

        if not arguments.continue_run or arguments.pipeline_cfg is None:
            if not arguments.arg_cfg:
                raise Exception('-arg_cfg is needed!')
            if not os.path.exists(arguments.arg_cfg):
                raise Exception('arg_cfg file not exist')

    def init_workflow_dict(self):
        commands = configparser.ConfigParser()
        commands.optionxform = str
        commands['mode'] = dict(
            threads=self.workflow_arguments.threads,
            retry=self.workflow_arguments.retry,
            monitor_resource=not self.workflow_arguments.no_monitor_resource,
            monitor_time_step=self.workflow_arguments.monitor_time_step,
            check_resource_before_run=not self.workflow_arguments.no_check_resource_before_run,
        )
        return commands

    def cmd_dict(self, cmd, **kwargs):
        args = dict(
            retry=self.workflow_arguments.retry,
            monitor_resource=not self.workflow_arguments.no_monitor_resource,
            monitor_time_step=self.workflow_arguments.monitor_time_step,
            check_resource_before_run=not self.workflow_arguments.no_check_resource_before_run,
            cpu=1, mem=1024 * 1024 * 1024,
            cmd=cmd
        )
        if kwargs:
            args.update(kwargs)
        return args

    @staticmethod
    def mkdir(path):
        if not os.path.exists(path):
            os.mkdir(path)

    @staticmethod
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

    def run_existed_pipeline(self, steps=''):
        if self.workflow_arguments.pipeline_cfg is None or not os.path.exists(self.workflow_arguments.pipeline_cfg):
            raise Exception('Please provide valid pipeline.ini file')
        workflow = RunCommands(
            self.workflow_arguments.pipeline_cfg,
            outdir=self.project_dir, logger=self.logger
        )

        if self.workflow_arguments.only_show_steps:
            pprint('----Pipeline has the following main steps----')
            tmp_list = [x.split('_', 1)[0] for x in workflow.names()]
            main_steps = list()
            _ = [main_steps.append(x) for x in tmp_list if x not in main_steps]
            pprint(main_steps)
            return
        elif self.workflow_arguments.only_show_detail_steps:
            pprint('----Pipeline has the following steps----')
            pprint(workflow.names())
            return
        if self.workflow_arguments.continue_run:
            workflow.continue_run(steps=steps)
        else:
            workflow.parallel_run()

    def show_cmd_example(self, cmd_name):
        """
        :param cmd_name: cmd_generator中的函数名,也是arguments.ini中的section名
        :return: None
        """
        print('This function should be overrided')
        if not self.workflow_arguments.arg_cfg:
            raise Exception("please first input arg_cfg ")
        if cmd_name not in self.arg_pool:
            raise Exception('please provide valid cmd_name, refer --list_cmd_names')
        exec("print(cmdx.{}(**self.arg_pool['{}']))".format(cmd_name, cmd_name))

    def list_cmd_names(self):
        if not self.workflow_arguments.arg_cfg:
            raise Exception("please first input arg_cfg ")
        print(list(self.arg_pool.keys())[1:])

    def skip_some_steps(self):
        commands = self.workflow
        skip_steps = self.workflow_arguments.skip
        project_dir = self.project_dir
        if skip_steps:
            with open(os.path.join(project_dir, 'pipeline.ini'), 'w') as configfile:
                commands.write(configfile)
            workflow = RunCommands(os.path.join(project_dir, 'pipeline.ini'),
                                   outdir=project_dir, logger=self.logger)
            for each in skip_steps:
                skips = [x for x in commands if x == each or x.startswith(each + '_')]
                if not skips:
                    exit('Step {} was Not found, please refer --only_show_steps or --only_show_detail_steps'.format(
                        each))
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
                self.logger.warning("Warning: the following steps are also skipped for depending relationship")
                self.logger.warning(set(x.split('_', 1)[0] for x in total_deduced_skips))

    def run(self):
        self.skip_some_steps()
        commands = self.workflow
        arguments = self.workflow_arguments
        project_dir = self.project_dir
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
        workflow = RunCommands(os.path.join(project_dir, 'pipeline.ini'),
                               logger=self.logger,
                               outdir=project_dir,
                               timeout=arguments.wait_resource_time)
        workflow.parallel_run()


def basic_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-arg_cfg', required=False, help="config file containing all parameters info")
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
    return parser
