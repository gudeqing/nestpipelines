# coding=utf-8
import time
import os
import configparser
import psutil
import queue
import logging
from subprocess import PIPE
import threading
import matplotlib
matplotlib.use('agg')
from threading import Timer, Lock
import weakref
import atexit
import signal
import pygraphviz as pgv
__author__ = 'gdq and dp'


PROCESS_SET = weakref.WeakKeyDictionary()


@atexit.register
def _kill_processes_when_exit():
    print("....Ending....")
    for proc, cmd_name in PROCESS_SET.items():
        if psutil.pid_exists(proc.pid):
            print('Shutting down running tasks {}:{}'.format(proc.pid, cmd_name))
            proc.kill()


def shutdownFunction(signum, frame):
    print('Killing main process, thus its derived processes will also be killed')
    exit(0)

# kill signal will be captured
signal.signal(signal.SIGTERM, shutdownFunction)
signal.signal(signal.SIGINT, shutdownFunction)


def set_logger(name='workflow.log', logger_id='x'):
    logger = logging.getLogger(logger_id)
    logger.propagate = False
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(name, mode='w+')
    fh.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setLevel(logging.WARNING)
    # fmt = '%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
    fmt = '%(asctime)s: %(message)s'
    format_str = logging.Formatter(fmt)  # 设置日志格式
    fh.setFormatter(format_str)
    logger.addHandler(sh)
    logger.addHandler(fh)
    return logger


class Command(object):
    def __init__(self, cmd, name, timeout=604800, outdir=os.getcwd(),
                 monitor_resource=True, monitor_time_step=2, logger=None, **kwargs):
        self.name = name
        self.cmd = cmd
        self.proc = None
        self.stdout = None
        self.stderr = None
        self.timeout = int(timeout)
        self.used_time = 0
        self.max_mem = 0
        self.max_cpu = 0
        self.monitor = monitor_resource
        self.monitor_time_step = int(monitor_time_step)
        self.outdir = outdir
        if not logger:
            self.logger = set_logger(name=os.path.join(self.outdir, 'command.log'))
        else:
            self.logger = logger

    def _monitor_resource(self):
        while self.proc.is_running():
            try:
                # if os.name == 'posix':
                #     cpu_num = self.proc.cpu_num()
                # elif os.name == 'nt':
                #     cpu_num = psutil.cpu_count()
                # else:
                #     cpu_num = 0
                cpu_percent = self.proc.cpu_percent(self.monitor_time_step)
                used_cpu = round(cpu_percent*0.01, 4)
                if used_cpu > self.max_cpu:
                    self.max_cpu = used_cpu
                memory_obj = self.proc.memory_full_info()
                # memory = (memory_obj.vms - memory_obj.shared)/1024/1024
                # memory = round(memory_obj.vms/1024/1024, 4)
                memory = round(memory_obj.uss/1024/1024, 4)
                if memory > self.max_mem:
                    self.max_mem = memory
            except Exception as e:
                # print('Failed to capture cpu/mem info for: ', e)
                break

    def run(self):
        start_time = time.time()
        self.logger.warning("RunStep: {}".format(self.name))
        self.logger.info("RunCmd: {}".format(self.cmd))
        self.proc = psutil.Popen(self.cmd, shell=True, stderr=PIPE, stdout=PIPE)
        PROCESS_SET[self.proc] = self.name
        if self.monitor:
            thread = threading.Thread(target=self._monitor_resource, daemon=True)
            thread.start()
        timer = Timer(self.timeout, self.proc.kill)
        try:
            timer.start()
            self.stdout, self.stderr = self.proc.communicate()
            if self.monitor:
                thread.join()
        finally:
            timer.cancel()
        self._write_log()
        end_time = time.time()
        self.used_time = round(end_time - start_time, 4)

    def _write_log(self):
        log_dir = os.path.join(self.outdir, 'logs')
        if not os.path.exists(log_dir):
            os.mkdir(log_dir)
        prefix = os.path.join(self.outdir, 'logs', self.name+'.'+str(self.proc.pid))
        if self.stderr:
            with open(prefix+'.stderr.txt', 'wb') as f:
                f.write(self.stderr)
        if self.stdout:
            with open(prefix+'.stdout.txt', 'wb') as f:
                f.write(self.stdout)
        if self.max_cpu or self.max_mem:
            with open(prefix+'.resource.txt', 'w') as f:
                f.write('max_cpu: {}\n'.format(self.max_cpu))
                f.write('max_mem: {}M\n'.format(round(self.max_mem, 4)))


class CommandNetwork(object):
    def __init__(self, cmd_config):
        self.parser = configparser.ConfigParser()
        self.parser.read(cmd_config, encoding='utf-8')

    def names(self):
        sections = self.parser.sections()
        # mode section is not cmd
        sections.pop(sections.index('mode'))
        return sections

    def orphans(self):
        independent_cmds = list()
        for name in self.names():
            if 'depend' not in self.parser[name]:
                independent_cmds.append(name)
            else:
                depend = self.parser[name]['depend'].strip()
                if not depend:
                    independent_cmds.append(name)
        return independent_cmds

    def get_dependency(self, name):
        if 'depend' not in self.parser[name]:
            return []
        else:
            depend = self.parser[name]['depend'].strip()
            if not depend:
                return []
            else:
                return [x.strip() for x in depend.split(',')]

    def get_cmd_description_dict(self, name):
        tmp_dict = dict(self.parser[name])
        tmp_dict['name'] = name
        if 'cpu' not in tmp_dict:
            tmp_dict['cpu'] = 0
        if 'mem' not in tmp_dict:
            tmp_dict['mem'] = 0
        if 'depend' not in tmp_dict:
            tmp_dict['depend'] = None
        if 'retry' not in tmp_dict:
            tmp_dict['retry'] = self.parser.getint('mode', 'retry')
        else:
            tmp_dict['retry'] = self.parser.getint(name, 'retry')
        if 'monitor_resource' not in tmp_dict:
            tmp_dict['monitor_resource'] = self.parser.getboolean('mode', 'monitor_resource')
        else:
            tmp_dict['monitor_resource'] = self.parser.getboolean(name, 'monitor_resource')
        if 'timeout' not in tmp_dict:
            tmp_dict['timeout'] = 3600*24*7
        else:
            tmp_dict['timeout'] = self.parser.getint(name, 'timeout')
        if 'monitor_time_step' not in tmp_dict:
            tmp_dict['monitor_time_step'] = self.parser.getint('mode', 'monitor_time_step')
        else:
            tmp_dict['monitor_time_step'] = self.parser.getint(name, 'monitor_time_step')
        if 'check_resource_before_run' not in tmp_dict:
            tmp_dict['check_resource_before_run'] = self.parser.getboolean('mode', 'check_resource_before_run')
        else:
            tmp_dict['check_resource_before_run'] = self.parser.getboolean(name, 'check_resource_before_run')
        return tmp_dict


class CheckResource(object):
    @staticmethod
    def available_mem():
        return psutil.virtual_memory().free

    @staticmethod
    def available_cpu():
        total = psutil.cpu_count()
        return int(total - total*psutil.cpu_percent()*0.01)

    def is_enough(self, cpu, mem, timeout=10):
        start_time = time.time()
        enough_num = 0
        while True:
            if float(cpu) <= self.available_cpu() \
                    and float(mem) <= self.available_mem():
                enough_num += 1
                if enough_num >= 3:
                    return True
                if enough_num >= 1 and timeout <= 10:
                    return True
            if time.time() - start_time >= timeout:
                return False
            time.sleep(3)


class StateGraph(object):
    def __init__(self, state):
        self.state = state
        self.graph = pgv.AGraph(directed=True, rankdir='LR')

    def _add_nodes(self):
        for node, cmd_info in self.state.items():
            status = cmd_info['state']
            node_detail = node.split('_', 1)
            if status == 'success':
                color = '#7FFF00'
            elif status == 'failed':
                color = '#FFD700'
            elif status == 'running':
                color = '#9F79EE'
            elif status == 'queueing':
                color = '#87CEFF'
            else:
                color = '#A8A8A8'
            used_time = cmd_info['used_time']
            if isinstance(used_time, str):
                if used_time == 'unknown':
                    pass
                else:
                    try:
                        float(used_time)
                        node_detail.append(used_time+'s')
                    except ValueError:
                        node_detail.append(used_time)
            elif float(used_time) <= 0:
                pass
            else:
                node_detail.append(str(used_time) + 's')
            self.graph.add_node(
                node,
                # 谷歌浏览器可以正常显示tooltip
                tooltip=cmd_info['cmd'].replace(' ', '\n'),
                shape="box",
                style="rounded, filled",
                fillcolor=color,
                color="mediumseagreen",
                label='\n'.join(node_detail)
            )

    def _add_edges(self):
        for target in self.state:
            sources = self.state[target]['depend'].strip()
            if sources:
                sources = sources.split(',')
                edges = zip(sources, [target]*len(sources))
                if self.state[target]['state'] == 'success':
                    color = 'green'
                elif self.state[target]['state'] == 'running':
                    color = '#836FFF'
                else:
                    color = '#4D4D4D'
                self.graph.add_edges_from(edges, color=color)
            else:
                self.graph.add_edge('Input', target, color='green')

    def draw(self, img_file='state.svg'):
        self._add_nodes()
        self._add_edges()
        img_fmt = os.path.splitext(img_file)[1][1:]
        self.graph.draw(path=img_file, format=img_fmt, prog='dot')


class RunCommands(CommandNetwork):
    __LOCK__ = Lock()

    def __init__(self, cmd_config, outdir=os.getcwd(), timeout=10, logger=None):
        super().__init__(cmd_config)
        self.ever_queued = set()
        self.queue = self.__init_queue()
        self.state = self.__init_state()
        self.outdir = outdir
        self._draw_state()
        self.timeout = timeout
        if not logger:
            self.logger = set_logger(name=os.path.join(self.outdir, 'workflow.log'))
        else:
            self.logger = logger

    def __init_queue(self):
        cmd_pool = queue.Queue()
        for each in self.orphans():
            cmd_pool.put(each)
            self.ever_queued.add(each)
        return cmd_pool

    def __init_state(self):
        state_dict = dict()
        for name in self.names():
            state_dict[name] = dict()
            fields = ['state', 'used_time', 'mem', 'cpu', 'pid', 'depend', 'cmd']
            for each in fields:
                state_dict[name][each] = 'unknown'
            state_dict[name]['cmd'] = self.parser[name]['cmd']
            state_dict[name]['depend'] = ','.join(self.get_dependency(name))
        return state_dict

    def _update_queue(self):
        success = set(x for x in self.state if self.state[x]['state'] == 'success')
        failed = set(x for x in self.state if self.state[x]['state'] == 'failed')
        waiting = set(self.names()) - self.ever_queued
        if not waiting:
            self.queue.put(None)
        for each in waiting:
            dependency = set(self.get_dependency(each))
            if dependency & failed:
                self.ever_queued.add(each)
                self.state[each]['state'] = 'failed'
                self.state[each]['used_time'] = 'FailedDependencies'
                self.logger.warning(each, 'cannot be started for some failed dependencies!')
            if not (dependency - success):
                self.ever_queued.add(each)
                self.queue.put(each, block=True)

    def _update_state(self, cmd=None):
        if cmd is not None:
            cmd_state = self.state[cmd.name]
            if cmd.proc is None:
                cmd_state['state'] = 'failed'
                cmd_state['used_time'] = 'NotEnoughResource'
                self.logger.warning(cmd.name, 'cannot be started for not enough resource!')
            else:
                cmd_state['state'] = 'success' if cmd.proc.returncode == 0 else 'failed'
                cmd_state['used_time'] = cmd.used_time
                cmd_state['mem'] = cmd.max_mem
                cmd_state['cpu'] = cmd.max_cpu
                cmd_state['pid'] = cmd.proc.pid
        success = set(x for x in self.state if self.state[x]['state'] == 'success')
        failed = set(x for x in self.state if self.state[x]['state'] == 'failed')
        running = self.ever_queued - success - failed
        waiting = set(self.names()) - self.ever_queued
        tmp_dict = {y: x for x, y in PROCESS_SET.items()}
        for each in running:
            if each in tmp_dict and psutil.pid_exists(tmp_dict[each].pid):
                if tmp_dict[each].is_running():
                    self.state[each]['state'] = 'running'
                else:
                    if tmp_dict[each].returncode == 0:
                        self.state[each]['state'] = 'success'
                    else:
                        self.state[each]['state'] = 'failed'
            else:
                self.state[each]['state'] = 'queueing'
        for each in waiting:
            self.state[each]['state'] = 'outdoor'

    def _write_state(self):
        with open(os.path.join(self.outdir, 'cmd_state.txt'), 'w') as f:
            fields = ['name', 'state', 'used_time', 'mem', 'cpu', 'pid', 'depend', 'cmd']
            f.write('\t'.join(fields)+'\n')
            for name in self.state:
                content = '\t'.join([str(self.state[name][x]) for x in fields[1:]])
                f.write(name+'\t'+content+'\n')

    def _draw_state(self):
        StateGraph(self.state).draw(os.path.join(self.outdir, 'state.svg'))

    def single_run(self):
        while True:
            if self.queue.empty():
                time.sleep(1)
                with self.__LOCK__:
                    self._update_queue()
                    self._write_state()
                    self._draw_state()
                continue
            name = self.queue.get(block=True)
            if name is None:
                self.queue.put(None)
                break
            tmp_dict = self.get_cmd_description_dict(name)

            try_times = 0
            cmd = Command(**tmp_dict, outdir=self.outdir, logger=self.logger)
            while try_times <= int(tmp_dict['retry']):
                try_times += 1
                enough = True
                if tmp_dict['check_resource_before_run']:
                    if not CheckResource().is_enough(tmp_dict['cpu'],
                                                     tmp_dict['mem'],
                                                     timeout=self.timeout):
                        enough = False
                if enough:
                    if try_times > 1:
                        self.logger.warning('{}th run {}'.format(try_times, cmd.name))
                    cmd.run()
                    if cmd.proc.returncode == 0:
                        break
            with self.__LOCK__:
                self._update_state(cmd)
                self._update_queue()
                self._write_state()
                self._draw_state()

    def parallel_run(self):
        pool_size = self.parser.getint('mode', 'threads')
        threads = list()
        for _ in range(pool_size):
            thread = threading.Thread(target=self.single_run, daemon=True)
            threads.append(thread)
            thread.start()
        # update state
        time.sleep(1)
        with self.__LOCK__:
            self._update_state()
            self._write_state()
            self._draw_state()
        # join threads
        _ = [x.join() for x in threads]
        self.logger.warning('Finished all tasks!')

    def continue_run(self, steps=''):
        detail_steps = []
        if steps:
            for each in steps:
                detail_steps += [x for x in self.names() if x == each or x.startswith(each + '_')]

        self.ever_queued = set()
        with open(os.path.join(self.outdir, 'cmd_state.txt'), 'r') as f:
            _ = f.readline()
            for line in f:
                line_lst = line.strip().split('\t')
                fields = ['state', 'used_time', 'mem', 'cpu', 'pid', 'depend', 'cmd']
                if line_lst[1] == 'success':
                    if line_lst[0] in detail_steps:
                        continue
                    self.ever_queued.add(line_lst[0])
                    self.state[line_lst[0]] = dict(zip(fields, line_lst[1:]))
        failed = set(self.names()) - self.ever_queued
        if failed:
            self.logger.warning('Continue to run: {}'.format(failed))
        else:
            self.logger.warning('Nothing to continue run')
        self.queue = queue.Queue()
        self._update_queue()
        self._draw_state()
        self.parallel_run()


if __name__ == '__main__':
    workflow = RunCommands('sample.ini', timeout=10)
    # workflow.single_run()
    # workflow.parallel_run()
    workflow.continue_run()
