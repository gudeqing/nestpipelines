# coding=utf-8
import time
import os
import configparser
import psutil
import queue
from subprocess import PIPE
import threading
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import networkx as nx
from threading import Timer, Lock
import weakref
import atexit
__author__ = 'gdq and dp'


PROCESS_SET = weakref.WeakSet()


@atexit.register
def _kill_processes_when_exit():
    print("....Ending....")
    for proc in PROCESS_SET:
        if psutil.pid_exists(proc.pid):
            print("Shutting down running tasks...")
            print('killing process {}:{}'.format(proc.pid, proc.name()))
            proc.kill()


class Command(object):
    def __init__(self, cmd, name, timeout=604800,
                 monitor_resource=True, monitor_time_step=2, **kwargs):
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
                print('Failed to capture cpu/mem info for: ', e)
                break

    def run(self):
        start_time = time.time()
        print("Run {}: ".format(self.name), self.cmd)
        self.proc = psutil.Popen(self.cmd, shell=True, stderr=PIPE, stdout=PIPE)
        PROCESS_SET.add(self.proc)
        if self.monitor:
            thread = threading.Thread(target=self._monitor_resource, daemon=True)
            thread.start()
        timer = Timer(self.timeout, self.proc.kill)
        try:
            timer.start()
            self.stdout, self.stderr = self.proc.communicate()
        finally:
            timer.cancel()
        self._write_log()
        end_time = time.time()
        self.used_time = round(end_time - start_time, 4)

    def _write_log(self):
        if not os.path.exists('logs'):
            os.mkdir('logs')
        prefix = os.path.join('logs', self.name+'.'+str(self.proc.pid))
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
        """
        drawing
        :param state: state dict from RunCommands.state
        """
        self.state = state
        self.graph = nx.DiGraph()

    def add_edges(self):
        for target in self.state:
            sources = self.state[target]['depend'].strip()
            if sources:
                sources = sources.split(',')
                edges = zip(sources, [target]*len(sources))
                self.graph.add_edges_from(edges)

    def get_color_dict(self):
        colors = list()
        all_nodes = self.graph.nodes()
        for node in all_nodes:
            state = self.state[node]['state']
            if state == 'success':
                colors.append('lightgreen')
            elif state == 'failed':
                colors.append('orange')
            elif state == 'waiting':
                colors.append('lightgray')
            else:
                colors.append('lightblue')
        return dict(zip(all_nodes, colors))

    def get_label_dict(self):
        node_label_dict = dict()
        for each in self.graph.nodes():
            used_time = self.state[each]['used_time']
            if isinstance(used_time, str):
                if used_time == 'unknown':
                    node_label_dict[each] = ''
                else:
                    node_label_dict[each] = used_time
            elif float(used_time) <= 0:
                node_label_dict[each] = ''
            else:
                node_label_dict[each] = str(used_time) + 's'
        return node_label_dict

    def draw(self):
        self.add_edges()
        pos = nx.kamada_kawai_layout(self.graph)
        # pos = nx.spring_layout(self.graph)
        node_label_dict = self.get_label_dict()
        color_dict = self.get_color_dict()
        tmp_dict = dict()
        for k, v in color_dict.items():
            tmp_dict.setdefault(v, list())
            tmp_dict[v].append(k)
        for color, group in tmp_dict.items():
            state = self.state[group[0]]['state']
            nx.draw(self.graph, pos=pos, nodelist=group, with_labels=True,
                    node_color=color, label=state, alpha=1, width=0.7, style='dashed')

        pos_attrs = {}
        for node, coords in pos.items():
            pos_attrs[node] = (coords[0], coords[1] - 0.08)

        nx.draw_networkx_labels(self.graph, pos=pos_attrs, labels=node_label_dict,
                                font_size=8, alpha=0.8)
        # nx.draw_networkx_edges(self.graph, pos=pos, style='dashed', width=0.8,)
        plt.axis('off')
        plt.legend(loc='best', fontsize='small', markerscale=0.7, frameon=False)
        plt.savefig('state_graph.png', dpi=150, bbox_inches='tight')
        plt.close()


class RunCommands(CommandNetwork):
    __LOCK__ = Lock()

    def __init__(self, cmd_config):
        super().__init__(cmd_config)
        self.ever_queued = set()
        self.queue = self.__init_queue()
        self.state = self.__init_state()

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
                print(each, 'cannot be started for some failed dependencies!')
            if not (dependency - success):
                self.ever_queued.add(each)
                self.queue.put(each, block=True)

    def _update_state(self, cmd: Command):
        cmd_state = self.state[cmd.name]
        if cmd.proc is None:
            cmd_state['state'] = 'failed'
            cmd_state['used_time'] = 'NotEnoughResource'
            print(cmd.name, 'cannot be started for not enough resource!')
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
        for each in running:
            self.state[each]['state'] = 'running'
        for each in waiting:
            self.state[each]['state'] = 'waiting'

    def _write_state(self):
        with open('cmd_state.txt', 'w') as f:
            fields = ['name', 'state', 'used_time', 'mem', 'cpu', 'pid', 'depend', 'cmd']
            f.write('\t'.join(fields)+'\n')
            for name in self.state:
                content = '\t'.join([str(self.state[name][x]) for x in fields[1:]])
                f.write(name+'\t'+content+'\n')

    def _draw_state(self):
        StateGraph(self.state).draw()

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
            cmd = Command(**tmp_dict)
            while try_times <= int(tmp_dict['retry']):
                try_times += 1
                enough = True
                if tmp_dict['check_resource_before_run']:
                    if not CheckResource().is_enough(tmp_dict['cpu'], tmp_dict['mem']):
                        enough = False
                if enough:
                    if try_times > 1:
                        print('{}th run {}'.format(try_times, cmd.name))
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
        _ = [x.join() for x in threads]

    def continue_run(self):
        self.ever_queued = set()
        with open('cmd_state.txt', 'r') as f:
            header = f.readline()
            for line in f:
                line_lst = line.strip().split('\t')
                fields = ['state', 'used_time', 'mem', 'cpu', 'pid', 'depend', 'cmd']
                if line_lst[1] == 'success':
                    self.ever_queued.add(line_lst[0])
                    self.state[line_lst[0]] = dict(zip(fields, line_lst[1:]))
        failed = set(self.names()) - self.ever_queued
        print('continue to run: ', failed)
        self.queue = queue.Queue()
        self._update_queue()
        self.parallel_run()


if __name__ == '__main__':
    workflow = RunCommands('cmds.ini')
    # workflow.single_run()
    # workflow.parallel_run()
    workflow.continue_run()
