# coding=utf-8
import time
import os
import configparser
import psutil
import queue
from subprocess import PIPE
import threading
from threading import Timer, Lock
from concurrent.futures import ThreadPoolExecutor
__author__ = 'gdq and dp'


class Command(object):
    def __init__(self, cmd, name, timeout=604800,
                 monitor_resource=True, monitor_time_step=2, **kwargs):
        self.name = name
        self.cmd = cmd
        self.proc = None
        self.stdout = None
        self.stderr = None
        self.timeout = timeout
        self.used_time = 0
        self.max_mem = 0
        self.max_cpu = 0
        self.monitor = monitor_resource
        self.monitor_time_step = int(monitor_time_step)

    def _monitor_resource(self):
        if not isinstance(self.proc, psutil.Process):
            raise Exception('Please provide valid process instance of psutil')
        mr = MonitorResource(self.proc, time_step=self.monitor_time_step)
        mr.monitoring()
        self.max_cpu = mr.max_cpu
        self.max_mem = mr.max_mem

    def run(self):
        start_time = time.time()
        print("Run {}: ".format(self.name), self.cmd)
        self.proc = psutil.Popen(self.cmd, shell=True, stderr=PIPE, stdout=PIPE)
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
        if self.stderr:
            with open(self.name+'.'+str(self.proc.pid)+'.stderr.txt', 'wb') as f:
                f.write(self.stderr)
        if self.stdout:
            with open(self.name+'.'+str(self.proc.pid)+'.stdout.txt', 'wb') as f:
                f.write(self.stdout)
        if self.max_cpu or self.max_mem:
            with open(self.name+'.'+str(self.proc.pid)+'.resource.txt', 'w') as f:
                f.write('max_cpu: {}\n'.format(self.max_cpu))
                f.write('max_mem: {}\n'.format(self.max_mem))


class CommandNetwork(object):
    def __init__(self, cmd_config):
        self.parser = configparser.ConfigParser()
        self.parser.read(cmd_config)

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
        if 'monitor' not in tmp_dict:
            tmp_dict['monitor'] = self.parser.getboolean('mode','monitor_resource')
        if 'timeout' not in tmp_dict:
            tmp_dict['timeout'] = 3600*24*7
        if 'monitor_time_step' not in tmp_dict:
            tmp_dict['monitor_time_step'] = self.parser.getint('mode', 'monitor_time_step')
        if 'check_resource_before_run' not in tmp_dict:
            tmp_dict['check_resource_before_run'] = self.parser.getboolean('mode', 'check_resource_before_run')
        return tmp_dict


class MonitorResource(object):
    def __init__(self, proc, time_step=1):
        if not isinstance(proc, psutil.Process):
            raise Exception('Please provide valid process instance of psutil')
        self.proc = proc
        self.time_step = int(time_step)
        self.max_mem = 0
        self.max_cpu = 0

    def monitoring(self):
        while self.proc.is_running():
            try:
                if os.name == 'posix':
                    cpu_num = self.proc.cpu_num()
                elif os.name == 'nt':
                    cpu_num = psutil.cpu_count()
                else:
                    cpu_num = 0
                cpu_percent = self.proc.cpu_percent(self.time_step)
                used_cpu = round(cpu_num*cpu_percent, 4)
                if used_cpu > self.max_cpu:
                    self.max_cpu = used_cpu
                memory_obj = self.proc.memory_info()
                # memory = (memory_obj.vms - memory_obj.shared)/1024/1024
                memory = memory_obj.vms/1024/1024
                if memory > self.max_mem:
                    self.max_cpu = memory
            except Exception as e:
                print('Failed to capture cpu/mem info for: ', e)
                break


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
                if enough_num >=1 and timeout <= 15:
                    return True
            if time.time() - start_time >= timeout:
                return False
            time.sleep(3)


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
                self.state[each]['state'] = 'failed'
            if not (dependency - success):
                self.ever_queued.add(each)
                self.queue.put(each, block=True)

    def _update_state(self, cmd:Command):
        cmd_state = self.state[cmd.name]
        if cmd.proc is None:
            cmd_state['state'] = 'failed'
        else:
            cmd_state['state'] = 'success' if cmd.proc.returncode==0 else 'failed'
        cmd_state['used_time'] = cmd.used_time
        cmd_state['mem'] = cmd.max_cpu
        cmd_state['cpu'] = cmd.max_cpu
        cmd_state['pid'] = cmd.proc.pid if cmd.proc else 'unknown'

    def _write_state(self):
        with open('cmd_state.txt', 'w') as f:
            fields = ['name', 'state', 'used_time', 'mem', 'cpu', 'pid', 'depend', 'cmd']
            f.write('\t'.join(fields)+'\n')
            for name in self.state:
                content = '\t'.join([str(self.state[name][x]) for x in fields[1:]])
                f.write(name+'\t'+content+'\n')

    def single_run(self):
        while True:
            name = self.queue.get(block=True)
            if name is None:
                self.queue.put(None)
                break
            tmp_dict = self.get_cmd_description_dict(name)

            try_times = 0
            cmd = Command(**tmp_dict)
            while try_times <= tmp_dict['retry']:
                try_times += 1
                enough = True
                if tmp_dict['check_resource_before_run']:
                    if not CheckResource().is_enough(tmp_dict['cpu'], tmp_dict['mem']):
                        enough = False
                if enough:
                    cmd.run()
                    if cmd.proc.returncode == 0:
                        try_times = tmp_dict['retry'] + 2
            with self.__LOCK__:
                self._update_state(cmd)
                self._update_queue()
                self._write_state()

    def parallel_run(self):
        pool_size = self.parser.getint('mode', 'threads')
        with ThreadPoolExecutor(pool_size) as pool:
            for i in range(pool_size):
                pool.submit(self.single_run)

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
        self.queue = queue.Queue()
        self._update_queue()
        self.parallel_run()


if __name__ == '__main__':
    workflow = RunCommands('cmds.ini')
    # workflow.parallel_run()
    workflow.continue_run()







