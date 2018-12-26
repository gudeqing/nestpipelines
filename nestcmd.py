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
__author__ = 'gdq'


class Command(object):
    def __init__(self, cmd, name, timeout=604800,
                 monitor=True, monitor_time_step=2, **kwargs):
        self.name = name
        self.cmd = cmd
        self.proc = None
        self.stdout = None
        self.stderr = None
        self.timeout = timeout
        self.used_time = 0
        self.max_mem = 0
        self.max_cpu = 0
        self.monitor = monitor
        self.monitor_time_step = int(monitor_time_step)

    def _monitor_resource(self):
        if not isinstance(self.proc, psutil.Process):
            raise Exception('Please provide vaild process instance of psutil')
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
        self.mode = dict(self.parser['mode'])

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
            return None
        else:
            depend = self.parser[name]['depend'].strip()
            if not depend:
                return None
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
        if 'monitor' not in tmp_dict:
            tmp_dict['monitor'] = self.mode['monitor_resource']
        if 'timeout' not in tmp_dict:
            tmp_dict['timeout'] = 3600*24*7
        if 'monitor_time_step' not in tmp_dict:
            tmp_dict['monitor_time_step'] = self.mode['monitor_time_step']
        if 'check_resource_before_run' not in tmp_dict:
            tmp_dict['check_resource_before_run'] = self.mode['check_resource_before_run']
        return tmp_dict


class MonitorResource(object):
    def __init__(self, proc, time_step=1):
        self.proc = proc
        self.time_step = time_step
        self.max_mem = 0
        self.max_cpu = 0

    def monitoring(self):
        while self.proc.is_running():
            try:
                cpu_num = self.proc.cpu_num()
                if cpu_num > self.max_cpu:
                    self.max_cpu = cpu_num
                memory_obj = self.proc.memory_info()
                # memory = (memory_obj.vms - memory_obj.shared)/1024/1024
                memory = memory_obj.vms/1024/1024
                if memory > self.max_mem:
                    self.max_cpu = memory
            except Exception as e:
                print(self.proc.name(),'failed to capture resource for: ', e)
                break
            time.sleep(self.time_step)


class Resource(object):
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


class RunNestedCommands(CommandNetwork):
    def __init__(self, cmd_config):
        super(CommandNetwork).__init__(cmd_config)
        self.mode = dict(self.parser['mode'])
        self.queue = self.__init_queue()
        self.state = self.__init_state()
        self.__LOCK__ = Lock()

    def __init_queue(self):
        cmd_pool = queue.Queue()
        for each in self.orphans():
            cmd_pool.put(each)
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
        with self.__LOCK__:
            success = set(name for name in self.state if self.state[name]['state']=='success')
            failed = set(name for name in self.state if self.state[name]['state']=='failed')
            waiting = set(self.names()) - success - failed
            for each in waiting:
                dependency = set(self.get_dependency(each))
                if dependency & failed:
                    self.state[each]['state'] = 'failed'
                if not (dependency - success):
                    self.queue.put(each)

    def _update_state(self, cmd:Command):
        with self.__LOCK__:
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
        with self.__LOCK__:
            with open('cmd_state.txt', 'w') as f:
                fields = ['name', 'state', 'used_time', 'mem', 'cpu', 'pid', 'depend', 'cmd']
                f.write('\t'.join(fields)+'\n')
                for name in self.state:
                    content = '\t'.join(self.state[name][x] for x in fields[1:])
                    f.write(name+'\t'+content+'\n')

    def _send_end_signal(self):
        with self.__LOCK__:
            success = set(name for name in self.state if self.state[name]['state']=='success')
            failed = set(name for name in self.state if self.state[name]['state']=='failed')
            if len(self.names()) - len(success) - len(failed) == 0:
                self.queue.put(None)

    def _run(self):
        while 1:
            name = self.queue.get()
            if name is None:
                self.queue.put(None)
                break
            tmp_dict = self.get_cmd_description_dict(name)
            enough = True
            if tmp_dict['check_resource_before_run']:
                if not Resource().is_enough(tmp_dict['cpu'], tmp_dict['mem']):
                    enough = False
            cmd = Command(**tmp_dict)
            if enough:
                cmd.run()
            self._update_state(cmd)
            self._update_queue()
            self._send_end_signal()

    def run(self):
        pool_size = int(self.mode['threads'])
        with ThreadPoolExecutor(pool_size) as pool:
            for i in range(pool_size):
                pool.submit(self._run)









