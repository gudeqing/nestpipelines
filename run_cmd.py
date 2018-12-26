# coding=utf-8
import time
import os
import configparser
import psutil
import queue
from subprocess import PIPE
import threading
from threading import Timer, Semaphore
from concurrent.futures import ThreadPoolExecutor
__author__ = 'gdq'

semaphore = Semaphore(1)


class Command(object):
    def __init__(self, cmd, name, timeout=3600*24*7,
                 monitor=True, monitor_time_step=2, **kwargs):
        self.proc = None
        self.cmd = cmd
        self.name = name
        self.stderr = None
        self.stdout = None
        self.returncode = None
        self.timeout = int(timeout)
        self.max_mem = 0
        self.max_cpu = 0
        self.monitor = monitor
        self.monitor_time_step = int(monitor_time_step)
        self.used_time = 0

    def run(self):
        start_time = time.time()
        if self.monitor:
            thread = threading.Thread(target=self.monitor_resource, daemon=True)
            thread.start()
        self.proc = psutil.Popen(self.cmd, shell=True, stderr=PIPE, stdout=PIPE)
        print("Run {}: ".format(self.name), self.cmd)
        timer = Timer(self.timeout, self.proc.kill)
        try:
            timer.start()
            self.stdout, self.stderr = self.proc.communicate()
        except Exception as e:
            print(e)
        finally:
            timer.cancel()
        self.returncode = self.proc.returncode
        self.write_log()
        end_time = time.time()
        self.used_time = round(end_time - start_time, 4)
        return self.returncode, self.max_mem, self.max_cpu, self.used_time

    def monitor_resource(self):
        max_cpu = 0
        max_mem = 0
        while True:
            if self.proc is None:
                continue
            if self.proc.is_running():
                try:
                    cpu_num = self.proc.cpu_num()
                    if cpu_num > max_cpu:
                        max_cpu = cpu_num
                except Exception as e:
                    print(self.name, ': failed to capture cpu info, maybe for being too fast')

                try:
                    memory_obj = self.proc.memory_info()
                    # memory = (memory_obj.vms - memory_obj.shared)/1024/1024
                    memory = memory_obj.vms/1024/1024
                    if memory > max_mem:
                        max_cpu = memory
                except Exception as e:
                    print(self.name, ': failed to capture memory info, maybe for being too fast')

                if time.time() - self.proc.create_time() > self.timeout:
                    self.max_cpu = max_cpu
                    self.max_mem = max_mem
                    break
            else:
                self.max_cpu = max_cpu
                self.max_mem = max_mem
                break
            time.sleep(self.monitor_time_step)

    def write_log(self):
        if self.stderr:
            with open(self.name+'.'+str(self.proc.pid)+'.stderr', 'wb') as f:
                f.write(self.stderr)
        if self.stdout:
            with open(self.name+'.'+str(self.proc.pid)+'.stdout', 'wb') as f:
                f.write(self.stdout)
        if self.max_cpu or self.max_mem:
            with open(self.name+'.'+str(self.proc.pid)+'.resource', 'w') as f:
                f.write('max_cpu: {}\n'.format(self.max_cpu))
                f.write('max_mem: {}\n'.format(self.max_mem))


class Resource(object):
    @staticmethod
    def available_mem():
        return psutil.virtual_memory().free*0.9

    @staticmethod
    def available_cpu():
        total = psutil.cpu_count()
        return int(total - total*psutil.cpu_percent())-1

    def is_enough(self, cpu, mem, timeout=10):
        start_time = time.time()
        enough_num = 0
        while True:
            if float(cpu) <= self.available_cpu() \
                    and float(mem) <= self.available_mem():
                enough_num += 1
                if enough_num >= 3:
                    return True
                if enough_num >=1 and timeout <= 10:
                    return True
            if time.time() - start_time >= timeout:
                return False
            time.sleep(3)


class RunNestedCmd():
    def __init__(self, cmd_config):
        self.parser = configparser.ConfigParser()
        self.parser.read(cmd_config)
        self.all = self.parser.sections()
        self.all.pop(self.all.index('mode'))
        self.depend = dict()
        self.queue = self._init_queue()
        self.success = list()
        self.failed = list()
        self.retry = self.parser.getboolean('mode', 'retry')
        self.threads = self.parser.getint('mode', 'threads')
        self.check_resource = self.parser.getboolean('mode', 'check_resource')
        self.monitor_resource = self.parser.getboolean('mode', 'monitor_resource')
        self.checked_list = list()
        self.ever_queued = set()
        self.pid = dict()
        self.time = dict()
        self.cpu =  dict()
        self.mem = dict()

    def _init_queue(self):
        cmd_pool = queue.Queue()
        depend_dict = dict()
        for name in self.all:
            if 'depend' not in self.parser[name]:
                cmd_pool.put(name)
                depend_dict[name] = None
            else:
                depend = self.parser[name]['depend'].strip()
                if not depend:
                    cmd_pool.put(name)
                    depend_dict[name] = None
                else:
                    depend_dict[name] = [x.strip() for x in depend.split(',')]
        self.depend = depend_dict
        if cmd_pool.empty():
            raise Exception("Please make sure that at least one CMD has no dependency")
        return cmd_pool

    def _update_queue(self):
        for each in set(self.all) - set(self.success) - set(self.failed):
            if not (set(self.depend[each]) - set(self.success)):
                if set(self.depend[each]) & set(self.failed):
                    self.failed.append(each)
                    continue
                if each not in self.ever_queued:
                    print(each, ": It's dependencies finished!")
                    self.ever_queued.add(each)
                    self.queue.put(each)

    def _log_state(self):
        run_or_wait = set(self.all) - set(self.success) - set(self.failed)
        running = set()
        waiting = set()
        for each in run_or_wait:
            if each in self.pid:
                if psutil.pid_exists(self.pid[each]):
                    running.add(each)
            else:
                waiting.add(each)
        with open('state_detail.txt', 'w') as f:
            f.write('name\tstate\tused_time\tmemory\tcpu\tpid\tcmd\n')
            for each in self.all:
                if each in self.success:
                    state = 'success'
                elif each in self.failed:
                    state = 'failed'
                elif each in running:
                    state = 'running'
                elif each in waiting:
                    state = 'waiting'
                else:
                    state = 'unknown'
                used_time = self.time[each] if each in self.time else 'unknown'
                mem = self.mem[each] if each in self.mem else 'unknown'
                cpu = self.cpu[each] if each in self.cpu else 'unknown'
                pid = self.pid[each] if each in self.pid else 'unknown'
                cmd = self.parser.get(each, 'cmd')
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(each, state, used_time, mem, cpu, pid, cmd))

    def run_cmd(self, name):
        tmp_dict = dict(self.parser[name])
        tmp_dict['name'] = name
        if 'monitor' not in tmp_dict:
            tmp_dict['monitor'] = self.monitor_resource
        if 'timeout' not in tmp_dict:
            tmp_dict['timeout'] = 3600*24*7
        if 'monitor_time_step' not in tmp_dict:
            tmp_dict['monitor_time_step'] = 2
        cmd = Command(**tmp_dict)
        cmd.run()
        with semaphore:
            if cmd.returncode == 0:
                self.success.append(name)
            else:
                if self.retry:
                    cmd.run()
                if cmd.returncode == 0:
                    self.success.append(name)
                else:
                    self.failed.append(name)
            self.pid[name] = cmd.proc.pid
            self.time[name] = cmd.used_time
            self.cpu[name] = cmd.max_cpu
            self.mem[name] = cmd.max_mem
            self._update_queue()
            self._log_state()

    def run_all(self):
        def run(check_resource=self.check_resource):
            while True:
                name = self.queue.get(block=True)
                if name is None:
                    self.queue.put(None)
                    break
                if check_resource:
                    cpu = self.parser.get(name, 'cpu')
                    mem = self.parser.get(name, 'mem')
                    if not Resource().is_enough(cpu, mem):
                        with semaphore:
                            if self.checked_list.count(name) >= 5:
                                self.failed.append(name)
                                raise Exception(name+" cannot be started after checking resource!")
                            else:
                                self.checked_list.append(name)
                                print(name, ': get back to wait for resource')
                                self.queue.put(name, block=True)
                                continue
                self.run_cmd(name)
                time.sleep(1)
                with semaphore:
                    if len(self.success) + len(self.failed) == len(self.all):
                        self.queue.put(None)
                        break

        with ThreadPoolExecutor(self.threads) as pool:
            for i in range(self.threads):
                pool.submit(run)

    def continue_run(self):
        if not os.path.exists('state_detail.txt'):
            exit('Cannot find previous log file: state_detail.txt')
        with open('state_detail.txt', 'r') as f:
            success = list()
            _ = f.readline()
            for line in f:
                name, state, _ = line.split('\t', 3)
                if state == 'success':
                    success.append(name)
        self.success = success
        self.queue = queue.Queue()
        self._update_queue()
        if self.queue.empty():
            exit('No need to continue for all steps are completed!')
        self.run_all()


if __name__ == '__main__':
    a = RunNestedCmd('cmds.ini')
    a.run_all()
    # a.continue_run()
