#! /data/users/dqgu/anaconda3/bin/python
import psutil
import time
from subprocess import PIPE
import os
import threading
from threading import Timer

TMP_DIR = 'nestcmd.tmp'
if not os.path.exists(TMP_DIR):
    os.mkdir(TMP_DIR)


class MyThread(threading.Thread):
    def __init__(self, func, args=(), daemon=True):
        super().__init__()
        self.func = func
        self.args = args
        self.result = None
        self.setDaemon(daemon)

    def run(self):
        self.result = self.func(*self.args)

    def get_result(self):
        return self.result


def run(cmd, marker, timeout=3600*24*10, no_monitor=False, monitor_time_step=3):
    proc = psutil.Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    out_file = os.path.join(TMP_DIR, '{}.pid.{}'.format(marker, proc.pid))
    with open(out_file, 'w') as f:
        # in future we will get pid from file name
        # f.write(str(proc.pid))
        pass
    timer = Timer(timeout, proc.kill)
    try:
        timer.start()
        max_cpu, max_mem = 0, 0
        if not no_monitor:
            thread = MyThread(func=monitor_resource, args=(proc, monitor_time_step))
            thread.start()
        stdout, stderr = proc.communicate()
        if not no_monitor:
            thread.join()
            max_cpu, max_mem = thread.get_result()
    finally:
        timer.cancel()
    print(proc.returncode, max_cpu, max_mem, stdout, stderr)


def monitor_resource(proc, time_step=3):
    if type(proc) == str or type(proc) == int:
        proc = psutil.Process(int(proc))
    max_mem, max_cpu = 0, 0
    while proc.is_running():
        try:
            cpu_percent = proc.cpu_percent(time_step)
            used_cpu = round(cpu_percent*0.01, 4)
            if used_cpu > max_cpu:
                max_cpu = used_cpu
            memory_obj = proc.memory_full_info()
            # memory = (memory_obj.vms - memory_obj.shared)/1024/1024
            # memory = round(memory_obj.vms/1024/1024, 4)
            memory = round(memory_obj.uss/1024/1024, 4)
            if memory > max_mem:
                max_mem = memory
        except Exception as e:
            # print('Failed to capture cpu/mem info for: ', e)
            break
    # print(max_cpu, max_mem)
    return max_cpu, max_mem


def get_pid(marker):
    start = time.time()
    success = False
    while not success:
        tmp = os.listdir(TMP_DIR)
        for each in tmp:
            if each.startswith(marker):
                success = True
                os.remove(os.path.join(TMP_DIR, each))
                print(each.split('.pid.')[1])
                break
        if time.time() - start > 1000:
            raise Exception('Failed to capture pid')


def is_running(pid):
    print(psutil.Process(int(pid)).is_running())


def kill(pid):
    psutil.Process(int(pid)).kill()


def pid_exists(pid):
    print(psutil.pid_exists(int(pid)))


def available_mem():
    return psutil.virtual_memory().free


def available_cpu():
    total = psutil.cpu_count()
    return int(total - total*psutil.cpu_percent()*0.01)


def resource_is_enough(cpu, mem, timeout=10):
    start_time = time.time()
    enough_num = 0
    while True:
        if float(cpu) <= available_cpu() \
                and float(mem) <= available_mem():
            enough_num += 1
            if enough_num >= 3:
                print(True)
                return True
            if enough_num >= 1 and timeout <= 10:
                print(True)
                return True
        if time.time() - start_time >= timeout:
            print(False)
            return False
        time.sleep(3)


if __name__ == '__main__':
    class Func2Command(object):
        def __init__(self, callable_dict):
            self.call(callable_dict)

        @staticmethod
        def introduce_command(func):
            import argparse
            import inspect
            import json
            import time
            if isinstance(func, type):
                description = func.__init__.__doc__
            else:
                description = func.__doc__
            if description:
                _ = [print(x.strip()) for x in description.split('\n') if x.strip()]
                parser = argparse.ArgumentParser(add_help=False)
            else:
                parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=description)
            func_args = inspect.getfullargspec(func)
            arg_names = func_args.args
            if not arg_names:
                func()
                return
            arg_defaults = func_args.defaults
            if not arg_defaults:
                arg_defaults = []
            arg_defaults = ['None']*(len(arg_names) - len(arg_defaults)) + list(arg_defaults)
            sig = inspect.signature(func)
            for arg, value in zip(arg_names, arg_defaults):
                arg_type = sig.parameters[arg].annotation
                if arg == 'self':
                    continue
                if value == 'None':
                    if arg_type in [list, tuple, set]:
                        parser.add_argument('-' + arg, nargs='+', required=True, metavar=arg)
                    elif arg_type in [int, str, float]:
                        parser.add_argument('-' + arg, type=arg_type, required=True, metavar=arg)
                    else:
                        parser.add_argument('-'+arg, required=True, metavar=arg)
                elif type(value) == bool:
                    if value:
                        parser.add_argument('--'+arg, action="store_false", help='default: True')
                    else:
                        parser.add_argument('--'+arg, action="store_true", help='default: False')
                elif value is None:
                    parser.add_argument('-' + arg, default=value, metavar='Default:' + str(value), )
                else:
                    if sig.parameters[arg].annotation in [list, tuple]:
                        parser.add_argument('-' + arg, default=value, nargs='+', type=type(value), metavar='Default:' + str(value), )
                    else:
                        parser.add_argument('-' + arg, default=value, type=type(value), metavar='Default:' + str(value), )
            if func_args.varargs is not None:
                print("warning: *varargs is not supported, and will be neglected! ")
            if func_args.varkw is not None:
                print("warning: **keywords args is not supported, and will be neglected! ")
            args = parser.parse_args().__dict__
            # try:
            #     with open("Argument_detail.json", 'w') as f:
            #         json.dump(args, f, indent=2, sort_keys=True)
            # except IOError:
            #     print('Current Directory is not writable, thus argument log is not written !')
            # start = time.time()
            func(**args)
            # print("total time: {}s".format(time.time() - start))

        def call(self, callable_dict):
            import sys
            excludes = ['introduce_command', 'Func2Command']
            _ = [callable_dict.pop(x) for x in excludes if x in callable_dict]
            if len(callable_dict) >= 2:
                if len(sys.argv) <= 1:
                    print("The tool has the following sub-commands: ")
                    _ = [print(x) for x in callable_dict]
                    return
                sub_cmd = sys.argv[1]
                sys.argv.remove(sub_cmd)

                if sub_cmd in callable_dict:
                    self.introduce_command(callable_dict[sub_cmd])
                else:
                    print('sub-command: {} is not defined'.format(sub_cmd))
            elif len(callable_dict) == 1:
                self.introduce_command(callable_dict.pop(list(callable_dict.keys())[0]))

    callable_dict = {x: y for x, y in locals().items() if callable(y)}
    _ = [callable_dict.pop(x) for x in {'Func2Command'} if x in callable_dict]
    Func2Command(callable_dict)
