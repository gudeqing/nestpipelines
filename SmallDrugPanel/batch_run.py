import os
import sys
import configparser
from concurrent.futures import ThreadPoolExecutor
import subprocess
import psutil
import weakref
from threading import Timer
import atexit
PROCESS_local = weakref.WeakKeyDictionary()


def generate_cmds():
    if len(sys.argv) < 3:
        print('usage: python batch_run.py <pipeline_template.ini> <fastq.info> [rerun 续跑时加上]')
        exit()
    if len(sys.argv) >= 4:
        rerun = True
    else:
        rerun = False

    pipeline_ini = sys.argv[1]
    fastq_info = sys.argv[2]
    cfp = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    cfp.read(pipeline_ini, encoding='utf-8')
    cmds = []
    with open(fastq_info) as f:
        for line in f:
            lst = line.strip().split()
            sample = lst[0]
            cfp['mode']['sample'] = sample
            cfp['mode']['fq'] = lst[1]
            cfp['mode']['fq2'] = lst[2]
            os.makedirs(sample, exist_ok=True)
            abspath = os.path.abspath(sample)
            pipeline_path = f'{sample}/{lst[0]}.pipeline.ini'
            cfp.write(open(pipeline_path, 'w'))
            pipeline_path = os.path.abspath(pipeline_path)
            if rerun:
                cmd = f'python /data/users/dqgu/PycharmProjects/nestcmd/basic/nestcmd.py -cfg {pipeline_path} --rerun'
            else:
                cmd = f'python /data/users/dqgu/PycharmProjects/nestcmd/basic/nestcmd.py -cfg {pipeline_path}'
            cmds.append([sample, abspath, cmd])
    return cmds


def run_cmd(cmd):
    os.chdir(cmd[1])
    proc = psutil.Popen(cmd[2], shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    PROCESS_local[proc] = cmd[0]
    timer = Timer(3600*36, proc.kill)
    try:
        timer.start()
        print(f'Start running {cmd[0]}')
        stdout, stderr = proc.communicate()
    finally:
        timer.cancel()


@atexit.register
def _kill_processes_when_exit():
    print("....Ending....")
    for proc, cmd_name in PROCESS_local.items():
        if psutil.pid_exists(proc.pid):
            print('Shutting down running tasks {}:{}'.format(proc.pid, cmd_name))
            proc.kill()


if __name__ == '__main__':
    cmds = generate_cmds()
    with ThreadPoolExecutor(3) as pool:
        pool.map(run_cmd, cmds)
