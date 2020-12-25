分析流程概况：
    本流程仅完成call变异分析，分三个线路call变异
    1.使用mutscan搜索67个热点突变
    2.使用gencore + vardict 检测变异
    3.使用内部开发脚本consensus_reads.py 检测变异，并输出consensus fastq

一键化分析命令，需准备两个文件的输入，进入工作目录：
python = /data/users/dqgu/anaconda3/bin/python
一键化命令：usage: python batch_run.py <pipeline_template.ini> <fastq.info> [rerun 续跑时加上]
示例如下：
nohup python /data/users/dqgu/PycharmProjects/nestcmd/SmallDrugPanel/batch_run.py /data/users/dqgu/PycharmProjects/nestcmd/SmallDrugPanel/rawpipeline_for_QIA.ini fastq.info &
* 其中fastq.info为样本信息，第一列为样本名称，第二列为read1路径，第三列为read2路径
* rawpipeline_for_QIA.ini 为单样本分析流程模板

只有一个样本，且不想制作fastq.info时，可以复制和修改单样本分析流程模板，进入pipeline所在目录，然后跑如下命令：
python ~/PycharmProjects/nestcmd/basic/nestcmd.py -cfg pipeline.ini

某个样本的某些分析步骤失败时，找到原因并修改pipeline.ini后，进入pipeline所在目录，可以用如下命令续跑：
python ~/PycharmProjects/nestcmd/basic/nestcmd.py -cfg pipeline.ini --rerun

