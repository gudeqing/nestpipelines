RNASeqPipeline使用说明

1. quick start
python: /data/users/dqgu/anaconda3/bin/python
先查看流程帮助信息:
python /data/users/dqgu/PycharmProjects/nestcmd/RnaSeqPipeline -h
usage: RnaSeqPipeline.py [-h] [-arg_cfg ARG_CFG] [-o O]
                         [-skip SKIP [SKIP ...]] [--only_show_steps]
                         [--only_show_detail_steps] [--only_write_pipeline]
                         [-threads THREADS] [-retry RETRY] [--continue_run]
                         [-rerun_steps RERUN_STEPS [RERUN_STEPS ...]]
                         [-pipeline_cfg PIPELINE_CFG] [--list_cmd_names]
                         [-show_cmd_example SHOW_CMD_EXAMPLE]
                         [--no_monitor_resource]
                         [--monitor_time_step MONITOR_TIME_STEP]
                         [-wait_resource_time WAIT_RESOURCE_TIME]
                         [--no_check_resource_before_run] [-hostname HOSTNAME]
                         [-port PORT] [-username USERNAME]
                         [-password PASSWORD] [-fastq_info FASTQ_INFO]
                         [-group GROUP] [-compare COMPARE]
optional arguments:
  -h, --help            show this help message and exit
  -arg_cfg ARG_CFG      参数配置文件, 包含流程所有软件的需要的参数, 默认使用该脚本当前目录的argument.ini
  -o O                  分析目录或结果目录
  -skip SKIP [SKIP ...]
                        指定要跳过的步骤名, 空格分隔,程序会自动跳过依赖他们的步骤,
                        --only_show_steps可查看步骤名; 注意: 对于rnaseq_pipeline,
                        如果跳过trim步骤, 则用原始数据进行后续分析; 如果跳过assembl或mergeTranscript,
                        则仅对参考基因定量
  --only_show_steps     仅仅显示当前流程包含的主步骤, 且已经排除指定跳过的步骤; 你还可以通过--
                        list_cmd_names查看当前流程包含哪些命令行
  --only_show_detail_steps
                        仅仅显示当前流程包含的详细步骤, 且已经排除指定跳过的步骤
  --only_write_pipeline
                        仅仅生成流程pipeline.ini
  -threads THREADS      允许的最大并行的cmd数目, 默认5
  -retry RETRY          某步骤运行失败后再尝试运行的次数, 默认1次. 如需对某一步设置不同的值,
                        可在运行流程前修改pipeline.ini
  --continue_run        流程运行结束后, 从失败的步骤续跑, 记得要用-o指定之前的结果目录,
                        用-pipeline_cfg指定pipeline.ini; 如果顺便还想重新跑已经成功运行的步骤, 可通过-
                        rerun_steps指定, 或者在状态表cmd_stat.txt中将其修改为failed即可
  -rerun_steps RERUN_STEPS [RERUN_STEPS ...]
                        使用--continue_run有效, 通过该参数指定重跑已经成功的步骤, 空格分隔,
                        这样做的可能原因可以是: 你重新设置了参数
  -pipeline_cfg PIPELINE_CFG
                        已有的pipeline.ini, 续跑时必须需提供; 如提供该参数, 则此时无视
                        arg_cfg,fastq_info,group,cmp,skip 等用于生成项目流程的参数
  --list_cmd_names      仅输出参数配置文件里的包含的cmd名称
  -show_cmd_example SHOW_CMD_EXAMPLE
                        提供一个cmd名称,输出该cmd的样例
  --no_monitor_resource
                        是否监控每一步运行时的资源消耗, 如需对某一步设置不同的值, 可在运行流程前修改pipeline.ini
  --monitor_time_step MONITOR_TIME_STEP
                        监控资源时的时间间隔, 默认3秒, 如需对某一步设置不同的值, 可在运行流程前修改pipeline.ini
  -wait_resource_time WAIT_RESOURCE_TIME
                        等待资源的时间上限, 默认1500秒, 等待时间超过这个时间时,资源不足时判定任务失败
  --no_check_resource_before_run
                        指示运行某步骤前检测指定的资源是否足够, 如不足, 则该步骤失败; 如果设置该参数, 则运行前不检查资源.
                        如需对某一步设置不同的值,可运行前修改pipeline.ini. 如需更改指定的资源,
                        可在运行流程前修改pipeline.ini
  -hostname HOSTNAME    其他服务器地址,默认 10.60.2.133
  -port PORT            服务器端口号, 默认 22
  -username USERNAME    登录其他服务器的用户名
  -password PASSWORD    登录密码, 使用前请确定两个服务器具备共享文件的能力
  -fastq_info FASTQ_INFO
                        第一列为样本名,分析结果使用名; 第二列为read1的fastq路径,如有多个,需要分号分隔;
                        第三列为可选, read2的fastq路径,单端测序时则无第三列
  -group GROUP          样本分组信息文件,至少两列,第一列样本名,第二列为分组名,其他列也是分组名
  -compare COMPARE      比较信息文件,两列,第1列是对照组名,第2列是实验组名
  注意：带有'--'前缀的参数表示该参数是bool类型，你无需给他提供值，即如果提供该参数，则该参数值自动转为True

2. 运行流程cmd示例
python /data/users/dqgu/PycharmProjects/nestcmd/rnaseq_pipeline/RnaSeqPipeline.py -fastq fastq.info -group group.info -compare cmp.info -skip NGSCheckMate FastqToSam ArribaAlign MapSplice MergeBamAlignment Assembly
示例说明：
    （1）fastq.info是raw fastq的路径信息，通常包含3列，第一列是样本名，也是后续分析所用名， 第二列和第三列分别是fq1和fq2的绝对路径
    该文件可以由python /data/users/dqgu/PycharmProjects/biodev/smallScriptsAtDunwill/findFastqAndLink.py get_all_fastq_abs_path 命令搜索某个文件夹并匹配指定的fastq文件而生成。
    （2）group.info是样本分组信息，第一列为样本名，第二列为分组信息，主要用于差异分析，需要header
    （3）cmp.info 是比较信息，第一列为对照组名，第二列为实验组名，用于差异分析，无需header
    （4）skip参数，可以赋予多个步骤名，一般分析跳过NGSCheckMate FastqToSam ArribaAlign MapSplice MergeBamAlignment Assembly这些步骤
    （5）如果你想调整参数，复制默认参数配置文件并修改后传给流程即可：/data/users/dqgu/PycharmProjects/nestcmd/rnaseq_pipeline/arguments.ini

3. 查看流程运行状态
进入流程工作目录，找到cmd_state.txt文件，可以看到详细的cmd信息和运行状态，另外也可以通过state.svg文件查看，该文件需要通过浏览器打开，最后还可以通过日志文件查看

4. 流程续跑
假设流程工作/结果目录为Result/
续跑失败的步骤，示例如下：
python /data/users/dqgu/PycharmProjects/nestcmd/rnaseq_pipeline/RnaSeqPipeline.py -pipeline Result/pipeline.ini --continue_run
续跑失败的步骤，并且重跑部分成功的步骤（这样做的原因可能是想进行参数调整，你可以在Result/pipeline.ini直接进行参数调整，如果想批量调整参数，可以通过--only_write_pipeline参数重新生成流程pipeline.ini）,示例如下：
python /data/users/dqgu/PycharmProjects/nestcmd/rnaseq_pipeline/RnaSeqPipeline.py -pipeline Result/pipeline.ini --continue_run -rerun DiffGene KeggEnrichGene

5. 结果整理及报告
假设结果目录为Result/（使用绝对路径）, 则运行如下命令，默认生成Report目录
python /data/users/dqgu/PycharmProjects/nestcmd/rnaseq_report/report.py -project_outdir Result/
进入Report目录，修改或增加必要信息，如：
    确认1.Workflow的信息是否正确，
    在2.QCSummary增加或修改样本信息，删除不必要的质控结果统计表，
    在3.QCFigures中删除不必要的质控图，由于网页版的图占很大内存，影响报告发生，可以考虑替换图片，尤其是饱和度分析的图等，
    可以考虑删除inner distribution图
    删除1.ExpPCA和2.ExpDensity里不必要的数据表，保留图片即可
    建议保留：1.GeneBodyCoverage/         2.ReadDistribution/         3.ChrReadDistribution/      4.ExpSaturation/
    手动增加4.ExpAnalysis/4.ExpressionMatrix，可参考项目80011077_RNA/Report/4.ExpAnalysis
    删除报告流程生成的日志文件等信息
    其他

6.生成网页报告
整理和确认好Report目录后，回到report目录，执行如下命令生成网页版报告
python /data/users/dqgu/PycharmProjects/nestcmd/rnaseq_report/plotly/makeHtmlReport.py make_report -cfg_from ./
执行成功后，会生成index.html文件，把整个报告目录复制到window，双击index.html即可浏览网页报告，最后的结果目录可能如下
1.Workflow  2.QCSummary  3.QCFigures  4.ExpAnalysis  5.DiffGene  6.GoEnrichment  7.KEGGEnrichment  html.utils  index.html
可手动删除report.yml文件

Questions
1. 如何知道该流程中包含哪些步骤？
使用--only_show_steps参数查看
2. 如何获取和修改流程参数信息？
复制/data/users/dqgu/PycharmProjects/nestcmd/rnaseq_pipeline/arguments.ini到本地即可查看和修改
3. 如何仅仅生成流程而不运行流程？
输入fastq等必须参数后，加上参数--only_write_pipeline即可,在结果目录会有一个pipeline.ini的文件，这个文件就是所谓的流程。
如果需要运行该流程，可用--continue_run 进行续跑，用-pipeline_cfg参数指定该流程
4. 如何获得每个运行步骤对应的命令行
查看结果文件中的cmd_state参数，或者查看流程pipeline.ini

