0. 本流程主要是对vdjtools的包装
1. 从https://10.62.2.16/ir/secure/analyses.html下载数据，可以选择多个样本后批量下载
2. unzip下载的数据，放在目录FromIR中，每个样本占一个子目录，里面包含了IR对每个样本的分析结果
3. 合并不同样本的分析数据：
python ~/PycharmProjects/nestcmd/tcr_pipeline/utils/merge_tcr.py merge_metric_matrix -file ../FromIR/*/RESULTS/*.metrics.csv -group sample.info.txt -new_name_col sample
该部分生成的结果可直接用于报告
4. 把IR的数据转换为vdjtools的输入, 也就是生成一个metadata.txt文件作为流程的输入
python ~/PycharmProjects/nestcmd/tcr_pipeline/utils/merge_tcr.py convert2vdjtools -file FromIR/*/RESULTS/*.clone_summary.csv -g sample.info.txt -o vdjtools_input
5. 跑流程：
cp ~/PycharmProjects/nestcmd/tcr_pipeline/arguments.ini .
修改arguments.ini
nohup python ~/PycharmProjects/nestcmd/tcr_pipeline/TcrPipeline.py -arg arguments.ini &


