0. 本流程主要是对vdjtools的包装

1. 从https://10.62.2.16/ir/secure/analyses.html下载数据，可以选择多个样本后批量下载

2. unzip下载的数据，放在目录FromIR中，每个样本占一个子目录，里面包含了IR对每个样本的分析结果：先解压出每个样本的压缩文件，再对循环对每个压缩文件进行解压
cd FromIR
for i in `ls |sed 's/_All.zip//'`; do unzip ${i}_All.zip -d $i; done
rm *_All.zip

3. 合并不同样本的分析数据：
cd ../
mkdir MergeMetrics
cd MergeMetrics
python ~/PycharmProjects/nestcmd/tcr_pipeline/utils/merge_tcr.py merge_metric_matrix -file ../FromIR/*/RESULTS/*.metrics.csv -group sample.info.txt -new_name_col sample -factor_col group
* 该脚本将自动使用'_'作为分隔符从目录名中提取出第一个字符串作为样本id
* 输入文件中的sample.info.txt的第一列必须是样本id，如果不确定id长啥样，可以先不使用该参数得到结果后查看。
* -new_name_col: 指定group_info中的一列, 作为样本别名，将放在第二列，这列信息将作为分析用名
* -factor_col: 指定group_info中的一列, 作为factor列，也就是分组信息, 该列信息要带入到Report.diversity.summary.csv
* 输出4个文件，有metric,count, frequency， metric.report
* 该步骤生成的目录 1.SampleInfo和2.SampleQC和1.DiversitySummary用于后续报告整理

4. 把IR的数据转换为vdjtools的输入, 也就是生成一个metadata.txt文件作为流程的输入
python ~/PycharmProjects/nestcmd/tcr_pipeline/utils/merge_tcr.py convert2vdjtools -file FromIR/*/RESULTS/*.clone_summary.csv -g sample.info.txt -o vdjtools_input
* 输入的sample.info.txt和第三步输入的一样
* 该步骤的输出的目录可直接用于后续报告整理

5. 跑流程：
cp ~/PycharmProjects/nestcmd/tcr_pipeline/arguments.ini .
修改arguments.ini
nohup python ~/PycharmProjects/nestcmd/tcr_pipeline/TcrPipeline.py -arg arguments.ini &

6. 整理报告目录, 主要分析结果还是以IR的分析结果为主，因此很多vdjtools的分析结果并不用于汇报。
进入分析结果目录
mkdir ../Report
cp -r 1.QC/ 2.Diversity/ 3.VJ-Usage/ 4.CDR3Stat/ 5.Cluster/ ../Report/
cd ../Report

cd 1.QC
删除已有目录，移动步骤3的 1.SampleInfo  2.SampleQC 到该目录

cd ../2.Diversity
rm -r 1.DiversitySummary/
移动步骤3的1.DiversitySummary到该目录
移动步骤4的输出目录3.CloneDetail到该目录

cd ../5.Cluster
rm -r 2.PairwiseDistances/
结束

7. 生成报告
cd Report
cp /data/users/dqgu/PycharmProjects/nestcmd/tcr_report/content.ini .
修改content.ini，对结果目录请使用相对路径，以免在报告中带入的路径是linux下的绝对路径
python /data/users/dqgu/PycharmProjects/nestcmd/tcr_report/generate.py content.ini
得到word版报告后，查看world内容，尤其是表格，可以手动调整布局，最后使用world的“引用|目录”功能在第一页末尾创建文件目录
world转换成pdf版用于报告

8. 特殊报告内容
手动整理，用ppt或其他形式展示