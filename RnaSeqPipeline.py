# coding=utf-8
import os
from rna_tools import NestedCmd
from workflow_basic import basic_arg_parser

parser = basic_arg_parser()
parser.add_argument('-fastq_info', required=False,
                    help="第一列为样本名,分析结果使用名; 第二列为read1的fastq路径,如有多个,需要分号分隔; "
                         "第三列为可选, read2的fastq路径,单端测序时则无第三列")
parser.add_argument('-group', help="样本分组信息文件,至少两列,第一列样本名,第二列为分组名,其他列也是分组名")
parser.add_argument('-compare', help="比较信息文件,两列,第1列是对照组名,第2列是实验组名")
args = parser.parse_args()

if args.only_show_steps:
    script_path = os.path.abspath(__file__)
    if args.fastq_info is None:
        test_data_dir = os.path.join(os.path.dirname(script_path), 'testdata')
        if os.path.exists(test_data_dir):
            args.fastq_info = os.path.join(test_data_dir, 'fastq_info.txt')
            args.compare = os.path.join(test_data_dir, 'compare')
            args.group = os.path.join(test_data_dir, 'group')

nc = NestedCmd(args)
terminate = nc.do_some_pre_judge(nc.workflow_arguments)
if terminate:
    exit(0)
if (args.pipeline_cfg is None) and (not args.fastq_info):
    raise Exception('-fastq_info or -pipeline_cfg is needed! Use -h for help')


def pipeline():
    """
    **为了能正常跳过一些步骤,步骤名即step_name不能包含'_',且保证最后生成的步骤名不能有重复**
    """
    # fastqc and trimmomatic
    fastq_info_dict = nc.parse_fastq_info(nc.workflow_arguments.fastq_info)
    nc.fastqc_raw_data_cmds(fastq_info_dict, step_name='RawDataQC')
    trim_cmds = nc.trimmomatic_cmds(fastq_info_dict, step_name='Trim')
    fastqc_cmds = nc.fastqc_trimmed_data_cmds(trimming_cmds=trim_cmds, step_name='TrimmedDataQC')

    # align and bam index
    star_indexing = nc.star_index_cmd(step_name='AlignIndex')
    
    if list(trim_cmds.keys())[0].split('_', 1)[0] in nc.workflow_arguments.skip:
        align_cmds = nc.star_align_with_rawdata_cmds(fastq_info_dict, index_cmd=star_indexing, step_name='Align')
    else:
        align_cmds = nc.star_align_cmds(trimming_cmds=trim_cmds, index_cmd=star_indexing, step_name='Align')
    
    bam_indexing_cmds = nc.bam_index_cmds(align_cmds, step_name='IndexBam')

    # run some RseQC cmds
    gbc_cmds = nc.gene_body_coverage_cmds(bam_indexing_cmds, step_name='GeneBodyCoverage')
    inner_dist_cmds = nc.inner_distance_cmds(bam_indexing_cmds, step_name='InnerDistance')
    rd_cmds = nc.read_distribution_cmds(bam_indexing_cmds, step_name='ReadDistribution')
    rdup_cmds = nc.read_duplication_cmds(bam_indexing_cmds, step_name='ReadDuplication')
    frag_size_cmds = nc.rna_fragment_size_cmds(bam_indexing_cmds, step_name='FragmentSize')

    # run TPM saturation
    saturation_cmds = nc.tpm_saturation_cmds(bam_indexing_cmds, step_name='TPMSaturation')

    # stat chromosome read distribution
    chr_read_distribution_cmds = nc.chromosome_read_distribution_cmds(bam_indexing_cmds)

    # 根据star比对结果log文件和使用bedtools intersect 结合 samtools flagstat 统计比对结果
    # alignment_summary_cmds = nc.get_alignment_summary_cmds(bam_indexing_cmds, step_name='AlignmentSummary')

    # run some picard tools 获得大量的比对情况的统计结果，包括insert_size 和 gene_body_coverage
    collect_alignment_summary_cmds = nc.CollectAlignmentSummaryMetrics_cmds(bam_indexing_cmds)
    collect_insert_size_cmds = nc.CollectInsertSizeMetrics_cmds(bam_indexing_cmds)
    collect_target_info_cmds = nc.CollectTargetedPcrMetrics_cmds(bam_indexing_cmds)
    collect_rna_metric_cmds = nc.CollectRnaSeqMetrics_cmds(bam_indexing_cmds)

    # assemble and merge
    assembly_cmds = nc.scallop_cmds(align_cmds=align_cmds, step_name='Assembly')
    merge_trans_cmd = nc.merge_scallop_transcripts_cmd(assembly_cmds, step_name='MergeTranscript')

    # quant with salmon and merge result
    if list(merge_trans_cmd.keys())[0] in nc.workflow_arguments.skip or list(assembly_cmds.keys())[0].split('_', 1)[0] in nc.workflow_arguments.skip:
        salmon_indexing = nc.salmon_index_cmd(merge_transcript_cmd=None, step_name='QuantIndex')
    else:
        salmon_indexing = nc.salmon_index_cmd(merge_transcript_cmd=merge_trans_cmd, step_name='QuantIndex')
    
    if list(trim_cmds.keys())[0].split('_', 1)[0] not in nc.workflow_arguments.skip:
        quant_cmds = nc.salmon_quant_with_clean_data_cmds(trim_cmds, index_cmd=salmon_indexing, step_name='Quant')
    else:
        quant_cmds = nc.salmon_quant_with_raw_data_cmds(fastq_info_dict, index_cmd=salmon_indexing, step_name='Quant')
    
    merge_gene_exp_cmd = nc.merge_quant_cmd(quant_cmds=quant_cmds, step_name="MergeQuant")

    # expression analysis
    gene_exp_cluster_cmd = nc.exp_analysis_cmd(merge_gene_exp_cmd, step_name='ExpAnalysis', level='gene')
    merge_trans_exp_cmd = nc.merge_quant_cmd(quant_cmds=quant_cmds, level='transcript', step_name='MergeQuant')
    trans_exp_cluster_cmd = nc.exp_analysis_cmd(merge_trans_exp_cmd, step_name='ExpAnalysis', level='transcript')

    # diff and enrich
    if nc.workflow_arguments.group and nc.workflow_arguments.compare:
        diff_gene_cmd = nc.diff_exp_cmd(merge_gene_exp_cmd, level='gene', step_name='Diff')
        gene_go_enrich_cmds = nc.go_enrich_cmds(diff_gene_cmd, level='gene', step_name='GoEnrich')
        gene_kegg_enrich_cmds = nc.kegg_enrich_cmds(diff_gene_cmd, level='gene', step_name='KeggEnrich')
        
        diff_trans_cmd = nc.diff_exp_cmd(merge_trans_exp_cmd, level='transcript', step_name='Diff')
        trans_go_enrich_cmds = nc.go_enrich_cmds(diff_trans_cmd, level='transcript', step_name='GoEnrich')
        trans_kegg_enrich_cmds = nc.kegg_enrich_cmds(diff_trans_cmd, level='transcript', step_name='KeggEnrich')

    # skip, show some help or run
    nc.run()


if __name__ == '__main__':
    pipeline()
