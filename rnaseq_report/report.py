# coding=utf-8
import sys
import os
from os.path import join
from glob import glob

script_path = os.path.abspath(__file__)
if os.path.islink(script_path):
    script_path = os.readlink(script_path)
sys.path.append(os.path.dirname(os.path.dirname(script_path)))

from rnaseq_report.batch_cmd_generator import NestedCmd
from basic.workflow_basic import basic_arg_parser

parser = basic_arg_parser()
parser.add_argument('-result_dir', required=False)
parser.add_argument('-group_dict', required=False)
args = parser.parse_args()
args.script_path = script_path
if args.o == join(os.getcwd(), 'Result'):
    args.o = join(os.getcwd(), 'Report')

nc = NestedCmd(args)
if (args.pipeline_cfg is None) and (args.result_dir is None):
    raise Exception('-result_dir or -pipeline_cfg is needed! Use -h for help')


def pipeline():
    result_dir = args.result_dir
    gene_cov_files = glob(join(result_dir, 'CollectRnaSeqMetrics', '*.CollectRnaSeqMetrics.xls'))
    fragment_length = glob(join(result_dir, 'FragmentSize', '*.fragment_size.txt'))
    inner_distance = glob(join(result_dir, 'InnerDistance', '*.inner_distance.txt'))
    read_distribution = glob(join(result_dir, 'ReadDistribution', '*.read_distribution.txt'))
    read_dup = glob(join(result_dir, 'ReadDuplication', '*.pos.DupRate.xls'))
    tpm_saturation = glob(join(result_dir, 'TPMSaturation', '*.tpm.xls'))
    chr_read_distribution = glob(join(result_dir, 'ChrReadDistribution', '*.chromosome.alignment.stat.txt'))
    exp_tables = glob(join(result_dir, 'MergeQuantGene', 'gene.isoform.TMM.EXPR.matrix'))
    casm = glob(join(result_dir, 'CollectAlignmentSummaryMetrics', '*.CollectAlignmentSummaryMetrics.xls'))
    cism = glob(join(result_dir, 'CollectInsertSizeMetrics', '*.CollectInsertSizeMetrics.xls'))
    crsm = glob(join(result_dir, 'CollectRnaSeqMetrics', '*.CollectRnaSeqMetrics.xls'))
    ctpm = glob(join(result_dir, 'CollectTargetedPcrMetrics', '*.CollectTargetedPcrMetrics.xls'))

    if gene_cov_files:
        gene_cov_cmd = nc.gene_body_coverage_cmd(gene_cov_files)
        gene_cov_slider = nc.make_gene_body_cov_slider(gene_cov_cmd)
    else:
        print('geneBodyCoverage stat file not found')

    if fragment_length:
        fragment_length_cmd = nc.fragment_length_cmd(fragment_length)
        fragment_length_slider = nc.make_fragment_size_slider(fragment_length_cmd)
    else:
        print('fragment_length stat file not found')

    if inner_distance:
        inner_distance_cmd = nc.inner_distance_cmd(inner_distance)
        inner_distance_slider = nc.make_inner_size_slider(inner_distance_cmd)
    else:
        print('inner_distance stat file not found')

    if read_distribution:
        read_distribution_cmd = nc.read_distribution_cmd(read_distribution)
        read_distribution_slider = nc.make_read_distribution_slider(read_distribution_cmd)
    else:
        print('read_distribution stat file not found')

    if read_dup:
        read_dup_cmd = nc.read_duplication_cmd(read_dup)
        read_dup_slider = nc.make_read_duplication_slider(read_dup_cmd)
    else:
        print('read_dup stat file not found')

    if tpm_saturation:
        tpm_saturation_cmd  = nc.exp_saturation_cmd(tpm_saturation)
        tpm_saturation_slider = nc.make_exp_saturation_slider(tpm_saturation_cmd)
    else:
        print('tpm_saturation stat file not found')

    if chr_read_distribution:
        chr_read_distribution_cmd  = nc.chromosome_read_distribution_cmd(chr_read_distribution)
        chr_read_distribution_slider = nc.make_chr_read_distribution_slider(chr_read_distribution_cmd)
    else:
        print('chr_read_distribution stat file not found')

    if exp_tables:
        exp_table = exp_tables[0]
        exp_density_cmd = nc.exp_density_cmd(exp_table)
        exp_density_slider = nc.make_exp_distribution_slider(exp_density_cmd)
        exp_pca_cmd = nc.exp_pca_cmd(exp_table)
        exp_pca_slider = nc.make_pca_slider(exp_pca_cmd)
        sample_corr_cmd = nc.sample_correlation_cmd(exp_table)
        sample_corr_slider = nc.make_corr_cluster_slider(sample_corr_cmd)
    else:
        print('expression matrix file not found')

    if casm:
        casm_cmd = nc.CollectAlignmentSummaryMetrics(casm)
        casm_slider = nc.make_AlignmentSummary_slider(casm_cmd)
    else:
        print('CollectAlignmentSummaryMetrics result not found')

    if cism:
        cism_cmd = nc.CollectInsertSizeMetrics(cism)
        cism_slider = nc.make_InsertSizeSlider_slider(cism_cmd)
    else:
        print('CollectInsertSizeMetrics result not found')

    if crsm:
        crsm_cmd = nc.CollectRnaSeqMetrics(crsm)
        crsm_slider = nc.make_RnaSeqMetrics_slider(crsm_cmd)
    else:
        print('CollectRnaSeqMetrics result not found')

    if ctpm:
        ctpm_cmd = nc.CollectTargetedPcrMetrics(ctpm)
        ctpm_slider = nc.make_TargetedSummary_slider(ctpm_cmd)
    else:
        print('CollectTargetedPcrMetrics result not found')

    nc.run()


if __name__ == '__main__':
    pipeline()






