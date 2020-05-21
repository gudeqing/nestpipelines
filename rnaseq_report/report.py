# coding=utf-8
import sys
import os
from os.path import join
from glob import glob
import shutil

script_path = os.path.abspath(__file__)
if os.path.islink(script_path):
    script_path = os.readlink(script_path)
sys.path.append(os.path.dirname(os.path.dirname(script_path)))

from rnaseq_report.batch_cmd_generator import NestedCmd
from basic.workflow_basic import basic_arg_parser

parser = basic_arg_parser()
parser.add_argument('-project_outdir', required=False)
parser.add_argument('-group_dict', required=False)
args = parser.parse_args()
args.script_path = script_path
if args.o == join(os.getcwd(), 'Result'):
    args.o = join(os.getcwd(), 'Report')

nc = NestedCmd(args)
if (args.pipeline_cfg is None) and (args.project_outdir is None):
    raise Exception('-result_dir or -pipeline_cfg is needed! Use -h for help')


def pipeline():
    result_dir = args.project_outdir
    gene_cov_files = glob(join(result_dir, 'GeneBodyCoverage', '*.geneBodyCoverage.txt'))
    inner_distance = glob(join(result_dir, 'InnerDistance', '*.inner_distance.txt'))
    read_distribution = glob(join(result_dir, 'ReadDistribution', '*.read_distribution.txt'))
    read_dup = glob(join(result_dir, 'ReadDuplication', '*.pos.DupRate.xls'))
    tpm_saturation = glob(join(result_dir, 'TPMSaturation', '*.tpm.xls'))
    chr_read_distribution = glob(join(result_dir, 'ChrReadDistribution', '*.chromosome.alignment.stat.txt'))
    exp_tables = glob(join(result_dir, 'MergeQuantGene', 'gene.isoform.TMM.EXPR.matrix'))
    crsm = glob(join(result_dir, 'CollectRnaSeqMetrics', '*.CollectRnaSeqMetrics.xls'))
    ctpm = glob(join(result_dir, 'CollectTargetedPcrMetrics', '*.CollectTargetedPcrMetrics.xls'))
    star_log = glob(join(result_dir, 'Align', '*', '*.Log.final.out'))
    diff_tables = glob(join(result_dir, 'DiffGene', '*_vs_*.*.xls'))
    go_tables = glob(join(result_dir, 'GoEnrichGene', '*_vs_*.goea.xls'))
    kegg_tables = glob(join(result_dir, 'KeggEnrichGene', '*_vs_*.*.kegg*.xls'))

    if gene_cov_files:
        nc.gene_body_coverage_cmd(gene_cov_files)
    else:
        print('geneBodyCoverage stat file not found')

    if inner_distance:
        nc.inner_distance_cmd(inner_distance)
    else:
        print('inner_distance stat file not found')

    if read_distribution:
        nc.read_distribution_cmd(read_distribution)
    else:
        print('read_distribution stat file not found')

    if read_dup:
        nc.read_duplication_cmd(read_dup)
    else:
        print('read_dup stat file not found')

    if tpm_saturation:
        nc.exp_saturation_cmd(tpm_saturation)
    else:
        print('tpm_saturation stat file not found')

    if chr_read_distribution:
        nc.chromosome_read_distribution_cmd(chr_read_distribution)
    else:
        print('chr_read_distribution stat file not found')

    if exp_tables:
        exp_table = exp_tables[0]
        nc.exp_density_cmd(exp_table)
        nc.exp_pca_cmd(exp_table)
        nc.sample_correlation_cmd(exp_table)
    else:
        print('expression matrix file not found')

    for each in [crsm, ctpm, star_log]:
        if not each:
            print('Star alignment log file not found')
    if any([crsm, ctpm, star_log]):
        nc.merge_qc_metrics_cmd(result_dir)

    if diff_tables:
        nc.volcano(diff_tables)
    else:
        print('diff table not found')

    if go_tables:
        nc.go_bubble(go_tables)
    else:
        print('go enrichment not found')

    if kegg_tables:
        nc.kegg_bubble(kegg_tables)
    else:
        print('kegg enrichment not found')

    nc.run()
    shutil.copytree(join(os.path.dirname(script_path), '1.Workflow'), join(args.o, '1.Workflow'))


if __name__ == '__main__':
    pipeline()






