import os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import venn
import pandas as pd
from upsetplot import from_contents
from upsetplot import plot


def merge_clone_seq(metadata_file, label_field, factor_field, seq_type='CDR3nt', outdir=None):
    meta_df = pd.read_csv(metadata_file, sep='\t')
    seq_set_dict = dict()
    for each_file, sample in zip(meta_df['files'], meta_df[label_field]):
        clone_info = pd.read_csv(each_file, sep='\t')
        seq_set_dict[sample] = set(clone_info[seq_type])
    if outdir is None:
        outdir = os.getcwd()
    group = os.path.join(outdir, f'{factor_field}.group.txt')
    meta_df[[label_field, factor_field]].to_csv(group, index=False, header=False, sep='\t')
    return seq_set_dict, group


def venn_plot(venn_set_dict, set_group, out_prefix='result', graph_format='png'):
    # plot venn
    with open(set_group) as f:
        group_dict = dict(x.strip().split()[:2] for x in f)
        tmp_dict = dict()
        for k, v in group_dict.items():
            tmp_dict.setdefault(v, set())
            tmp_dict[v].add(k)
    venn_list = []
    venn_names = []
    for k, v in tmp_dict.items():
        venn_list.append(','.join(v))
        venn_names.append(k)

    for group, name in zip(venn_list, venn_names):
        groups = group.split(',')
        tmp_dict = {x: y for x, y in venn_set_dict.items() if x in groups}
        if 2 <= len(tmp_dict) <= 6:
            venn.venn(tmp_dict, cmap="tab10", fmt="{size}\n{percentage:.2f}%", fontsize=9)
            out_name = out_prefix + '.{}.venn.{}'.format(name, graph_format)
            plt.savefig(out_name, dpi=300)
            plt.close()
        else:
            print('venn for {}?'.format(groups))
            print('venn only support 2-6 sets')

    # intersection plot
    if venn_list is None:
        if len(venn_set_dict) <= 8:
            plot(from_contents(venn_set_dict), sum_over=False, sort_categories_by=None, show_counts=True)
            plt.savefig('{}.upSet.{}'.format(out_prefix, graph_format), dpi=300)
            plt.close()
    else:
        for group, name in zip(venn_list, venn_names):
            groups = group.split(',')
            tmp_dict = {x: y for x, y in venn_set_dict.items() if x in groups}
            if len(tmp_dict) > 1:
                plot(from_contents(tmp_dict), sum_over=False, sort_categories_by=None, show_counts=True)
                plt.savefig('{}.{}.upSet.{}'.format(out_prefix, name, graph_format), dpi=300)
                plt.close()


def tcr_venn(metadata, label_field, factor_field, outdir=None):
    seq_types = ['CDR3nt', 'CDR3aa']
    for seq_type in seq_types:
        set_dict, group = merge_clone_seq(metadata, label_field, factor_field, seq_type=seq_type, outdir=outdir)
        venn_plot(set_dict, group, out_prefix=os.path.join(outdir, seq_type))


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['tcr_venn'], log=False)
