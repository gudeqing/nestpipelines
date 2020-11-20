import os
import pandas as pd
import scipy.stats as stats
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import seaborn as sns
from shutil import copyfile


def merge_vj_matrix(file_list:list, column_sep='D', out='merged.counts.csv', group_info:str=None):
    table = pd.read_csv(file_list[0], index_col=[0, 1, 2, 3], header=0, sep=None, engine='python')
    for each in file_list[1:]:
        each_table = pd.read_csv(each, index_col=[0, 1, 2, 3], header=0)
        table = table.join(each_table, how='outer')
    table.columns = [x.strip().split(column_sep, 1)[0] for x in table.columns]
    table = table.fillna(0)
    table.index.name = 'SampleID'
    group_df = pd.read_csv(group_info, index_col=0, header=0, sep=None, engine='python')
    table = group_df.transpose().append(table, sort=False)
    # write out result
    table.to_csv(out)
    return table


def merge_metric_matrix(file_list:list, name_sep='_', out='merged.metric.csv',
                        group_info:str=None, new_name_col=None, factor_col=None):
    """
    合并从IR（https://10.62.2.16/ir/secure/home.html）到处得结果文件*metric.csv
    :param file_list:
    :param name_sep: 用于分割字符串提取出样本名，如metric文件名ZD909232002R1L1_135_RNA_v1_*.metrics.csv,
        则用'_'分隔得到ZD909232002R1L1作为样本名或建库编号
    :param out: 结果文件名
    :param group_info: 分组信息, 第一列必须是建库编号或样本名，第二列是新的样本名称，样本名称必须唯一，其他列可以是各种信息或分组信息
    :param new_name_col: 指定group_info中的一列, 作为样本别名，将放在第二列
    :param factor_col: 指定group_info中的一列, 作为factor列，也就是分组信息, 该列信息要带入到输出文件Report.diversity.summary.csv
    :return: 输出4个文件，有metric,count, frequency， metric.report
    """
    new_names = [
        'productive_reads', 'rescued_productive_reads', 'unproductive_reads',
        'off_target_reads', 'plus_strand_counts', 'minus_strand_counts'
    ]
    old_names = [
        'Productive_reads', 'Rescued_productive_reads', 'Unproductive_reads',
        'Off_target_reads', 'Plus_strand_counts', 'Minus_strand_counts',
    ]
    new2old = dict(zip(new_names, old_names))
    table = pd.read_csv(file_list[0], index_col=None, header=0, sep=None, engine='python')
    table.columns = [new2old[x.strip()] if x.strip() in new2old else x.strip() for x in table.columns]
    for each in file_list[1:]:
        each_table = pd.read_csv(each, index_col=0, header=0)
        each_table.columns = [new2old[x.strip()] if x.strip() in new2old else x.strip() for x in each_table.columns]
        table = table.append(each_table, sort=False)
    samples = [x.strip().split(name_sep, 1)[0] for x in table['Sample_Name']]
    # table = table.iloc[:, 4:]
    table.index = samples
    table.index.name = 'SampleID'
    # print(table.head())
    if group_info is not None:
        # for report
        os.makedirs('1.SampleInfo/', exist_ok=True)
        # copyfile(group_info, '1.SampleInfo/sample.info.txt')
        if group_info.endswith('xlsx'):
            group_df = pd.read_excel(group_info, index_col=0, header=0)
        else:
            group_df = pd.read_csv(group_info, index_col=0, header=0, sep=None, engine='python')
        group_df.to_csv('1.SampleInfo/sample.info.txt', sep='\t')
        samples_found = group_df.index.intersection(table.index)
        print('samples:', len(samples_found))
        print('these samples has no sample information', table.index.difference(group_df.index))
        group_df = group_df.loc[samples_found]
        table = group_df.join(table, how='right', sort=False)
        table = table.loc[list(group_df.index)]
    table.columns = [x.strip() for x in table.columns]
    table.to_csv(out, encoding='utf_8_sig')
    # metric for report
    if new_name_col is not None:
        target_cols = [
            new_name_col, 'Productive_reads', 'Rescued_productive_reads', 'Unproductive_reads',
            'Off_target_reads', 'Plus_strand_counts', 'Minus_strand_counts',
        ]
    else:
        target_cols = [
            'Productive_reads', 'Rescued_productive_reads', 'Unproductive_reads',
            'Off_target_reads', 'Plus_strand_counts', 'Minus_strand_counts',
        ]
    # print(list(table.columns))
    os.makedirs('2.SampleQC', exist_ok=True)
    outfile = os.path.join('2.SampleQC', 'Report.qc.summary.csv')
    table.loc[:, target_cols].to_csv(outfile)

    # clone summary for report
    target_cols = []
    if factor_col is not None:
        target_cols += [factor_col]
    if new_name_col is not None:
        target_cols += [new_name_col]
    target_cols += [
        'Clones', 'Evenness', 'Reads', 'Shannon_diversity',
        'convergent_TCR_frequency', 'clone_gini_index'
    ]
    # target_cols += [x for x in table.columns if x.endswith('downsample')]
    os.makedirs('1.DiversitySummary', exist_ok=True)
    outfile = os.path.join('1.DiversitySummary', 'Report.diversity.summary.csv')
    table.loc[:, target_cols].to_csv(outfile)

    # set new_name as index
    if new_name_col is not None:
        table.set_index(new_name_col, inplace=True)
    # extract_vj_frequency matrix
    vj_freq_cols = [x for x in table.columns if x.endswith('_frequency') and x.startswith('TR')]
    freq_data = table[vj_freq_cols]
    freq_data = freq_data.transpose()
    freq_data.index.name = 'Gene'
    freq_data.index = [x.split('_frequency')[0] for x in freq_data.index]
    freq_data.to_csv(out.rsplit('.', 1)[0]+'.frequency.csv')
    # _counts
    count_cols = [x for x in table.columns if x.endswith('_counts') and x.startswith('TR')]
    count_data = table[count_cols]
    count_data = count_data.transpose()
    count_data.index.name = 'Gene'
    count_data.index = [x.split('_counts')[0] for x in count_data.index]
    count_data.to_csv(out.rsplit('.', 1)[0] + '.count.csv')
    return table


def violin_plot(df, data_col, group_cols:list, hue_cols:list=None, index_col=None, split=False, orient=None,
                exchange_xy=False, out=None, scale='width', style='darkgrid', target_index=None, inner=None):
    sns.set(style=style)
    sns.set(font_scale=0.5)
    init_inner = inner
    # output name
    if out is None:
        if hue_cols is not None:
            out = '{data_col}.{group_cols}.{hue_cols}.violin.pdf'.format(
                data_col=data_col, group_cols='_'.join(group_cols), hue_cols='_'.join(hue_cols)
            )
        else:
            out = '{data_col}.{group_cols}.violin.pdf'.format(
                data_col=data_col, group_cols='_'.join(group_cols)
            )
    # prepare data
    if type(df) is str:
        data = pd.read_csv(df, index_col=0, header=0, sep=None, engine='python')
    else:
        data = df
    if index_col is not None:
        data.set_index(index_col, inplace=True)
    if target_index is not None:
        targets = [x.strip().split()[0] for x in open(target_index)]
        data = data.loc[targets]
    data.columns = [x.strip() for x in data.columns]
    # plot
    if hue_cols is None:
        hue_cols = [None] * len(group_cols)
    if len(group_cols) > 1:
        fig, axes = plt.subplots(len(group_cols), 1)
        for ind, group_col, hue_col in zip(range(len(group_cols)), group_cols, hue_cols):
            print(data.groupby(group_col).size())
            print(data.groupby(group_col).size().mean())
            if init_inner is None:
                if data.groupby(group_col).size().mean() > 50:
                    inner = 'quartile'
                else:
                    inner = 'stick'
            if hue_col is not None:
                if hue_col.lower() == 'none':
                    hue_col = None
            if not exchange_xy:
                ax = sns.violinplot(x=group_col, y=data_col, data=data, hue=hue_col, scale=scale, orient=orient,
                                    width=0.8, linewidth=0.5, inner=inner, split=split, ax=axes[ind])
            else:
                ax = sns.violinplot(x=data_col, y=group_col, data=data, hue=hue_col, scale=scale, orient=orient,
                                    width=0.8, linewidth=0.5, inner=inner, split=split, ax=axes[ind])
            plt.setp(ax.collections, linewidth=0.3)
    else:
        for ind, group_col, hue_col in zip(range(len(group_cols)), group_cols, hue_cols):
            if inner is None:
                if data.groupby(group_col).size().mean() > 50:
                    inner = 'quartile'
                else:
                    inner = 'stick'
            if hue_col is not None:
                if hue_col.lower() == 'none':
                    hue_col = None
            if not exchange_xy:
                ax = sns.violinplot(x=group_col, y=data_col, data=data, hue=hue_col, scale=scale, orient=orient,
                                    width=0.8, linewidth=0.5, inner=inner, split=split)
            else:
                ax = sns.violinplot(x=data_col, y=group_col, data=data, hue=hue_col, scale=scale, orient=orient,
                                    width=0.8, linewidth=0.5, inner=inner, split=split)
            plt.setp(ax.collections, linewidth=0.3)

    plt.savefig(out, dpi=300)


def scatter(df, data_col, group_cols:list, hue_cols:list=None, style_cols:list=None, index_col=None, exchange_xy=False,
            plot_style='whitegrid', target_index=None, out=None, legend='brief', alpha:float=None):
    sns.set_style(style=plot_style)
    sns.set(font_scale=0.5)
    # output name
    if out is None:
        out_name = 'scatter.' + data_col
        out_name += '.' + '-'.join(group_cols)
        if hue_cols is not None:
            out_name += '.' + '-'.join(hue_cols)
        if style_cols is not None:
            out_name += '.' + '-'.join(style_cols)
        out_name += '.pdf'
        print('Output:', out_name)
    else:
        out_name = out
    # prepare data
    if type(df) is str:
        data = pd.read_csv(df, index_col=0, header=0, sep=None, engine='python')
    else:
        data = df
    if index_col is not None:
        data.set_index(index_col, inplace=True)
    if target_index is not None:
        targets = [x.strip().split()[0] for x in open(target_index)]
        data = data.loc[targets]
    data.columns = [x.strip() for x in data.columns]
    # plot
    if hue_cols is None:
        hue_cols = [None] * len(group_cols)
    if style_cols is None:
        style_cols = [None] * len(group_cols)

    axes = None
    if len(group_cols) > 1:
        fig, axes = plt.subplots(len(group_cols), 1)
    if alpha is None:
        alpha='auto'

    for ind, group_col, hue_col, style_col in zip(range(len(group_cols)), group_cols, hue_cols, style_cols):
        print(data.groupby(group_col).size())
        print(data.groupby(group_col).size().mean())
        if hue_col is not None:
            if hue_col.lower() == 'none':
                hue_col = None
        if style_col is not None:
            if style_col.lower() == 'none':
                style_col = None
        if axes is not None:
            tmp_ax = axes[ind]
        else:
            tmp_ax = None
        if not exchange_xy:
            x_data = group_col
            y_data = data_col
        else:
            y_data = group_col
            x_data = data_col
        ax = sns.scatterplot(x=x_data, y=y_data, data=data, hue=hue_col, style=style_col, ax=tmp_ax,
                             legend=legend, alpha=alpha)
        plt.setp(ax.collections, linewidth=0.3)

    plt.savefig(out_name, dpi=300)


def diff_diversity_test(df, data_col, group_df, cmp_info, index_col=None, target_index=None, paired=False, out=None):
    # prepare data
    if type(df) is str:
        data = pd.read_csv(df, index_col=0, header=0, sep=None, engine='python')
    else:
        data = df
    if index_col is not None:
        data.set_index(index_col, inplace=True)
    if target_index is not None:
        targets = [x.strip().split()[0] for x in open(target_index)]
        data = data.loc[targets]
    data.columns = [x.strip() for x in data.columns]

    # get group information
    if type(group_df) == str:
        group_df = pd.read_csv(group_df, header=0, index_col=0, sep=None, engine='python')
    group_dict = dict()
    for sample, cols in group_df.iterrows():
        for group_name in cols:
            group_dict.setdefault(group_name, list())
            group_dict[group_name].append(sample)
    # print(list(group_dict.keys()))

    # get compare information
    cmp_list = [x.strip().split()[:2] for x in open(cmp_info)]
    print('compare_info:', cmp_list)
    # check group if matches compare info
    for cmp_groups in cmp_list:
        for group in cmp_groups:
            if group not in group_dict:
                raise Exception(f'{group} is not defined in group info!')
    # test
    if out is None:
        out = 'test.result.txt'
    with open(out, 'w') as f:
        f.write('ctrl_group\ttest_group\ttest_method\tpvalue\tstatistic\tctrl_size\ttest_size\tctrl_data\ttest_data\n')
        for ctrl, test in cmp_list:
            ctrl_samples = [x for x in group_dict[ctrl] if x in data.index]
            test_samples = [x for x in group_dict[test] if x in data.index]
            ctrl_data = data.loc[ctrl_samples, data_col]
            test_data = data.loc[test_samples, data_col]
            if len(ctrl_samples) < 2 and len(test_samples) < 2:
                print(f'ctrl or test group has too few samples for testing, skip {ctrl} vs {test}')
                continue
            ctrl_data_str = ','.join(str(round(x, 3)) for x in ctrl_data)
            test_data_str = ','.join(str(round(x, 3)) for x in test_data)
            data_detail = f'{len(ctrl_data)}\t{len(test_data)}\t{ctrl_data_str}\t{test_data_str}'
            statistic, pvalue = stats.ranksums(ctrl_data, test_data)
            f.write(f'{ctrl}\t{test}\tWilcoxon_rank_sum\t{pvalue}\t{statistic}\t{data_detail}\n')
            statistic, pvalue = stats.mannwhitneyu(ctrl_data, test_data)
            f.write(f'{ctrl}\t{test}\tMann_Whitney_U\t{pvalue}\t{statistic}\t{data_detail}\n')
            if paired:
                if len(ctrl_data) == len(test_data):
                    statistic, pvalue = stats.wilcoxon(ctrl_data, test_data)
                f.write(f'{ctrl}\t{test}\t Wilcoxon_signed_rank\t{pvalue}\t{statistic}\t{data_detail}\n')
    df = pd.read_csv(out, header=0, index_col=None, sep='\t')
    df.to_excel(out.rsplit('.', 1)[0]+'.xlsx', index=False)


def mean_duplicated(df, dup_col:list, target_index=None, index_col=None, out=None):
    if type(df) is str:
        data = pd.read_csv(df, index_col=0, header=0, sep=None, engine='python')
    else:
        data = df
    if index_col is not None:
        data.set_index(index_col, inplace=True)
    if target_index is not None:
        targets = [x.strip().split()[0] for x in open(target_index)]
        data = data.loc[targets]
    data.columns = [x.strip() for x in data.columns]
    dup_id = data[dup_col[0]].apply(str)
    if len(dup_col) >= 2:
        for each in dup_col[1:]:
            dup_id += '|' + data[each].apply(str)
    data['dup_id'] = dup_id
    data = data.groupby('dup_id').mean()
    if len(dup_col) > 1:
        group_info = pd.DataFrame([x.split('|') for x in data.index], columns=dup_col, index=data.index)
        data = group_info.join(data)
    data.index = [x.replace('|', '') for x in data.index]
    if out is None:
        out = 'uniq.{}.mean.csv'.format('-'.join(dup_col))
    data.to_csv(out)
    data.to_excel(out.rsplit('.', 1)[0] + '.xlsx')


def convert2vdjtools(files:list, out_dir='Vdjtools_input', group_info=None):
    """
    输入从IR导出的*clone_summary.txt文件
    :param files: *clone_summary.txt
    :param out_dir: 输出结果的绝对路径
    :param group_info: 样本分组信息，第一列为样本id，需和推断的os.path.basename(each).split('_', 1)[0]一致。
    :return:
    """
    samples = list()
    out_path = list()
    clone_detail_dir = os.path.join(out_dir, '3.CloneDetail')
    os.makedirs(clone_detail_dir, exist_ok=True)
    for each in files:
        sample = os.path.basename(each).split('_', 1)[0]
        copyfile(each, os.path.join(clone_detail_dir, f'{sample}.clone_summary.csv'))
        # samples.append(sample[:-4])
        samples.append(sample)
        out_name = os.path.join(out_dir, '{}.clone_summary.txt'.format(sample))
        out_path.append(out_name)
        table = pd.read_csv(each, index_col=None, header=0, sep=None, engine='python')
        target_cols = ['Total Counts', 'Frequency', 'CDR3 NT', 'CDR3 AA', 'Variable', 'Diversity', 'Joining']
        new_cols = ['count', 'frequency', 'CDR3nt', 'CDR3aa', 'V', 'D', 'J']
        table = table[target_cols]
        table.columns = new_cols
        table.to_csv(out_name, sep='\t', index=False)
    # make matadata
    if group_info:
        if group_info.endswith('xlsx'):
            group_df = pd.read_excel(group_info, index_col=0, header=0)
        else:
            group_df = pd.read_csv(group_info, index_col=0, header=0, sep=None, engine='python')
        # group_df = group_df.applymap(lambda x: x.replace(' ', '_'))
        if len(set(samples) & set(group_df.index)) < 1:
            print(samples)
            print(set(group_df.index))
            raise Exception('sample id does not match!')
        ori_order = group_df.index
        group_df = group_df.loc[samples]
        group_df['files'] = out_path
        group_df['sample_id'] = group_df.index
        cols = ['files', 'sample_id'] + [x for x in group_df.columns if x not in ['files', 'sample_id']]
        group_df = group_df[cols]
        group_df = group_df.loc[[x for x in ori_order if x in samples]]
        group_df.to_csv('metadata.txt', sep='\t', index=False)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), exclude=['pd'])

