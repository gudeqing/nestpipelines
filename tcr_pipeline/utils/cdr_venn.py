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
    seq_df = pd.DataFrame()
    for each_file, sample in zip(meta_df['files'], meta_df[label_field]):
        clone_info = pd.read_csv(each_file, sep='\t')
        seq_df[sample] = clone_info[seq_type]
    if outdir is None:
        outdir = os.getcwd()
    out = os.path.join(outdir, seq_type + '.merged.xls')
    seq_df.to_csv(out, sep='\t')
    group = os.path.join(outdir, f'{factor_field}.group.txt')
    meta_df[[label_field, factor_field]].to_csv(group, index=False, header=False, sep='\t')
    return out, group


def venn_plot(files: list, exp=None, out_prefix='result', has_header=False,
        intersect_only=True, intersect_xoy=1, union_only=False,
        set_names: list = None, venn_list: list = None, venn_names: list = None, graph_format='png'):
    """
    根据文件内容构建集合, 并按指定规则进行运算, 默认计算所有集合的交集
    :param files: 当仅提供一个文件时, 文件的各列被当作是集合, 集合的元素是单元格的内容;
    提供多个文件时, 每个文件内容被当作一个集合, 集合的元素为一整行。
    :param exp: 表达式, 字符串的形式, 如's1-s2'表示第一个集合减去第二个集合, 集合顺序与文件提供的顺序一一对应
    :param out_prefix: 指定集合运算结果的文件名前缀
    :param has_header: 指定文件是否包含header, 默认无, 如有header, header不参与计算
    :param intersect_only: 如提供, 不考虑exp指定的运算, 而是计算所有集合的交集, 即交集结果的所有元素在集合中出现的频数等于集合数
    :param intersect_xoy: 如提供, 不考虑exp指定的运算, 而是计算所有集合的交集, 而且输出交集结果的元素
    在所有集合中出现的频数大于或等于该参数指定的阈值.
    :param union_only: 计算各个集合的并集
    :param set_names: 用于画venn图, 对各个集合进行命名, 与文件名顺序应一致, 默认对文件名进行'.'分割获取第一个字符串作为集合名
    :param venn_list: 用于画venn图, 如 'A,B,C' 'B,C,D'表示画两个韦恩图, 第一个韦恩图用ABC集合, 而第二个韦恩图用BCD集合,
    默认None, 用所有集合画一个韦恩图; 另外, 可以给该参数输入一个文件, 第一列为集合名, 第二列为分组信息, 后续画图将按照此分组信息分别进行
    :param venn_names: 与venn_list一一对应, 用于分别命名venn图文件
    :param graph_format: output figure format, default png
    :return: None
    """
    venn_set_dict = dict()
    set_number = len(files)
    if len(files) >= 2:
        for ind, each in enumerate(files, start=1):
            exec('s{}=set(open("{}").readlines())'.format(ind, each))
            if set_names is None:
                name = os.path.basename(each).split('.', 1)[0]
                exec('venn_set_dict["{}"] = s{}'.format(name, ind))
            else:
                exec('venn_set_dict["{}"] = s{}'.format(set_names[ind - 1], ind))
    else:
        import pandas as pd
        table = pd.read_table(files[0], header=0 if has_header else None)
        set_number = table.shape[1]
        set_names = table.columns if set_names is None else set_names
        for i in range(table.shape[1]):
            exec('s{}=set(table.iloc[:, {}])'.format(i + 1, i))
            exec('venn_set_dict["{}"] = s{}'.format(set_names[i], i + 1))

    result = list()
    count_dict = dict()
    if exp:
        print("do as you say in exp")
        result = eval(exp)
    elif intersect_xoy > 1:
        print('do intersect_xoy')
        union = eval('|'.join(['s' + str(x) for x in range(1, set_number + 1)]))
        result = set()
        for each in union:
            varspace = dict(locals())
            in_times = sum(eval("each in s{}".format(x), varspace) for x in range(1, set_number + 1))
            if in_times >= intersect_xoy:
                result.add(each)
                count_dict[each] = in_times
    elif union_only:
        print('do union only')
        result = eval('|'.join(['s' + str(x) for x in range(1, set_number + 1)]))
    elif intersect_only:
        print('do intersect only')
        result = eval('&'.join(['s' + str(x) for x in range(1, set_number + 1)]))
    if not result:
        print('result is empty!')
    else:
        print('result size: {}'.format(len(result)))
    with open(out_prefix + '.list', 'w') as f:
        if not count_dict:
            _ = [f.write(x) for x in result]
        else:
            data = ([x, count_dict[x]] for x in result)
            _ = [f.write(x.strip() + '\t' + str(count_dict[x]) + '\n') for x in result]

    # plot venn
    if venn_list is None:
        if 2 <= len(venn_set_dict) <= 6:
            venn.venn(venn_set_dict, cmap="tab10")
            plt.savefig(out_prefix + f'.venn.{graph_format}')
    else:
        if len(venn_list) == 1 and ',' not in venn_list[0]:
            with open(venn_list[0]) as f:
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

        if venn_names is None:
            venn_names = []
            for group in venn_list:
                venn_names.append(group.replace(',', '-'))

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
            plt.savefig('all.upSet.{}'.format(graph_format), dpi=300)
            plt.close()
    else:
        for group, name in zip(venn_list, venn_names):
            groups = group.split(',')
            tmp_dict = {x: y for x, y in venn_set_dict.items() if x in groups}
            if len(tmp_dict) > 1:
                plot(from_contents(tmp_dict), sum_over=False, sort_categories_by=None, show_counts=True)
                plt.savefig('all.{}.upSet.{}'.format(name, graph_format), dpi=300)
                plt.close()


def tcr_venn(metadata, label_field, factor_field, outdir=None):
    seq_types = ['CDR3nt', 'CDR3aa']
    for seq_type in seq_types:
        data, group = merge_clone_seq(metadata, label_field, factor_field, seq_type=seq_type, outdir=outdir)
        venn_plot([data], out_prefix=os.path.join(outdir, seq_type), venn_list=[group])


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['tcr_venn'])
