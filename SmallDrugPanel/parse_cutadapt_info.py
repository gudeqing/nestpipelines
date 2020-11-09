import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
# import matplotlib._color_data as mcd
# import matplotlib.patches as mpatch
# import pandas as pd
# import itertools
# from cycler import cycler


def parse_info_file(file, out_prefix='primer', primer_fasta=None):
    no_primer = 0
    prefix_primer = 0
    middle_primer = 0
    # perfect_primer仅允许2个错配, 且长度一致相差在3bp以内
    perfect_primer = 0
    # 错误>=3 或者和预期primer长度相差超过3bp
    dubious_primer = 0
    read_number = 0
    primer_dict = dict()
    umi_dict = dict()
    # 统计primer测序错误的情况
    base_num = 0
    total_error = 0
    total_error_primer = 0  # 用于计算error的总primer数
    transition = dict()
    transition_base_num = 1
    umi_set = set()
    if primer_fasta:
        fasta = open(primer_fasta).readlines()
        names = [x[1:].strip() for x in fasta[0::2]]
        fa = [x.strip().strip('^').upper() for x in fasta[1::2]]
        primer_fa_dict = dict(zip(names, fa))
        # print(primer_fa_dict)
    else:
        primer_fa_dict = dict()

    with open(file) as f:
        for i, line in enumerate(f):
            read_number += 1
            lst = line.split('\t')
            read_name = lst[0]
            # if '/' in read_name[-2:]:
            #     is_read1 = read_name[-1] == '1'
            # else:
            #     is_read1 = read_name.split()[1][0] == '1'
            # if not is_read1:
            #     continue
            # 由于cutadapt默认只输出read1包含的primer匹配情况，无需上面的判断

            errors = int(lst[1])
            if errors == -1:
                no_primer += 1
                # print(line)
            else:
                start = int(lst[2])
                end = int(lst[3])
                umi = read_name.split()[0].split(':')[-1]
                primer_name = lst[7]
                primer_len = int(primer_name.split(':', 3)[2])
                if start == 0:
                    prefix_primer += 1
                else:
                    middle_primer += 1
                if errors <= 2 and abs(primer_len - (end - start)) <= 3:
                    perfect_primer += 1
                else:
                    # if errors >= 3 or abs(primer_len - end + start) >= 4:
                    dubious_primer += 1

                # 统计primer测序错误的情况
                if umi not in umi_set:
                    umi_set.add(umi)
                    base_num += primer_len  # 只要包含primer，则都将进入统计作为分母
                    if errors > 0:
                        total_error += errors
                        total_error_primer += 1
                    # 统计碱基之间的转换情况. 尽量排除indel带来的影响，所以需要下面的if
                    if primer_len == (end - start) and primer_fasta:
                        transition_base_num += primer_len
                        primer_fa = primer_fa_dict[primer_name]
                        # if primer_fa != lst[5].upper():
                        #     print(primer_fa)
                        #     print(lst[5])
                        for ref, alt in zip(primer_fa, lst[5].upper()):
                            if ref != alt:
                                transition.setdefault((ref, alt), 0)
                                transition[(ref, alt)] += 1

                # update stat dict
                primer_dict.setdefault(primer_name, 0)
                primer_dict[primer_name] += 1  # 相当于统计每个primer包含多少条reads
                umi_dict.setdefault(primer_name, set())
                umi_dict[primer_name].add(umi)

    # print seq error info
    print('统计时，只考虑UMI第一次出现所对应的那条reads, 其他包含相同UMI的read不进入统计，目的是要估计第一轮PCR导致的错误率')
    error_dict = dict()
    error_dict['overall_error_rate'] = total_error/base_num
    total_trans_times = sum(transition.values())
    overall_trans_rate = total_trans_times/transition_base_num
    overall_error_rate = total_error/base_num
    print('Read1 数量为', read_number)
    print('包含primer的read1数量为', read_number - no_primer)
    print('包含primer的read1所包含的Uniq UMI数量为', len(umi_set))
    print('平均每条primer包含碱基错误数量为', total_error/len(umi_set))
    print(f'primer总体错误率={total_error}/{base_num}={overall_error_rate:.4%}')
    print(f'Overall transition/transversion rate {overall_trans_rate:.4%}')
    print('转换类型: 转换发生率 | 转换类型占比率')
    for k in sorted(transition.keys()):
        num = transition[k]
        new_key = f'{k[0]}>{k[1]}'
        trans_rate = num/total_trans_times
        print(f'{new_key}: {num/transition_base_num:.4%} | {trans_rate:.4%}')
        error_dict[new_key] = trans_rate

    with open(out_prefix+'.PrimerErrorStat.txt', 'w') as f:
        for k in sorted(error_dict.keys()):
            v = error_dict[k]
            f.write(f'{k}\t{v}\n')

    primer_number = len(primer_dict)
    has_primer_number = read_number - no_primer
    # stat for each primer
    # chr2:29446286-29446312:len=27:strand=+:gene=ALK
    # order = lambda x: (x.split(':')[-1], x.split('strand=')[1].split(':')[0])
    order = lambda x: (x.split(':')[-1], x.split(':')[-2])
    keys = sorted(primer_dict.keys(), key=order)
    primer_dict = {k: primer_dict[k] for k in keys}
    umi_dict = {k: len(umi_dict[k]) for k in keys}
    # data = pd.DataFrame({'read_count': primer_dict, 'umi_count': umi_dict})
    # data['read_count/read_number'] = data['read_count']/read_number
    # data['umi_count/pair_number'] = data['umi_count']/read_number
    # data.index.name = 'primer'
    # data = data.round(3)
    with open(out_prefix+'.primer_stats.txt', 'w') as f:
        f.write(f'pair_number\t{read_number}\n')
        f.write(f'primer_number\t{primer_number}\n')
        f.write(f'has_primer\t{has_primer_number}\t{has_primer_number/read_number:.3%}\n')
        f.write(f'no_primer\t{no_primer}\t{no_primer/read_number:.3%}\n')
        f.write(f'perfect_primer\t{perfect_primer}\t{perfect_primer/read_number:.3%}\n')
        f.write(f'dubious_primer\t{dubious_primer}\t{dubious_primer/read_number:.3%}\n')
        f.write(f'prefix_primer\t{prefix_primer}\t{prefix_primer/read_number:.3%}\n')
        f.write(f'middle_primer\t{middle_primer}\t{middle_primer/read_number:.3%}\n')
        f.write(f'overall_primer_base_seq_error_rate\t{total_error}/{base_num}={total_error/base_num:.2%}\n')
        f.write(f'#primer\tread_count\tread_count/pair_number\tumi_count\tumi_count/pair_number\n')
        for k in keys:
            read_count = primer_dict[k]
            umi_count = umi_dict[k]
            read_ratio = read_count/read_number
            umi_ratio = umi_count/read_number
            f.write(f'{k}\t{read_count}\t{read_ratio:.3%}\t{umi_count}\t{umi_ratio:.3%}\n')

    # plot for each primer
    out_name = f'{out_prefix}.PrimerReadCounts.pdf'
    out_name2 = f'{out_prefix}.PrimerUMICounts.pdf'
    out_name3 = f'{out_prefix}.GenePrimerReadCounts.pdf'
    out_name4 = f'{out_prefix}.GenePrimerUMICounts.pdf'
    # new_key = lambda x: x.split('strand=')[1].replace(':gene=', '')
    get_gene = lambda x: x.split(':')[-1]
    genes = sorted({get_gene(x) for x in primer_dict.keys()})
    colors = get_color_pool(len(genes))
    color_dict = dict(zip(genes, colors))
    for info_dict, type_, out, out2 in zip([primer_dict, umi_dict], ['Read', 'UMI'], [out_name, out_name2], [out_name3, out_name4]):
        # plot for each primer
        colors = [color_dict[get_gene(x)] for x in info_dict.keys()]
        title = f'{type_} distribution for each primer'
        out = f'{out_prefix}.' + title.title().replace(' ', '') + '.pdf'
        plot_bar(info_dict, colors=colors, fontsize=3, rotation=90, out=out, title=title)

        # plot for each gene
        gene_stat_dict = dict()
        for k, read_count in info_dict.items():
            key = k.split(':', 3)[-1]
            gene_stat_dict.setdefault(key, 0)
            gene_stat_dict[key] += read_count

        colors = [color_dict[get_gene(x)] for x in gene_stat_dict.keys()]
        title = f'{type_} distribution for each gene'
        out2 = f'{out_prefix}.' + title.title().replace(' ', '') + '.pdf'
        plot_bar(gene_stat_dict, colors=colors, fontsize=5, rotation=90, out=out2, label_bar=True, title=title)


def plot_bar(info_dict, colors, fontsize, rotation, out, label_bar=False, title=''):
    fig, ax = plt.subplots()
    x_num = len(info_dict)
    x_pos = range(x_num)
    rects = ax.bar(x_pos, info_dict.values(), width=0.7, color=colors)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(info_dict.keys(), fontsize=fontsize, rotation=rotation)
    ax.set_xlim(-1, x_num)
    ax.tick_params(axis='y', labelsize=fontsize+1)
    # mean_number = sum(info_dict.values()) / len(info_dict)
    mean_number = np.mean(list(info_dict.values()))
    median_number = np.median(list(info_dict.values()))
    ax.axhline(y=mean_number, c="r", ls="--", lw=0.6, label=f'Mean={int(mean_number)}')
    ax.axhline(y=median_number, c="r", ls="--", lw=0.6, label=f'Median={int(median_number)}')
    ax.axhline(y=mean_number * 0.25, c="r", ls="--", lw=0.5, label='25%Mean')
    ax.axhline(y=median_number * 0.25, c="r", ls="--", lw=0.5, label='25%Median')
    ax.spines['right'].set_color('None')  # 右框不显示
    ax.spines['top'].set_color('None')  # 上框不显示
    if label_bar:
        def autolabel(rects):
            """Attach a text label above each bar in *rects*, displaying its height."""
            for rect in rects:
                height = rect.get_height()
                height_percent = f'{height/sum(info_dict.values()):.2%}'
                ax.annotate('{}'.format(height_percent),
                            xy=(rect.get_x() + rect.get_width() / 2, height),
                            xytext=(0, 2),  # 3 points vertical offset
                            textcoords="offset points",
                            ha='center', va='bottom', rotation=60, fontsize=5)
        autolabel(rects)
    ax.set_title(title, fontsize='small')
    plt.tight_layout()
    plt.legend(prop=dict(size=6))
    plt.savefig(out)
    plt.close()


def get_color_pool(n):
    # https://plot.ly/ipython-notebooks/color-scales/
    # import colorlover
    # if n <= 8:
    #     if n <= 3:
    #         n = 3
    #     return colorlover.scales[str(n)]['qual']['Set1']
    # if n <= 12:
    #     return colorlover.scales[str(n)]['qual']['Paired']
    # colorlover的颜色现在不能直接输给matplotlib，需要进行转换才可以

    import random
    random.seed(666)

    def get_random_color(pastel_factor=0.5):
        return [(x + pastel_factor) / (1.0 + pastel_factor) for x in [random.uniform(0, 1.0) for i in [1, 2, 3]]]

    def color_distance(c1, c2):
        return sum([abs(x[0] - x[1]) for x in zip(c1, c2)])

    def generate_new_color(existing_colors, pastel_factor=0.5):
        max_distance = None
        best_color = None
        for i in range(0, 100):
            color = get_random_color(pastel_factor=pastel_factor)
            # exclude some colors
            if np.absolute(np.array(color) - np.array([1, 1, 1])).sum() < 0.1:
                continue
            if not existing_colors:
                return color
            best_distance = min([color_distance(color, c) for c in existing_colors])
            if not max_distance or best_distance > max_distance:
                max_distance = best_distance
                best_color = color
        return best_color

    color_pool = []
    for i in range(0, n):
        color_pool.append(generate_new_color(color_pool, pastel_factor=0.3))
    # color_pool = [(int(x * 255), int(y * 255), int(z * 255)) for x, y, z in color_pool]
    color_pool = sorted(color_pool, key=lambda x: (x[0], x[1], x[2]))
    # color_pool = [matplotlib.colors.to_rgba(x) for x in color_pool]
    # print(color_pool)
    return color_pool


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['parse_info_file'])

