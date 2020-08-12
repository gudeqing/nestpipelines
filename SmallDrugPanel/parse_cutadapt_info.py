import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib._color_data as mcd
import matplotlib.patches as mpatch
import pandas as pd
import itertools
from cycler import cycler


def parse_info_file(file, out_prefix='primer'):
    no_primer = 0
    prefix_primer = 0
    middle_primer = 0
    # perfect_primer仅允许一个错配
    perfect_primer = 0
    # 错误>=4 或者和预期primer长度相差超过4bp
    dubious_primer = 0
    read_number = 0
    primer_dict = dict()
    umi_dict = dict()
    with open(file) as f:
        for i, line in enumerate(f):
            lst = line.split('\t')
            read_name = lst[0]
            errors = int(lst[1])
            if errors == -1:
                no_primer += 1
            else:
                start = int(lst[2])
                end = int(lst[3])
                primer_name = lst[7]
                primer_len = int(primer_name.split(':', 3)[2].split('=')[1])
                primer_dict.setdefault(primer_name, 0)
                primer_dict[primer_name] += 1
                # umi = read_name.split()[0].split(':')[-1][:11]
                umi = read_name.split()[0].split(':')[-1]
                umi_dict.setdefault(primer_name, set())
                umi_dict[primer_name].add(umi)

                if start == 0:
                    prefix_primer += 1
                else:
                    middle_primer += 1
                if errors <= 1 and primer_len == (end - start):
                    perfect_primer += 1
                else:
                    if errors >= 4 or abs(primer_len - end + start) >= 4:
                        dubious_primer += 1
        else:
            read_number = i + 1

    primer_number = len(primer_dict)
    has_primer_number = read_number - no_primer
    # stat for each primer
    # chr2:29446286-29446312:len=27:strand=+:gene=ALK
    order = lambda x: (x.split(':')[-1], x.split('strand=')[1].split(':')[0])
    keys = sorted(primer_dict.keys(), key=order)
    primer_dict = {k: primer_dict[k] for k in keys}
    umi_dict = {k: len(umi_dict[k]) for k in keys}
    # data = pd.DataFrame({'read_count': primer_dict, 'umi_count': umi_dict})
    # data['read_count/read_number'] = data['read_count']/read_number
    # data['umi_count/pair_number'] = data['umi_count']/read_number
    # data.index.name = 'primer'
    # data = data.round(3)
    with open(out_prefix+'.stats.txt', 'w') as f:
        f.write(f'pair_number\t{read_number}\n')
        f.write(f'primer_number\t{primer_number}\n')
        f.write(f'has_primer\t{has_primer_number}\t{has_primer_number/read_number:.3%}\n')
        f.write(f'no_primer\t{no_primer}\t{no_primer/read_number:.3%}\n')
        f.write(f'perfect_primer\t{perfect_primer}\t{perfect_primer/read_number:.3%}\n')
        f.write(f'dubious_primer\t{dubious_primer}\t{dubious_primer/read_number:.3%}\n')
        f.write(f'prefix_primer\t{prefix_primer}\t{prefix_primer/read_number:.3%}\n')
        f.write(f'middle_primer\t{middle_primer}\t{middle_primer/read_number:.3%}\n')
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
    new_key = lambda x: x.split('strand=')[1].replace(':gene=', '')
    genes = sorted({new_key(x) for x in primer_dict.keys()})
    colors = get_color_pool(len(genes))
    color_dict = dict(zip(genes, colors))
    for info_dict, type_, out, out2 in zip([primer_dict, umi_dict], ['reads', 'UMIs'], [out_name, out_name2], [out_name3, out_name4]):
        # plot for each primer
        colors = [color_dict[new_key(x)] for x in info_dict.keys()]
        title = f'Number of {type_} containing primer'
        plot_bar(info_dict, colors=colors, fontsize=3, rotation=90, out=out, title=title)

        # plot for each gene
        gene_stat_dict = dict()
        for k, read_count in info_dict.items():
            key = new_key(k)
            gene_stat_dict.setdefault(key, 0)
            gene_stat_dict[key] += read_count

        colors = [color_dict[x.split(':')[-1]] for x in gene_stat_dict.keys()]
        title = f'Number of {type_} along with primer for each gene'
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
    mean_number = sum(info_dict.values()) / len(info_dict)
    ax.axhline(y=mean_number, c="r", ls="--", lw=0.6, label=f'Mean={int(mean_number)}')
    ax.axhline(y=mean_number * 0.2, c="r", ls="--", lw=0.5, label='20%Mean')
    ax.axhline(y=mean_number * 0.3, c="r", ls="--", lw=0.5, label='30%Mean')
    ax.spines['right'].set_color('None')  # 右框不显示
    ax.spines['top'].set_color('None')  # 上框不显示
    if label_bar:
        def autolabel(rects):
            """Attach a text label above each bar in *rects*, displaying its height."""
            for rect in rects:
                height = rect.get_height()
                ax.annotate('{}'.format(height),
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

