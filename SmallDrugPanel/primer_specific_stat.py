import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam


"""
1.fastq处理时，把primer信息如在基因组上的位置信息带入read id
2.从bam信息中提取每一个read1比对到参考基因组的位置：
    a.比对位置和primer位置进行比较，假设primer的位置：chrX:a-b, strand=S1, read1的比对位置chrY:c-d, strand=S2
        if chrX != chrY or S1 != S2 or Un-mapped:
            primer 非特异性 + 1
        else:
            if S1 = +：
                if abs(c-b) > 15:
                    primer 非特异性 + 1
                else:
                    primer 特异性 + 1
            else:
                if abs(a-d) > 15:
                    primer 非特异性 + 1
                else:
                    primer 特异性 + 1
        当read1可以比对到多个位置，进统计primary alignment，
        这个必须在通读bam完成后才能知晓，除非bam用readName排序
        
3.如果使用umitools dedup后的bam统计，会有什么影响呢：
    1. primer的特异性结合分析结果是否更加可靠？应该是，因为排除了pcr的非均一性可能带来的偏差，dedup还原了扩增前的状况
    2. dedup后的bam不会包含ummapped的read，这可能导致分析结果存在一定的偏差，即特异性偏高
4. 基于3的考虑，直接采用dedup后的bam统计
5. 这里的特异性占比可以理解为mapping后的 off-target-ratio
6. 这里的统计仅针对genePrimer的设计有效
"""


def stat(bam, out_prefix):
    specific = dict()
    unspecific = dict()
    bam_obj = pysam.AlignmentFile(bam, "rb")
    for read in bam_obj.fetch():
        if (not read.is_read1) or read.is_secondary:
            continue
        read_name = read.query_name
        primer_lst = read_name.split(':', 5)[:5]
        primer = ':'.join(primer_lst)
        umi = read_name.rsplit(':', 1)[1]
        map_strand = '-' if read.is_reverse else '+'

        target = False
        if primer_lst[0] == read.reference_name and (map_strand == primer_lst[3]):
            ps, pe = primer_lst[1].split('-')
            if map_strand == '+':
                if abs(int(pe) - read.reference_start) < 15:
                    target = True
            else:
                if read.reference_end is not None:
                    if abs(int(ps) - read.reference_end) < 15:
                        target = True

        specific.setdefault(primer, set())
        unspecific.setdefault(primer, set())
        if target:
            specific[primer].add(umi)
        else:
            unspecific[primer].add(umi)
    else:
        bam_obj.close()
    print('Primer Number:', len(specific))
    order = lambda x: (x.split(':')[-1], x.split(':')[-2])
    keys = sorted(specific.keys(), key=order)
    specific = {k: len(specific[k]) for k in keys}
    unspecific = {k: len(unspecific[k]) for k in keys}

    data = pd.DataFrame({'on_target': specific, 'off_target': unspecific})
    data['on_target_ratio'] = data['on_target']/data.sum(axis=1)
    data.index.name = 'primer'
    data = data.round(3).to_csv(f'{out_prefix}.PrimerUMI.csv')

    # plot
    get_gene = lambda x: x.split(':')[-1]
    genes = sorted({get_gene(x) for x in specific.keys()})
    colors = get_color_pool(len(genes))
    color_dict = dict(zip(genes, colors))
    colors = [color_dict[get_gene(x)] for x in specific.keys()]
    out = f'{out_prefix}.PrimerUMI.bar.pdf'
    stacked_bar(specific, unspecific, colors, colors2='gray', fontsize=3, rotation=90, out=out)


def stacked_bar(info_dict, info_dict2, colors, colors2, fontsize, rotation, out, label_bar=False, title=''):
    fig, ax = plt.subplots()
    x_num = len(info_dict)
    x_pos = range(x_num)
    rects = ax.bar(x_pos, info_dict.values(), width=0.7, color=colors)
    matplotlib.rcParams['hatch.linewidth'] = 0.1
    rects2 = ax.bar(x_pos, info_dict2.values(), bottom=list(info_dict.values()),
                    width=0.7, color=colors2, label='Off Target', hatch="///")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(info_dict.keys(), fontsize=fontsize, rotation=rotation)
    ax.set_xlim(-1, x_num)
    ax.tick_params(axis='y', labelsize=fontsize+1)
    # mean_number = sum(info_dict.values()) / len(info_dict)
    mean_number = np.mean(list(info_dict.values()))
    median_number = np.median(list(info_dict.values()))
    # ax.axhline(y=mean_number, c="r", ls="--", lw=0.6, label=f'Mean={int(mean_number)}')
    ax.axhline(y=median_number, c="r", ls="--", lw=0.6, label=f'Median={int(median_number)}')
    # ax.axhline(y=mean_number * 0.25, c="r", ls="--", lw=0.5, label='25%Mean')
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
                            xytext=(0, -3) ,  # 3 points vertical offset
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
    xcmds.xcmds(locals(), include=['stat'])




