import json
from collections import Counter
import pandas as pd
import scipy.stats as stats
import pysam
import numpy as np


def get_seq_qual(contig, start, end, bam, min_bq=15):
    cols = bam.pileup(
        contig, start, end,
        stepper='samtools',
        truncate=True,
        min_base_quality=min_bq,
        ignore_orphans=False,
        ignore_overlaps=False,
    )
    seq_qual = [(col.get_query_sequences(), col.get_query_qualities()) for col in cols]
    return seq_qual


def get_center_seq(contig, pos, genome, size=1):
    return genome.fetch(contig, pos-size, pos+size+1)


def reverse_complement(seq):
    """
    :param seq: 输入序列
    :return: 返回reversed的序列
    """
    seq = seq.upper()
    complement = dict(zip(list("ATCG"), list("TAGC")))
    return ''.join(complement[base] if base in complement else base for base in seq[::-1])


def run(bed, bam, prefix, center_size=1, only_alt_site=False, exclude_sites=None,
        genome='/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta',):
    """
    对要统计的目标区域的每一个位点进行统计，统计每个位点发生测序错误的情况。
    :param bed: 要统计的目标区域
    :param bam: 比对结果文件
    :param prefix: 输出结果文件前缀
    :param center_size: 该值为统计时以当前位点为中心并向两边延申的长度，默认为1
    :param only_alt_site: 如果提供该参数，则仅仅对存在突变或测序错误的位点进行统计，和参考基因组完全一致的不参与统计。默认不提供该参数。
    :param exclude_sites: 统计时需要排除的位点列表文件，可以是vcf，要求前两列是染色体和坐标信息（1-based)，通常为已知突变位点信息
    :param genome: 参考基因组序列文件路径
    :return: 默认输出两个文件，一个文件为*.each_site.txt，记录每一个位点的信息；
        一个文件为centered{center_size}_site.json'，记录位点中心序列对应的测序错误频率信息。
    """
    gn =pysam.FastaFile(genome)
    bam = pysam.AlignmentFile(bam)
    mdict = dict()
    # excludes = {29497966,178922273,55228052,55249062,43613842,25386062,25391238,25400205}
    if exclude_sites:
        excludes = {tuple(x.strip().split()[:2]) for x in open(exclude_sites) if not x.startswith('#')}
    else:
        excludes = set()
    with open(bed) as f, open(f'{prefix}.each_site.txt', 'w') as fw:
        for line in f:
            r, s, t, _ = line.strip().split()
            s, t = int(s), int(t)
            seq_quals = get_seq_qual(r, s, t, bam)

            for pos, seq_qual in zip(range(s, t), seq_quals):
                if (r, str(pos+1)) in excludes:
                    # print('skip', pos)
                    continue
                seq_counter = Counter(x.upper() for x in seq_qual[0])
                qual_counter = Counter(seq_qual[1])
                center_seq = get_center_seq(r, pos, gn, size=center_size)
                ref = center_seq[len(center_seq)//2]

                alt_base = seq_counter.keys() - {ref.upper(), ref.lower()}
                alt_num = sum(seq_counter[x] for x in alt_base)
                total = sum(seq_counter.values())
                alt_freq = alt_num / total
                if alt_freq > 0.15 and sum(seq_counter.values()) > 10:
                    # 这里的判断可以较大程度的剔除真实突变的位点，目的是为了防止这些突变导致高估测序错误率
                    # 当然，对于测序比较浅的时候，这个过滤的作用不大
                    continue

                if not only_alt_site:
                    # 统计所有位点，即使没有发生任何突变的位点
                    alt_base = {True}
                else:
                    pass
                    # if 'C' not in alt_base:
                    #     alt_base = {}


                if alt_base:
                    mdict.setdefault(center_seq, Counter())
                    mdict[center_seq] += seq_counter
                    # alt_num = sum(seq_counter[x] for x in alt_base)
                    # mdict[center_seq] += {'alt':alt_num}
                    info = [
                        r, pos, ref,
                        center_seq,
                        len(seq_qual[0]),
                        dict(seq_counter.most_common()),
                        dict(qual_counter.most_common())
                    ]
                    fw.write('\t'.join(str(x) for x in info)+'\n')

    print(len(mdict))

    # 注释到下面的代码的原因, 因为要考虑反向互补的情况:由于观察到A->C时，可能是pcr时把A错配成C, 也有可能pcr时把互补链的T错配为G
    # mdict = dict(sorted(zip(mdict.keys(), mdict.values()), key=lambda x:x[0]))
    # with open(f'{prefix}.centered_site.txt', 'w') as fw:
    #     for k, v in mdict.items():
    #         # print(k, v)
    #         # total = sum(v.values()) - v['alt']
    #         total = sum(v.values())
    #         v = dict(v.most_common())
    #         freq = [(x, v[x]/total) for x in v]
    #         fw.write('\t'.join([k] + [f'{x}:{y:.2e}' for x, y in freq])+'\n')

    # 由于观察到A->C时，可能是pcr时把A错配成C, 也有可能pcr时把互补链的T错配为G
    # 当参考基因组为CAT时，中间的A测错为C的概率，可以由mdict中的CAT和GTA的值共同决定
    result = dict()
    for base in ['A', 'T', 'C', 'G', '']:
        for key in mdict:
            result.setdefault(key, Counter())
            result[key][base] = mdict[key][base]
            r_key = reverse_complement(key)
            r_base = reverse_complement(base)
            if r_key in mdict and r_base in mdict[r_key]:
                result[key][base] += mdict[r_key][r_base]

    mdict = dict(sorted(zip(result.keys(), result.values()), key=lambda x: x[0]))
    result = dict()
    # with open(f'{prefix}.centered{center_size}_site.txt', 'w') as fw:
    for k, v in mdict.items():
        # print(k, v)
        # total = sum(v.values()) - v['alt']
        total = sum(v.values())
        v = dict(v.most_common())
        freq = dict([(x, v[x] / total) for x in v])
        result[k] = freq
        # fw.write('\t'.join([k] + [f'{x}:{y:.2e}' for x, y in freq]) + '\n')
        # fw.write(f'{k}\t{freq}\n')
    with open(f'{prefix}.centered{center_size}_site.json', 'w') as fw:
        json.dump(result, fw, indent=4)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())


