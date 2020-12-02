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
    seq_qual = [[col.get_query_sequences(add_indels=True), col.get_query_qualities()] for col in cols]
    # A pattern
    # `\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion
    # between this reference position and the next reference
    # position. The length of the insertion is given by the
    # integer in the pattern, followed by the inserted
    # sequence. Similarly, a pattern `-[0-9]+[ACGTNacgtn]+'
    # represents a deletion from the reference. The deleted bases
    # will be presented as `*' in the following lines.
    for index in range(len(seq_qual)):
        seqs = []
        for i in seq_qual[index][0]:
            if '+' in i:
                # 当前碱基后面有插入的碱基
                seqs.append('I')
            elif '*' in i:
                # 当前碱基被删除
                seqs.append('D')
            elif '-' in i:
                # 当前碱基的后一个碱基发生删除
                seqs.append(i.split('-')[0])
            else:
                seqs.append(i)
        seq_qual[index][0] = seqs
    return seq_qual


def get_center_seq(contig, pos, genome, sizes=(1, 1)):
    return genome.fetch(contig, pos-sizes[0], pos+sizes[1]+1)


def reverse_complement(seq):
    """
    :param seq: 输入序列
    :return: 返回reversed的序列
    """
    seq = seq.upper()
    complement = dict(zip(list("ATCG"), list("TAGC")))
    return ''.join(complement[base] if base in complement else base for base in seq[::-1])


def run(bed, bam, prefix, center_size=(1, 1), exclude_sites=None,
        genome='/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta',):
    """
    对要目标区域的每一个位点进行统计，尝试得出测序错误条件下：碱基发生转换的条件概率
    分析思路：
    针对6种突变或测序错误情况'ref -> A/T/C/G/I/D' 进行6种统计，其中I:insertion D:deletion
    每一轮统计结果的意义：
    如关注的是ref突变为A，那么针对整个目标区域完成统计后，得到在已知ref突变为C的条件下，T/C/G突变为A的概率分别是多少,
    这里还可以假设前后碱基影响测序错误率，所以考虑把ref前后一个碱基纳入，如下示例，第1个概率表示ACA中的C被错误测为A的概率为0.003
    {
    "A": {
        "ACA": { "A": 0.003 },
        "AGC": { "A": 0.001 },
        "CTT": { "A": 0.0008 }，
        ...},
    ...}
    注:统计时会把测序DP>=10且AF>=10%的位点剔除，认为这些是真实突变，并不是测序错误导致
    :param bed: 要统计的目标区域
    :param bam: 比对结果文件
    :param prefix: 输出结果文件前缀
    :param center_size: 该值为统计时以当前位点为中心并向两边延申的长度，默认为1
    :param exclude_sites: 统计时需要排除的位点列表文件，可以是vcf，要求前两列是染色体和坐标信息（1-based)，通常为已知突变位点信息
    :param genome: 参考基因组序列文件路径, 用于获取参考序列
    :return: 默认输出两个文件。
        一个文件为*.each_site.txt，记录目标区域的每一个碱基的突变情况；
        一个文件为centered{center_size}_site.json'，记录位点中心序列对应的测序错误频率信息。
    """
    gn =pysam.FastaFile(genome)
    bam = pysam.AlignmentFile(bam)
    if exclude_sites:
        excludes = {tuple(x.strip().split()[:2]) for x in open(exclude_sites) if not x.startswith('#')}
    else:
        excludes = set()

    alt_raw_dict = {x:dict() for x in ['A', 'T', 'C', 'G', 'I', 'D']}
    # mdict = dict()
    with open(bed) as f, open(f'{prefix}.each_site.txt', 'w') as fw:
        header = ['contig', 'pos', 'ref', 'centered_ref', 'total_reads', 'alt_stat', 'base_qual_stat']
        fw.write('\t'.join(header)+'\n')
        for line in f:
            if line.startswith('track'):
                continue
            r, s, t, _ = line.strip().split()
            s, t = int(s), int(t)
            seq_quals = get_seq_qual(r, s, t, bam)

            for pos, seq_qual in zip(range(s, t), seq_quals):
                if (r, str(pos+1)) in excludes:
                    # print('skip', pos)
                    continue
                seq_counter = Counter(x.upper() for x in seq_qual[0])
                qual_counter = Counter(seq_qual[1])
                center_seq = get_center_seq(r, pos, gn, sizes=center_size).upper()
                ref = center_seq[len(center_seq)//2]

                alt_types = seq_counter.keys() - {ref.upper()}
                alt_num = sum(seq_counter[x] for x in alt_types)
                total = sum(seq_counter.values())
                alt_freq = alt_num / total

                if alt_freq >= 0.1 and sum(seq_counter.values()) >= 10:
                    # 这里的判断可以较大程度的剔除真实突变的位点，目的是为了防止这些突变导致高估测序错误率
                    # 当然，对于测序比较浅的时候，这个过滤的作用不大
                    continue

                for alt_type in alt_types:
                    alt_raw_dict[alt_type].setdefault(center_seq, Counter())
                    alt_raw_dict[alt_type][center_seq] += seq_counter

                if alt_types:
                    # mdict.setdefault(center_seq, Counter())
                    # mdict[center_seq] += seq_counter
                    info = [
                        r, pos, ref,
                        center_seq,
                        len(seq_qual[0]),
                        dict(seq_counter.most_common()),
                        dict(qual_counter.most_common())
                    ]
                    fw.write('\t'.join(str(x) for x in info)+'\n')

    print(dict(zip(
        alt_raw_dict.keys(), (len(v) for k, v in alt_raw_dict.items())
    )))

    alt_dict = dict()
    for alt_type in alt_raw_dict:
        # 合并：由于观察到A->C时，可能是pcr时把A错配成C, 也有可能pcr时把互补链的T错配为G
        result = dict()
        mdict = alt_raw_dict[alt_type]
        for base in ['A', 'T', 'C', 'G', 'I', 'D']:
            for key in mdict:
                result.setdefault(key, Counter())
                result[key][base] = mdict[key][base]
                r_key = reverse_complement(key)
                r_base = reverse_complement(base)
                if r_key in mdict and r_base in mdict[r_key]:
                    result[key][base] += mdict[r_key][r_base]

        # convert to freq
        freq_result = dict()
        for k in sorted(result.keys()):
            v = result[k]
            total = sum(v.values())
            v = dict(v.most_common())  # 排序
            # 仅仅输出目标突变的概率
            freq = {x: v[x]/total for x in v if x==alt_type}
            freq_result[k] = freq
        freq_result = dict(sorted(
            zip(freq_result.keys(), freq_result.values()),
            key=lambda x:sum(x[1].values()), reverse=True)
        )
        alt_dict[alt_type] = freq_result
        # output
        # with open(f'{prefix}.centered{center_size}_site.{alt_type}.json', 'w') as fw:
        #     json.dump(freq_result, fw, indent=4)

    # print(alt_dict.keys())
    with open(f'{prefix}.centered{center_size[0]}{center_size[1]}_site.json', 'w') as fw:
        json.dump(alt_dict, fw, indent=4)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())


