import os
import re
import statistics
from collections import Counter
from functools import partial
import logging
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Pool, Manager
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pysam
from umi_tools import UMIClusterer

"""
了解pysam中pileup表示突变的表达形式:
# A pattern
# `\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion
# between this reference position and the next reference
# position. The length of the insertion is given by the
# integer in the pattern, followed by the inserted
# sequence. Similarly, a pattern `-[0-9]+[ACGTNacgtn]+'
# represents a deletion from the reference. The deleted bases
# will be presented as `*' in the following lines.
# 对于indel，返回的qual并不是空的，可能是pysam根据前后碱基的质量直接填充过来的

"""


def set_logger(name='log.info', logger_id='x'):
    logger = logging.getLogger(logger_id)
    fh = logging.FileHandler(name, mode='w+')
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    return logger


def group_reads(bam, primer, contig, start: int, end: int, method='directional', out_bam=None):
    if type(bam) == str:
        bam = pysam.AlignmentFile(bam)
    group_umi_by_primer = dict()
    for read in bam.fetch(contig, start, end):
        # 只基于read1和primary alignment 统计
        if (not read.is_read1) or read.is_secondary or read.is_unmapped:
            continue
        if read.reference_id == read.next_reference_id and (read.next_reference_start is not None):
            if read.next_reference_start - read.reference_end <= 0:
                # overlap is true
                overlap_start = read.next_reference_start
                overlap_end = read.reference_end
            else:
                overlap_start = overlap_end = read.reference_end + 1
        else:
            overlap_start = overlap_end = read.reference_end + 1
        read_name = read.query_name
        primer_lst = read_name.split(':', 5)[:5]
        if ':'.join(primer_lst) == primer:
            umi = read_name.rsplit(':', 1)[1]
            group_umi_by_primer.setdefault(primer, list())
            group_umi_by_primer[primer].append((umi, read_name, overlap_start, overlap_end))

    cluster = UMIClusterer(cluster_method=method)
    # cluster = UMIClusterer(cluster_method='adjacency')
    result = dict()
    for primer, umi_read in group_umi_by_primer.items():
        group_lst = cluster(Counter(x[0].encode() for x in umi_read), threshold=2)
        for group in group_lst:
            represent_umi = group[0].decode()
            umi_set = set(group)
            for u, r, s, e in umi_read:
                if u.encode() in umi_set:
                    result[r] = (primer+':'+represent_umi, s, e)

    # 通读输出新的bam, 修改read_group
    if out_bam is not None:
        ob = pysam.AlignmentFile(out_bam, 'wb', template=bam)
        for read in bam.fetch(contig, start, end):
            if read.query_name in result:
                read.set_tag('RG', result[read.query_name][0])
                ob.write(read)
        ob.close()
        pysam.index(out_bam)
    return result


def consensus_base(bases, quals, insertions, depth, contig, position, ref_seq):
    # 返回3(相对depth以75%以上的比例票支持), 返回2(大部分选票支持），返回1（票选相当，靠qual取胜）
    # 输入的base中，'D'表示deletion, ‘I\d’表示insertion ‘S[ATCG]'表示clipped
    current_depth = len(bases)
    depth = current_depth if current_depth > depth else depth

    if not bases:
        # 当前位置完全没有read支持
        return 'X', [25], [3], current_depth, contig, position, ref_seq
    elif len(bases) < depth * 0.35:
        # 当前位点支持的reads数量少于当前考察范围的测序深度的35%时，认为当前位点无read支持
        return 'X', [25], [2], depth-current_depth, contig, position, ref_seq
    else:
        # 接下来的所有情况必须满足：当前位点测序深度不能低于depth的35%
        pass

    bases = [x.upper() for x in bases]
    base_counter = Counter(bases)
    top = base_counter.most_common(3)
    top1_ratio = top[0][1]/depth
    if len(top) == 1:
        represent = top[0][0]
        if depth <= 2:
            confidence = 1
        else:
            confidence = 3
    elif top1_ratio >= 0.75 and depth > 1:
        represent = top[0][0]
        confidence = 3
    else:
        if top[0][1]/top[1][1] > 1.5:
            # such as 5/3 > 1.5
            represent = top[0][0]
            confidence = 2
        else:
            # such as 3/2 < 1.5
            confidence = 1
            # 当第一名和第二名相差很接近时，根据碱基质量的平均值来选择base
            # 当只有2个read时，同样如此处理，虽然或许不大对
            # 也有可能只有3条reads，而且它们互不一致
            bqd = dict()
            for b, q in zip(bases, quals):
                bqd.setdefault(b, set())
                bqd[b].add(q)
            top1_mean_qual = statistics.mean(bqd[top[0][0]])
            top2_mean_qual = statistics.mean(bqd[top[1][0]])
            if top1_mean_qual >= top2_mean_qual:
                represent = top[0][0]
            else:
                represent = top[1][0]

    # 取最大支持的qual作为consensus结果的qual
    rep_qual = max(q for b, q in zip(bases, quals) if b == represent)
    support_depth = sum(b == represent for b in bases)
    # 针对insertion进行consensus
    if represent.startswith('I'):
        inserts = [i for b, i in zip(bases, insertions) if b == represent]
        represent = ''
        for insert_lst in zip(*inserts):
            i_counter = Counter(insert_lst)
            # 取出现次数最多的作为代表
            represent += i_counter.most_common(1)[0][0]
        # 用<>表示整体的插入序列, 用()表示deletion
        represent = '<' + represent + '>'
    rep_len = len(represent)
    return represent, [rep_qual]*rep_len, [confidence]*rep_len, support_depth, contig, position, ref_seq


def consensus_read(data):
    coverages = [len(set(v[2])) for k, v in data.items()]
    median_depth = statistics.median_high(coverages)
    consistent_bases = []
    for (pos, ref, c), (bases, quals, reads, insertions, overlap) in data.items():
        consistent_bases.append(consensus_base(bases, quals, insertions, median_depth, c, pos, ref))
    return consistent_bases, median_depth, Counter(coverages).most_common(3)


def format_consensus_bases(base_info):
    # base, qual, confidence, support_depth, contig, position, ref_seq
    b, q, c, s, _, p, r = base_info
    if b.startswith('<'):
        return b[1:-1], ''.join([chr(q[0]+33)] * (len(b) - 2)), ''.join([str(c[0])] * (len(b) - 2)), p
    elif b.startswith(('S', '(')):
        return b[1], chr(q[1]+33), str(c[1]), p
    elif b == '-':
        return '', '', '', p
    else:
        return b[0], chr(q[0]+33), str(c[0]), p


def consensus_reads(bam, primer, read_type=64, min_bq=0, fq_lst=None, ignore_overlaps=False,
                    genome='/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta'):
    # primer example:  chr1:115252197-115252227:31:+:NRAS, 坐标为1-based
    # logger = set_logger('consensus.log.txt', logger_id='consensus')
    primer_lst = primer.split(':')
    contig = primer_lst[0]
    pos = primer_lst[1]
    if read_type == 64:
        # 只解析read1
        extend = 200
    else:
        # 不仅仅解析read1, 还可能解析read2
        extend = 400
    if primer_lst[-2] == "+":
        start = int(pos.split('-')[1])
        end = start + extend
    else:
        end = int(pos.split('-')[0])
        start = end - extend
    print(f"parsing region of {contig}:{start}-{end}")

    if read_type == 64:
        print('getting consensus read1')
    elif read_type == 128:
        print('getting consensus read2')
    else:
        # Exception('Only support read_type of [64, 128]')
        print('getting consensus for both read1 and read2')

    if type(bam) == str:
        bam = pysam.AlignmentFile(bam)

    # group reads by primer and umi
    read2group = group_reads(bam, primer, contig, start, end, method='directional')
    if not read2group:
        print(f'primer {primer} has no corresponding reads!')
    group2read = dict()
    group2overlap = dict()
    for k, v in read2group.items():
        group2read.setdefault(v[0], set())
        group2read[v[0]].add(k)
        group2overlap.setdefault(v[0], [])
        if read_type == 128 or read_type == 64:
            # 单端read，直接判断为非overlap
            group2overlap[v[0]].append(0)
        else:
            # 判断是否overlap
            if v[1] == v[2]:
                group2overlap[v[0]].append(0)
            else:
                group2overlap[v[0]].append(1)
    else:
        print(f'there are {len(group2read)} groups for primer {primer}')
        group_size_dict = {k:len(v) for k,v in group2read.items()}
        group_sizes = group_size_dict.values()
        median_size = statistics.median_high(group_sizes)
        print('(min, median, max) group size', (min(group_sizes), median_size, max(group_sizes)))
        for k, v in group2overlap.items():
            if sum(v)/len(v) > 0.6:
                # 超过60%的read说当前read出现overlap
                group2overlap[k] = 1
            else:
                group2overlap[k] = 0

    # 获得每个分组的pileup列表, 返回的第一个位点不一定是start指定的，而是由实际比对情况决定的！
    # 对于clipped 的base， pileup不会纳入，因此后面有代码把这些base找回来。
    cols = bam.pileup(
        contig, start, end,
        stepper='samtools',
        truncate=True,
        min_base_quality=min_bq,
        ignore_orphans=False,
        ignore_overlaps=ignore_overlaps,  # set 为True则意味着取质量高的base作为代表
        # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
        max_depth=300000,
        fastafile=None,
        flag_require=read_type,
        flag_filter=4,
    )
    # 在pos处带入ref信息
    genome = pysam.FastaFile(genome)

    ## 初始化pileup存储容器
    final_used_read = dict()
    pileup_dict = dict()
    for group_name in group2read:
        pileup_dict[group_name] = dict()

    ## 逐一循环每一个位置，并按组分配
    for col in cols:
        ref = genome.fetch(contig, col.reference_pos, col.reference_pos+1).upper()
        for base, qual, read in zip(
                # 如不加add_indels参数，那么将无法知晓插入的碱基序列
                col.get_query_sequences(add_indels=True),
                col.get_query_qualities(),
                col.get_query_names()):
            # print(base, qual, read)
            # print(col.reference_pos)
            # print('FIND', read)
            if read not in read2group:
                continue
            insertion = ''
            if '+' in base:
                # such as 'G+2AT'
                # 对于insertion，也可能测错误，所以需要单独collapse，
                # 使用'I'代替是否有insertion，然后对insertion进行call consensus
                # 只能对相同长度的insertion进行校正, 下面的insertion包含ref
                insertion = re.sub('\+\d+',  '', base).upper()
                base = 'I'+str(len(insertion))
                # logger.info('insertion:' + read.to_string())
            elif '-' in base:
                # such as 'G-1N'
                # base = re.split('-\d+', base)[0]
                deletion_len = int(base.split('-')[1].upper().rstrip('N'))
                # 用'<'表示deletion
                base = '(' + genome.fetch(contig, col.reference_pos, col.reference_pos+deletion_len+1) + ')'
            elif '*' in base:
                # 这里deletion信息已经在deletion的上一个碱基进行了合并表示
                base = '-'
            elif base == '':
                # 如果该位点没有覆盖，该处的base为''
                base = 'X'
            group_name = read2group[read][0]
            final_used_read[read] = group_name
            data = pileup_dict[group_name]
            pos = (col.reference_pos, ref, contig)
            # data.setdefault(pos, [[base], [qual], [read], [insertion], [is_overlap])
            base_info = data.setdefault(pos, [[], [], [], [], []])
            base_info[0].append(base)
            base_info[1].append(qual)
            base_info[2].append(read)
            base_info[3].append(insertion)
            base_info[4].append(group2overlap[group_name])

     # final used read dict
    final_group_size_dict = dict()
    for k, v in final_used_read.items():
        final_group_size_dict.setdefault(v, 0)
        final_group_size_dict[v] += 1

    # 当同时进行read1和read2 consensus 且 read1和read2中间没有overlap，就会出现中间有些位置没有信息
    # pileup并不记录clipped情形，上面对于read中存在clip的情况无法处理，如果真实read存在clipped,
    # 基于上面的结果，最后的consenus reads将不包含clipped掉的reads
    # 下面尝试把clipped的序列找回来，使consenus read结果包含clipped碱基
    # 思路是根据cigar判断是否clipped，然后推断发生clipped对应的参考基因组位置，并标记S
    for read in bam.fetch(contig, start, end):
        if read.query_name in read2group and read.cigartuples is not None:
            if read.cigartuples[0][0] == 4:
                # left clipped
                aln_start = read.query_alignment_start
                ref_start = read.reference_start
                clip_pos = [ref_start-i-1 for i in range(aln_start)]
                clip_base = [x for x in read.query_sequence[:aln_start][::-1]]
            elif read.cigartuples[-1][0] == 4:
                # right clipped
                aln_end = read.query_alignment_end
                clip_num = read.query_length - 1 - aln_end
                ref_end = read.reference_end
                clip_pos = [ref_end+i+1 for i in range(clip_num)]
                clip_base = [x for x in read.query_sequence[aln_end+1:]]
                # logger.info('clipped:'+read.to_string())
            else:
                clip_pos = []
                clip_base = []
            tmp_dict = dict(zip(clip_pos, clip_base))
            for (pos, ref, c), data in pileup_dict[read2group[read.query_name][0]].items():
                if pos in tmp_dict:
                    for i, (b, r) in enumerate(zip(data[0], data[2])):
                        if r == read.query_name:
                            # 找某条read对应的位置发生clipped, 每个被clipped的base带入'S'标签，也许日后可以用到
                            data[0][i] = 'S'+tmp_dict[pos]
                            break
    # consensus
    result = dict()
    for group_name, data in pileup_dict.items():
        if len(data) == 0:
            # 某个分组在目标区域没有对应的reads?, 之前发现是pileup时参数max_depth不够大导致漏掉
            print(f'{group_name} is empty, but raw group_size is {len(group2read[group_name])}?')
            if len(group2read[group_name]) > 5:
                print(group2read[group_name])
            continue
        consistent_bases, median_cov, top = consensus_read(data)
        # print(f'>{group_name} overlap:{group2overlap[group_name]}')
        # print(f'median coverage is {median_cov} and base number for Top3 coverages is {top}')

        # 制作突变需要的字典
        # umi = group_name.rsplit(':', 1)[-1]
        complex = []  # 一条reads中连续的snv合并为complex,允许中间存在一个间隔
        gap = 0
        for base, qual, confidence, alt_depth, chr_name, pos, ref in consistent_bases:
            if base[0] not in ['S', 'X', '-']:
                # clipped 和没有read支持的位置不能call突变
                if base != ref:
                    if base.startswith('<') or base.startswith('('):
                        # insertion | deletion
                        key = (chr_name, pos, ref)
                        result.setdefault(key, [])
                        result[key].append((base, confidence, alt_depth, group_name))
                    else:
                        # SNP or Complex, 该位置的信息将在后续complex分析完成后存储
                        gap = 0  # gap 重新从0开始计数
                        complex.append([pos, base, ref, confidence, alt_depth])
                else:
                    # 对于没有突变的位点，正常更新每个位点信息
                    key = (chr_name, pos, ref)
                    result.setdefault(key, [])
                    result[key].append((base, confidence, alt_depth, group_name))

                    if complex:
                        gap += 1
                        if gap < 2:
                            # 直接进入下一轮循环, 目的是控制complex突变中允许存在一个未变异的gap
                            continue

                        # 复合突变延长过程时第二次碰到没有突变位点时，把complex容器中的突变取出并合并为一个复合突变
                        if len(complex) == 1:
                            # SNP
                            pos, base, ref, confidence, alt_depth = complex[0]
                            key = (chr_name, pos, ref)
                            result.setdefault(key, [])
                            result[key].append((base, confidence, alt_depth, group_name))
                        else:
                            # complex
                            p, b, r, c, a = list(zip(*complex))
                            # 以第一个突变位置为索引进行突变信息存储
                            key = (chr_name, p[0], r[0])
                            alt = ''.join(r) + ':'+ ''.join(b)  # 带入ref信息方便variant的输出
                            result.setdefault(key, [])
                            result[key].append((alt, max(c), max(a), group_name))
                            # 完成complex分析, 逐一更新之前没有更新的位点，并且标记为非突变，因为该突变信息已经在第一个位置存储
                            for p, b, r, c, a in complex[1:]:
                                key = (chr_name, p, r)
                                result.setdefault(key, [])
                                result[key].append((r, c, a, group_name))

                        # 清空之前得到的complex, 方便下一个complex突变的存储
                        complex = []
        else:
            # 循环结束, 清理最后一个可能的complex
            if complex:
                if len(complex) == 1:
                    # SNP
                    pos, base, ref, confidence, alt_depth = complex[0]
                    key = (chr_name, pos, ref)
                    result.setdefault(key, [])
                    result[key].append((base, confidence, alt_depth, group_name))
                else:
                    # complex
                    p, b, r, c, a = list(zip(*complex))
                    # 以第一个突变位置为索引进行突变信息存储
                    key = (chr_name, p[0], r[0])
                    alt = ''.join(r) + ':' + ''.join(b)  # 带入ref信息方便variant的输出
                    result.setdefault(key, [])
                    result[key].append((alt, max(c), max(a), group_name))
                    # 完成complex分析, 逐一更新之前没有更新的位点，并且标记为非突变，因为该突变信息已经在第一个位置存储
                    for p, b, r, c, a in complex[1:]:
                        key = (chr_name, p, r)
                        result.setdefault(key, [])
                        result[key].append((r, c, a, group_name))
                complex = []

        # 为输出consensus reads做准备, fq_lst是一个可以在多进程间共享的的list
        if fq_lst is not None:
            # 想办法去掉两端的X
            consensus_seq = ''.join(x[0] for x in consistent_bases)
            lm = re.match('X+', consensus_seq)
            left_start = 0 if lm is None else len(lm.group())
            rm = re.match('X+', consensus_seq[::-1])
            right_end = len(consistent_bases) if rm is None else len(consistent_bases) - len(rm.group())

            # 整理fastq的信息
            base_info = consistent_bases[left_start:right_end]
            # 处理indel和clipped的标签信息
            base_info = map(format_consensus_bases, base_info)
            # 形成fastq信息
            base_info_dict = {p:(b, q, c) for b, q, c, p in base_info}
            min_pos = min(base_info_dict.keys())
            max_pos = max(base_info_dict.keys())
            continuous_pos = range(min_pos, max_pos+1)
            # 把没有read覆盖的位置用X填补起来
            seqs = ''.join(base_info_dict[p][0] if p in base_info_dict else 'X' for p in continuous_pos)
            # print(seqs)
            quals = ''.join(base_info_dict[p][1] if p in base_info_dict else 'X' for p in continuous_pos)
            confs = '+'  # 原本打算存储confidences信息，但感觉用途不大，放弃
            header = f'@{group_name} read_number:N:{int(median_cov)}:{group2overlap[group_name]}'
            fq_lst.append([read_type, header, seqs, confs, quals])
    else:
        if fq_lst is not None:
            fq_lst.append([primer])

    return result, final_group_size_dict


def create_vcf(vcf_path, genome='hg19', chrom_name_is_numeric=False):
    vcf = pysam.VariantFile(vcf_path, 'w')
    contig_info = [
        '##fileformat=VCFv4.2',
        f'##assembly={genome}',
        "##contig=<ID=chr1,length=249250621>",
        "##contig=<ID=chr2,length=243199373>",
        "##contig=<ID=chr3,length=198022430>",
        "##contig=<ID=chr4,length=191154276>",
        "##contig=<ID=chr5,length=180915260>",
        "##contig=<ID=chr6,length=171115067>",
        "##contig=<ID=chr7,length=159138663>",
        "##contig=<ID=chr8,length=146364022>",
        "##contig=<ID=chr9,length=141213431>",
        "##contig=<ID=chr10,length=135534747>",
        "##contig=<ID=chr11,length=135006516>",
        "##contig=<ID=chr12,length=133851895>",
        "##contig=<ID=chr13,length=115169878>",
        "##contig=<ID=chr14,length=107349540>",
        "##contig=<ID=chr15,length=102531392>",
        "##contig=<ID=chr16,length=90354753>",
        "##contig=<ID=chr17,length=81195210>",
        "##contig=<ID=chr18,length=78077248>",
        "##contig=<ID=chr19,length=59128983>",
        "##contig=<ID=chr20,length=63025520>",
        "##contig=<ID=chr21,length=48129895>",
        "##contig=<ID=chr22,length=51304566>",
        "##contig=<ID=chrMT,length=16569>",
        "##contig=<ID=chrX,length=155270560>",
        "##contig=<ID=chrY,length=59373566>",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##INFO=<ID=Confidences,Number=1,Type=String,Description="Confidence level in range of [1, 2, 3], higher value indicates more reliable during consensusing process">',
        '##INFO=<ID=RawAlt,Number=1,Type=String,Description="Number of Raw reads that support the mutation in each umi cluster">',
        '##INFO=<ID=GroupName,Number=1,Type=String,Description="group_name list">',
        '##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant Type: SNV Insertion Deletion Complex">',
        '##INFO=<ID=RawAltSum,Number=1,Type=Integer,Description="Total number of raw reads that support the mutation">',
        '##INFO=<ID=RawAltMedian,Number=1,Type=Integer,Description="Median number of raw reads that support the mutation">',
        '##INFO=<ID=ConfidenceSum,Number=1,Type=Integer,Description="Sum of confidence levels ">',
        '##INFO=<ID=ConfidenceMedian,Number=1,Type=Integer,Description="Median confidence level">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        # '##INFO=<ID=END,Number=1,Type=Integer,Description="Chr End Position">',
        '##INFO=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype: 0/1 if AF < 0.9 else 1/1 ">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">',
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
    ]
    if chrom_name_is_numeric:
        contig_info = [x.replace('=chr', '=') for x in contig_info]
    for line in contig_info:
        vcf.header.add_line(line)
    return vcf


def call_variant(result, out='mutation.vcf', min_umi_depth=5, min_alt_num=2, min_conf=4, min_raw_alt_num=5):
    # call variant
    ordered_keys = sorted(result.keys(), key=lambda x: (x[0], x[1]))
    sample = out.split('.', 1)[0]
    vcf = create_vcf(out)
    vcf.header.add_sample(sample)
    variant_number = 0
    # 目前还没有实现对相邻突变进行合并的功能
    for key in ordered_keys:
        base_info = result[key]
        contig, position, ref_seq = key
        bases = [x[0] for x in base_info]
        depth = len(bases)
        base_counter = Counter(bases)
        for base, freq in base_counter.items():
            if base.upper() != ref_seq:
                # print(base_counter)
                # ad = (depth, freq)
                af = freq / depth
                confidences = [x[1][0] for x in base_info if x[0] == base]
                covs = [x[2] for x in base_info if x[0] == base]
                group_name = [x[3] for x in base_info if x[0] == base]
                # filtering
                umi_alt_depth = len(covs)
                confidence_sum = sum(confidences)
                raw_alt_depth = sum(covs)
                if umi_alt_depth >= min_alt_num and confidence_sum >= min_conf \
                        and raw_alt_depth >= min_raw_alt_num and depth >= min_umi_depth:
                    if raw_alt_depth/umi_alt_depth < 1.2 and af < 0.02:
                        # raw support read 太小 = UMI作用失效 = 测序错误率还是很大 = 低频突变不可靠
                        continue
                    # min_conf意味着至少要有2个一致性比较高的base支持
                    variant_number += 1
                    mut_type = 'SNP'
                    # format output
                    if base.startswith('('):
                        # deletion
                        ref = base[1:-1]
                        alt = ref_seq
                        mut_type = 'Deletion'
                    elif base.startswith('<'):
                        # insertion
                        ref = ref_seq
                        alt = base[1:-1]
                        mut_type = 'Insertion'
                    elif ':' in base:
                        ref, alt = base.split(':')
                        mut_type = 'Complex'
                    else:
                        ref = ref_seq
                        alt = base

                    info = dict(
                        TYPE=mut_type,
                        AF=af,
                        DP=depth,
                        VD=umi_alt_depth,
                        ConfidenceSum=confidence_sum,
                        RawAltSum=raw_alt_depth,
                        ConfidenceMedian=statistics.median_high(confidences),
                        RawAltMedian=statistics.median_high(covs),
                        Confidences=str(confidences),
                        RawAlt=str(covs),
                        # GroupName=str(group_name),
                    )
                    record = vcf.new_record()
                    record.contig = contig
                    record.start = position
                    record.ref = ref
                    record.alts = [alt]
                    record.qual = None
                    record.filter.add('PASS')
                    record.info.update(info)
                    record.samples[sample]['GT'] = '0/1' if af < 0.9 else '1/1'
                    record.samples[sample]['DP'] = depth
                    record.samples[sample]['VD'] = freq
                    record.samples[sample]['AF'] = af
                    vcf.write(record)
    else:
        vcf.close()
        print('variant number', variant_number)


def write_fastq(fq_lst, primer_number, r1='Consensus.R1.fq', r2='Consensus.R2.fq'):
    # {primer: {group_name:[[read_name,read_seq,+, read_qual], ...]}, ...}
    primers = set()
    read1_obj = set_logger(r1, 'read1')
    read2_obj = set_logger(r2, 'read2')
    while 1:
        if fq_lst:
            info = fq_lst.pop(0)
            if len(info) == 1:
                primers.add(info[0])
                if len(primers) == primer_number:
                    break
            else:
                t, *fq_info = info
                if t in [64, 128]:
                    if 'X' not in fq_info[1]:
                        if t == 64:
                            fq_info[0] = fq_info[0].replace('read_number', '1')
                            _ = [read1_obj.info(x) for x in fq_info]
                        else:
                            fq_info[0] = fq_info[0].replace('read_number', '2')
                            _ = [read2_obj.info(x) for x in fq_info]
                    else:
                        print('Discard', fq_info)
                else:
                    # fq_info里同时包含了read1和read2
                    # 认为X为read1和read2的分割标志
                    tmp = re.split('X+', fq_info[1])
                    if len(tmp) == 1:
                        # overlap, 则read1和read2一样
                        read1_info = fq_info.copy()
                        read1_info[0] = fq_info[0].replace('read_number', '1')
                        _ = [read1_obj.info(x) for x in read1_info]

                        read2_info = fq_info.copy()
                        read2_info[0] = fq_info[0].replace('read_number', '2')
                        _ = [read2_obj.info(x) for x in read2_info]
                    elif len(tmp) == 2:
                        left_ind = fq_info[1].find('X')
                        right_ind = fq_info[1].rfind('X')
                        read1_info = [fq_info[0].replace('read_number', '1'),
                                      fq_info[1][:left_ind], fq_info[2], fq_info[3][:left_ind]]
                        _ = [read1_obj.info(x) for x in read1_info]

                        read2_info = [fq_info[0].replace('read_number', '2'),
                                      fq_info[1][right_ind+1:], fq_info[2], fq_info[3][right_ind+1:]]
                        _ = [read2_obj.info(x) for x in read2_info]
                    else:
                        print('Discard', fq_info)


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
    ax.axhline(y=mean_number, c="k", ls="--", lw=0.6, label=f'Mean={int(mean_number)}')
    # ax.axhline(y=median_number, c="r", ls="--", lw=0.6, label=f'Median={int(median_number)}')
    ax.axhline(y=mean_number * 0.25, c="r", ls="--", lw=0.5, label='25%Mean')
    # ax.axhline(y=median_number * 0.25, c="r", ls="--", lw=0.5, label='25%Median')
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


def draw_primer_umi_bar(primer_umi_dict, out_prefix):
    info_dict = dict()
    for k, v in primer_umi_dict.items():
        primer, umi = k.rsplit(':', 1)
        info_dict.setdefault(primer, 0)
        info_dict[primer] += 1
    get_gene = lambda x: x.split(':')[-1]
    genes = sorted({get_gene(x) for x in info_dict.keys()})
    colors = get_color_pool(len(genes))
    color_dict = dict(zip(genes, colors))
    colors = [color_dict[get_gene(x)] for x in info_dict.keys()]
    title = 'Consensus UMI distribution for each primer'
    out = f'{out_prefix}.' + title.title().replace(' ', '') + '.pdf'
    plot_bar(info_dict, colors=colors, fontsize=3, rotation=90, out=out, title=title)

    # plot group size distribution
    group_sizes = primer_umi_dict.values()
    plt.hist(group_sizes, bins=50)
    plt.gca().set(title='Frequency Histogram', ylabel='Frequency', xlabel='Group size')
    plt.savefig(f'{out_prefix}.'+'UMIGroupSize.pdf')
    plt.close()


def run_all(primers, bam, read_type=0, cores=8, out_prefix='result',  min_bq=10,
            min_umi_depth=8, min_alt_num=2, min_conf=5, min_raw_alt_num=5,
            genome='/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta'):
    if not os.path.exists(bam):
        raise Exception(f'{bam} not found!')
    primers = [x.strip()[1:] for x in open(primers) if x.startswith('>')]
    cores = len(primers) if len(primers) <= cores else cores
    # 开拓进程之间的共享空间, 即使用一个进程间可以共享的list，
    # N-1个进程往list里添加信息，剩下一个进程从list清空信息并输出到同一个文件

    result = dict()
    primer_umi_group = dict()
    out_fq1 = f'{out_prefix}.R1.fq'
    out_fq2 = f'{out_prefix}.R2.fq'
    if cores <= 1:
        fq_lst = []
        for primer in primers:
            r, g = consensus_reads(bam, primer, read_type, min_bq, fq_lst, genome=genome)
            result.update(r)
            primer_umi_group.update(g)
            write_fastq(fq_lst, 1, r1=out_fq1, r2=out_fq2)
    else:
        manager = Manager()
        fq_lst = manager.list()
        get_consensus_reads = partial(
            consensus_reads, bam, read_type=read_type, min_bq=min_bq, fq_lst=fq_lst, genome=genome
        )
        with ProcessPoolExecutor(cores) as pool:
            tasks = [pool.submit(write_fastq, fq_lst, len(primers), out_fq1, out_fq2)]
            for each in primers:
                tasks.append(pool.submit(get_consensus_reads, each))
            futures = [x.result() for x in tasks[1:]]
            tasks[0].result()
            for r, g in futures:
                result.update(r)
                primer_umi_group.update(g)

    draw_primer_umi_bar(primer_umi_group, out_prefix=out_prefix)
    with open(f'{out_prefix}.primer_umi_stat.txt', 'w') as f:
        total_reads = sum(primer_umi_group.values())
        total_groups = len(primer_umi_group)
        median_group_size = statistics.median_high(primer_umi_group.values())
        group_size_over_5 = sum(x>=5 for x in primer_umi_group.values())
        group_size_over_10 = sum(x>=10 for x in primer_umi_group.values())
        group_size_over_20 = sum(x>=20 for x in primer_umi_group.values())
        group_size_over_50 = sum(x>=50 for x in primer_umi_group.values())
        f.write(f'##total_used_paired_reads: {total_reads}\n')
        f.write(f'##group_number(grouped by both Primer and UMI): {total_groups}\n')
        f.write(f'##median_group_size: {median_group_size}\n')
        f.write(f'##number_of_group_with_size_over_5: {group_size_over_5}\n')
        f.write(f'##number_of_group_with_size_over_10: {group_size_over_10}\n')
        f.write(f'##number_of_group_with_size_over_20: {group_size_over_20}\n')
        f.write(f'##number_of_group_with_size_over_50: {group_size_over_50}\n')
        f.write(f'#group_name\tgroup_size\n')
        for k, v in sorted(primer_umi_group.items(), key=lambda x:(x[0], -x[1])):
            f.write(f'{k}\t{v}\n')

    call_variant(result, f'{out_prefix}.mutation.vcf', min_umi_depth=min_umi_depth, min_alt_num=min_alt_num,
                 min_conf=min_conf, min_raw_alt_num=min_raw_alt_num)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
