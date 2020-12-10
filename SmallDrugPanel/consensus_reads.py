import re
import statistics
from collections import Counter
import logging
from concurrent.futures import ThreadPoolExecutor
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


def group_reads(bam, primer, contig, start:int, end:int, method='directional'):
    if type(bam) == str:
        bam = pysam.AlignmentFile(bam)
    group_umi_by_primer = dict()
    for read in bam.fetch(contig, start, end):
        # 只基于read1和primary alignment 统计
        if (not read.is_read1) or read.is_secondary or read.is_unmapped:
            continue
        if read.reference_id == read.next_reference_id and (read.next_reference_start is not None):
            if read.next_reference_start - read.reference_end <= 0:
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
    result = dict()
    for primer, umi_read in group_umi_by_primer.items():
        umi_lst = [x[0].encode() for x in umi_read]
        group_lst = cluster(Counter(umi_lst), threshold=1)
        for group in group_lst:
            represent_umi = group[0].decode()
            for umi in group:
                for u, r, s, e in umi_read:
                    if u == umi.decode():
                        result[r] = (primer+':'+represent_umi, s, e)
    return result


def consensus_base(bases, quals, insertions, depth):
    # 返回1(全票支持), 返回2(大部分选票支持），返回3（票选相当，靠qual取胜）
    if not bases:
        # 当前位置完全没有read支持
        return '', [25], [1]
    elif len(bases) < depth * 0.33:
        # 当前位点支持的reads数量少于当前考察范围的测序深度的33%时，认为当前位点无read支持
        return '', [25], [2]
    else:
        # 接下来的所有情况必须满足：当前位点测序深度不能低于depth的33%
        pass

    # read1和read2存在overlap时，肯定会导致overlap部分存在大小写的不一致，这不利于consensus
    bases = [x.upper() for x in bases]
    base_counter = Counter(bases)
    top = base_counter.most_common(3)
    if len(top) == 1:
        represent = bases[0]
        confidence = 1
    else:
        if top[0][1]/top[1][1] > 1.5:
            represent = top[0][0]
            confidence = 2
        else:
            confidence = 3
            # 当第一名和第二名相差很接近时，根据碱基质量的平均值来选择base
            # 当只有2个read时，同样如此处理，虽然或许不大对
            # 也有可能只有3条reads，而且它们互不一致
            # 所以当支持的reads比例少于一定比例时不应该进行consensus？
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

    # 针对insertion进行consensus
    if represent.startswith('I'):
        inserts = [i for b, i in zip(bases, insertions) if b == represent]
        represent = ''
        for i in range(len(inserts)):
            i_counter = Counter(x[i] for x in inserts)
            represent += i_counter.most_common(1)[0]
    rep_len = 1 if represent == '' else len(represent)
    return represent, [rep_qual]*rep_len, [confidence]*rep_len


def consensus_reads(bam, primer, read_type=64, min_bq=0):
    # primer example:  chr1:115252197-115252227:31:+:NRAS, 坐标为1-based
    # logger = set_logger('consensus.log.txt', logger_id='consensus')
    primer_lst = primer.split(':')
    contig = primer_lst[0]
    pos = primer_lst[1]
    extend = 350 if read_type == 128 else 200
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
        Exception('Only support read_type of [64, 128]')

    if type(bam) == str:
        bam = pysam.AlignmentFile(bam)
    # group reads by primer and umi
    read2group = group_reads(bam, primer, contig, start, end, method='directional')
    group2read = dict()
    group2overlap_start = dict()
    group2overlap_end = dict()
    for k, v in read2group.items():
        group2read.setdefault(v[0], set())
        group2read[v[0]].add(k)
        group2overlap_start.setdefault(v[0], list())
        group2overlap_end.setdefault(v[0], list())
        group2overlap_start[v[0]].append(v[1])
        group2overlap_end[v[0]].append(v[2])
    else:
        print(f'there are {len(group2read)} groups!')
        group_sizes = [len(v) for k, v in group2read.items()]
        median_size = statistics.median(group_sizes)
        print('median group size', median_size)

    # 获得每个分组的pileup列表
    cols = bam.pileup(
        contig, start, end,
        stepper='samtools',
        truncate=True,
        min_base_quality=min_bq,
        ignore_orphans=False,
        ignore_overlaps=False,
        max_depth=30000,
        fastafile=None,
        flag_require=read_type,
    )
    ## 初始化pileup存储容器
    pileup_dict = dict()
    for group_name in group2read:
        pileup_dict[group_name] = dict()

    ## 逐一循环每一个位置，并按组分配
    for col in cols:
        for base, qual, read in zip(
                # 如不加add_indels参数，那么将无法知晓插入的碱基序列
                col.get_query_sequences(add_indels=True),
                col.get_query_qualities(),
                col.get_query_names()):
            # print(base, qual, read)
            if read not in read2group:
                continue
            insertion = ''
            if '+' in base:
                # such as 'G+2AT'
                # 对于insertion，也可能测错误，所以需要单独collapse，
                # 使用'I'代替是否有insertion，然后对insertion进行call consensus
                # 只能对相同长度的insertion进行校正
                insertion = re.sub('\+\d+',  '', base)
                base = 'I'+str(len(insertion))
                # logger.info('insertion:' + read.to_string())
            elif '-' in base:
                # such as 'G-1N'
                base = re.split('-\d+', base)[0]
            elif '*' in base:
                base = ''

            data = pileup_dict[read2group[read][0]]
            pos = col.reference_pos
            # data.setdefault(pos, [[base], [qual], [read], [insertion]])
            data.setdefault(pos, [[], [], [], []])
            data[pos][0].append(base)
            data[pos][1].append(qual)
            data[pos][2].append(read)
            data[pos][3].append(insertion)

    # 上面对于read中存在clip的情况无法处理，如果真实read存在clipped,
    # 基于上面的结果，最后的consenus reads将不包含clipped掉的reads
    # 下面尝试把clipped的序列找回来，使consenus read结果包含clipped碱基
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
            for pos, data in pileup_dict[read2group[read.query_name][0]].items():
                if pos in tmp_dict:
                    for i, (b, r) in enumerate(zip(data[0], data[2])):
                        if r == read.query_name:
                            data[0][i] = tmp_dict[pos]
                            if b != '': # 测试后可以删除该判断
                                Exception('clip_contradiction:'+read.to_string())
                            break

    with ThreadPoolExecutor(1) as pool:
        tasks = []
        for group_name, data in pileup_dict.items():
            tasks.append([group_name, pool.submit(consensus_read, data)])

        for group_name, task in tasks:
            consistent_bases, median_cov, top = task.result()
            print(f'>{group_name}')
            print(f'median coverage is {median_cov}')
            print(f'top3 coverage frequency is {top}')
            print(''.join(x[0] for x in consistent_bases))
            print(consistent_bases)
            print(''.join(''.join(str(i) for i in x[2]) for x in consistent_bases if x[0] != '' or x[0] != 'X'))


def consensus_read(data):
    consistent_bases = []
    coverages = [len(v[0]) for k, v in data.items()]
    top = Counter(coverages).most_common(3)
    for pos, (bases, quals, reads, insertions) in data.items():
        coverages.append(len(bases))
        consistent_bases.append(consensus_base(bases, quals, insertions, top[0][0]))
    return consistent_bases, statistics.median(coverages), top


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
