import re
import statistics
from collections import Counter
from functools import partial
import logging
from concurrent.futures import ProcessPoolExecutor
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
    result = dict()
    for primer, umi_read in group_umi_by_primer.items():
        group_lst = cluster(Counter(x[0].encode() for x in umi_read), threshold=1)
        for group in group_lst:
            represent_umi = group[0].decode()
            umi_set = set(group)
            for u, r, s, e in umi_read:
                if u.encode() in umi_set:
                    result[r] = (primer+':'+represent_umi, s, e)
    return result


def consensus_base(bases, quals, insertions, depth, contig, position, ref_seq):
    # 返回3(相对depth以75%以上的比例票支持), 返回2(大部分选票支持），返回1（票选相当，靠qual取胜）
    # 输入的base中，'D'表示deletion, ‘I\d’表示insertion ‘S[ATCG]'表示clipped
    current_depth = len(bases)

    if not bases:
        # 当前位置完全没有read支持
        return 'X', [25], [3], current_depth, contig, position, ref_seq
    elif len(bases) < depth * 0.35:
        # 当前位点支持的reads数量少于当前考察范围的测序深度的35%时，认为当前位点无read支持
        # 这里有可能原本表示deletion的碱基表示成X，所以后续得到的序列如果中间出现X，可以考虑把X替换为D
        # 更有可能把clipped的位置表示未X，这是这段程序的目的之一
        return 'X', [25], [2], depth-current_depth, contig, position, ref_seq
    else:
        # 接下来的所有情况必须满足：当前位点测序深度不能低于depth的35%
        pass

    bases = [x.upper() for x in bases]
    base_counter = Counter(bases)
    top = base_counter.most_common(3)
    top1_ratio = top[0][1]/depth
    if top1_ratio >= 0.75:
        represent = top[0][0]
        confidence = 3
    elif len(top) == 1:
        # 只有一种碱基
        represent = top[0][0]
        confidence = 2
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
    support_depth = sum(b == represent for b in bases)
    # 针对insertion进行consensus
    if represent.startswith('I'):
        inserts = [i for b, i in zip(bases, insertions) if b == represent]
        represent = ''
        for insert_lst in zip(*inserts):
            i_counter = Counter(insert_lst)
            # 取出现次数最多的作为代表
            represent += i_counter.most_common(1)[0][0]
        # 用<>表示整体的插入序列
        represent = '<' + represent + '>'
    rep_len = len(represent)
    return represent, [rep_qual]*rep_len, [confidence]*rep_len, support_depth, contig, position, ref_seq


def consensus_reads(bam, primer, read_type=64, min_bq=0, genome='/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta'):
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
        extend = 350
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
        print('primer has no corresponding reads!')
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
        print(f'there are {len(group2read)} groups!')
        group_sizes = [len(v) for k, v in group2read.items()]
        median_size = statistics.median(group_sizes)
        print('median group size', median_size)
        for k, v in group2overlap.items():
            if sum(v)/len(v) > 0.6:
                # 超过60%的read说当前read出现overlap
                group2overlap[k] = 1
            else:
                group2overlap[k] = 0

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
        flag_filter=4,
    )
    # 在pos处带入ref信息
    genome = pysam.FastaFile(genome)
    refs = genome.fetch(contig, start, end)
    ## 初始化pileup存储容器
    pileup_dict = dict()
    for group_name in group2read:
        pileup_dict[group_name] = dict()

    ## 逐一循环每一个位置，并按组分配
    for col, ref in zip(cols, refs):
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
                # 这里使得每一个碱基的deletion都用D表示，所以长的deletion只能分开识别后最后合并
                base = 'D'
            group_name = read2group[read][0]
            data = pileup_dict[group_name]
            pos = (col.reference_pos, ref, contig)
            # data.setdefault(pos, [[base], [qual], [read], [insertion], [is_overlap])
            base_info = data.setdefault(pos, [[], [], [], [], []])
            base_info[0].append(base)
            base_info[1].append(qual)
            base_info[2].append(read)
            base_info[3].append(insertion)
            base_info[4].append(group2overlap[group_name])
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
            # 某个分组在目标区域没有对应的reads,
            # 当read_type=64或128时，由于fetch也会同时抓出read2和read1，下面这种情况经常发生
            # 当read_type=0时，read1和read2都要考虑的，那如果还出现这种情况？
            print(f'{group_name} is empty, but raw group_size is {len(group2read[group_name])}, ???')
            if len(group2read[group_name]) > 20:
                print(group2read[group_name])
            continue
        consistent_bases, median_cov, top = consensus_read(data)
        print(f'>{group_name}')
        print(f'median coverage is {median_cov}')
        print(f'top3 coverage frequency is {top}')
        consensus_seq = ''.join(x[0] for x in consistent_bases).strip('X')
        print(consensus_seq)
        print(''.join(''.join(str(i) for i in x[2]) for x in consistent_bases if x[0] != '' and x[0] != 'X'))
        # print(consistent_bases)
        for base, qual, confidence, depth, *key in consistent_bases:
            # base = base.replace('S', '').replace('<', '').replace('>', '')
            base = base.replace('S', '')
            if base != 'X':
                key = tuple(key)
                result.setdefault(key, [])
                result[key].append((base, confidence, depth))
    return result


def consensus_read(data):
    consistent_bases = []
    coverages = [len(v[0]) for k, v in data.items()]
    top = Counter(coverages).most_common(3)
    for (pos, ref, c), (bases, quals, reads, insertions, overlap) in data.items():
        coverages.append(len(bases))
        if Counter(overlap).most_common(1)[0][0] == 1:
            mean_depth = top[0][0]/2
        else:
            mean_depth = top[0][0]
        consistent_bases.append(consensus_base(bases, quals, insertions, mean_depth, c, pos, ref))
    return consistent_bases, statistics.median(coverages), top


def call_variant(result, out='mutation.txt', min_reads=2, min_conf=5, min_raw_reads=5):
    # call variant
    ordered_keys = sorted(result.keys(), key=lambda x:(x[0], x[1], x[2]))
    f = open(out, 'w')
    for key in ordered_keys:
        base_info = result[key]
        contig, position, ref_seq = key
        bases = [x[0] for x in base_info]
        depth = len(bases)
        base_counter = Counter(bases)
        for base, freq in base_counter.items():
            if base.upper() != ref_seq.upper():
                # print(base_counter)
                ad = (depth, freq)
                af = freq / depth
                confidences = [x[1][0] for x in base_info if x[0] == base]
                covs = [x[2] for x in base_info if x[0] == base]
                # 记录每个突变背后支持的证据情况
                if len(covs) >= min_reads and sum(confidences) >= min_conf and sum(covs) >= min_raw_reads:
                    # 5意味着至少要有2个一致性比较高的base支持
                    lst = (contig, position + 1, ref_seq, base, ad, af, confidences, covs)
                    f.write('\t'.join(str(x) for x in lst)+'\n')


def run_all(primers, bam, read_type=64, min_bq=5, cores=6, out='mutation.txt', min_reads=2, min_conf=5,
            genome='/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta'):
    primers = [x.strip() for x in open(primers)]
    cores = len(primers) if len(primers) <= cores else cores
    if cores <= 1:
        result = dict()
        for primer in primers:
            result.update(consensus_reads(bam, primer, read_type, min_bq, genome))
    else:
        get_consensus_reads = partial(consensus_reads, bam, read_type=read_type, min_bq=min_bq, genome=genome)
        with ProcessPoolExecutor(cores) as pool:
            result_list = pool.map(get_consensus_reads, primers)
            result = dict()
            _ = [result.update(x) for x in result_list]

    call_variant(result, out, min_reads, min_conf)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
