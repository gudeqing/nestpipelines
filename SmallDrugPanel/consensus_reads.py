import re
import statistics
from collections import Counter
from functools import partial
import logging
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Pool, Manager
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
    elif top1_ratio >= 0.75:
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
        # 用<>表示整体的插入序列, 用()表示deletion
        represent = '<' + represent + '>'
    rep_len = len(represent)
    return represent, [rep_qual]*rep_len, [confidence]*rep_len, support_depth, contig, position, ref_seq


def format_consensus_bases(bqc):
    b, q, c = bqc
    if b.startswith('<'):
        return b[1:-1], ''.join([chr(q[0]+33)] * (len(b) - 2)), ''.join([str(c[0])] * (len(b) - 2))
    elif b.startswith(('S', '(')):
        return b[1], chr(q[1]+33), str(c[1])
    else:
        return b[0], chr(q[0]+33), str(c[0])


def consensus_reads(bam, primer, read_type=64, min_bq=0, fq_lst=None,
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
        ignore_overlaps=False,
        # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
        max_depth=300000,
        fastafile=None,
        flag_require=read_type,
        flag_filter=4,
    )
    # 在pos处带入ref信息
    genome = pysam.FastaFile(genome)

    ## 初始化pileup存储容器
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
                insertion = re.sub('\+\d+',  '', base)
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
                base = ref
            elif base == '':
                # 如果该位点没有覆盖，该处的base为''
                base = 'X'
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
            # 某个分组在目标区域没有对应的reads?, 之前发现是pileup时参数max_depth不够大导致漏掉
            print(f'{group_name} is empty, but raw group_size is {len(group2read[group_name])}?')
            if len(group2read[group_name]) > 5:
                print(group2read[group_name])
            continue
        consistent_bases, median_cov, top = consensus_read(data)
        print(f'>{group_name}')
        print(f'median coverage is {median_cov}')
        print(f'top3 coverage frequency is {top}')
        consensus_seq = ''.join(x[0] for x in consistent_bases).strip('X')
        print(consensus_seq)
        # print(f'consensus sequence length is {len(consensus_seq)}')
        # print(f'we deem there is {group2overlap[group_name]}overlap between read1 and read2')
        # print(consistent_bases)

        # 制作突变需要的字典
        # umi = group_name.rsplit(':', 1)[-1]
        for base, qual, confidence, alt_depth, *key in consistent_bases:
            if base[0] not in ['S', 'X']:
                # clipped 和没有read支持的位置不能call突变
                key = tuple(key)
                result.setdefault(key, [])
                result[key].append((base, confidence, alt_depth))

        # 为输出consensus reads做准备, fq_lst是一个可以在多进程间共享的的list
        if fq_lst is not None:
            # 想办法去掉两端的X
            consensus_seq = ''.join(x[0] for x in consistent_bases)
            lm = re.match('X+', consensus_seq)
            left_start = 0 if lm is None else len(lm.group())
            rm = re.match('X+', consensus_seq[::-1])
            right_end = len(consistent_bases) if rm is None else len(consistent_bases) - len(rm.group())

            # 整理fastq的信息
            bqc = [x[:3] for x in consistent_bases[left_start:right_end]]
            # 处理indel和clipped的标签信息
            bqc = list(map(format_consensus_bases, bqc))
            # 形成fastq信息
            seqs = ''.join(x[0] for x in bqc)
            quals = ''.join(x[1] for x in bqc)
            # confs = '+' + ''.join(x[2] for x in bqc)[1:]
            confs = '+'
            header = f'@{group_name} read_number:N:{int(median_cov)}:{group2overlap[group_name]}'
            fq_lst.append([read_type, header, seqs, confs, quals])
    else:
        if fq_lst is not None:
            fq_lst.append([primer])

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


def create_vcf(vcf_path, genome='hg19', chrom_name_is_numeric=False):
    vcf = pysam.VariantFile(vcf_path, 'w')
    contig_info = [
        '##fileformat=VCFv4.3',
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
        '##INFO=<ID=Confidences,Number=1,Type=String,Description="convince of level of consuensus">',
        '##INFO=<ID=RawAltNumber,Number=1,Type=String,Description="raw read number that support the alt">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">',
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
    ]
    if chrom_name_is_numeric:
        contig_info = [x.replace('=chr', '=') for x in contig_info]
    for line in contig_info:
        vcf.header.add_line(line)
    return vcf


def call_variant(result, out='mutation.vcf', min_umi_depth=5, min_reads=2, min_conf=4, min_raw_reads=5):
    # call variant
    ordered_keys = sorted(result.keys(), key=lambda x: (x[0], x[1], x[2]))
    sample = out.split('.', 1)[0]
    vcf = create_vcf(out)
    vcf.header.add_sample(sample)
    variant_number = 0
    for key in ordered_keys:
        base_info = result[key]
        contig, position, ref_seq = key
        bases = [x[0] for x in base_info]
        depth = len(bases)
        base_counter = Counter(bases)
        for base, freq in base_counter.items():
            if base.upper() != ref_seq:
                # print(base_counter)
                ad = (depth, freq)
                af = freq / depth
                confidences = [x[1][0] for x in base_info if x[0] == base]
                covs = [x[2] for x in base_info if x[0] == base]
                # umis = [x[3] for x in base_info if x[0] == base]
                # filtering
                if len(covs) >= min_reads and sum(confidences) >= min_conf \
                        and sum(covs) >= min_raw_reads and depth >= min_umi_depth:
                    # min_conf意味着至少要有2个一致性比较高的base支持
                    variant_number += 1
                    # format output
                    if base.startswith('('):
                        # deletion
                        ref = base[1:-1]
                        alt = ref_seq
                    elif base.startswith('<'):
                        # insertion
                        ref = ref_seq
                        alt = base[1:-1]
                    else:
                        ref = ref_seq
                        alt = base
                    # lst = (contig, position + 1, ref, alt, ad, af, confidences, covs)
                    info = dict(Confidences=str(confidences), RawAltNumber=str(covs))
                    record = vcf.new_record()
                    record.contig = contig
                    record.start = position
                    record.ref = ref
                    record.alts = [alt]
                    record.qual = None
                    record.filter.add('PASS')
                    record.info.update(info)
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


def run_all(primers, bam, read_type=64, cores=8, out='mutation.vcf', min_reads=2, min_conf=5, min_bq=5,
            genome='/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta'):
    primers = [x.strip() for x in open(primers)]
    cores = len(primers) if len(primers) <= cores else cores
    # 开拓进程之间的共享空间, 即使用一个进程间可以共享的list，
    # N-1个进程往list里添加信息，剩下一个进程从list清空信息并输出到同一个文件
    manager = Manager()
    fq_lst = manager.list()

    result = dict()
    if cores <= 1:
        for primer in primers:
            result.update(consensus_reads(bam, primer, read_type, min_bq, fq_lst, genome))
            write_fastq(fq_lst, 1)
    else:
        get_consensus_reads = partial(
            consensus_reads, bam, read_type=read_type, min_bq=min_bq, fq_lst=fq_lst, genome=genome
        )
        with ProcessPoolExecutor(cores) as pool:
            tasks = [pool.submit(write_fastq, fq_lst, len(primers))]
            for each in primers:
                tasks.append(pool.submit(get_consensus_reads, each))
            futures = [x.result() for x in tasks[1:]]
            tasks[0].result()
            _ = [result.update(x) for x in futures]

    call_variant(result, out, min_reads, min_conf)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
