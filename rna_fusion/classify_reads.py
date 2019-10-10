"""
1. 从非clipped reads中中挑选可能包含junction位点或者融合断点的reads。
符合下述任何一个条件的read将作为2获得的junction位点或断点的辅助证据。
尤其，当spanning reads没有对应的split reads组合，他们就是主要证据。
    a. 从左到右增加窗口大小，窗口大小起始值4，从1开始计数，每增加一次，判断错配个数大于等于窗口增长次数，则停止滑动，收录该reads。
    b. 同上，从右向左滑窗
2. 挑选 split read即clipped reads，这些reads极有可能包含junction位点或者融合断点的reads
3. 定义spanning reads：r1和r2比对到junction或者断点两边的reads
4. 如果spanning reads比对到不同的基因，可以认为是支持基因间融合的reads
5. 如果spanning reads比对到同一个基因内，可以认为是支持某种可变剪接的reads

已知junction位点应该用两个坐标表示

"""
import pysam
import re
import itertools
import logging
import difflib
import statistics
import time


def timer(func):
    def inner(*args, **kw):
        start_time = time.time()
        result = func(*args,**kw)
        print(f'Runtime of {func.__name__} is', round(time.time()-start_time), 's')
        return result
    return inner


def set_logger(name='log.info', logger_id='x'):
    logger = logging.getLogger(logger_id)
    fh = logging.FileHandler(name, mode='w+')
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    return logger


logger = set_logger('jc_from_pair.log', logger_id='global')


def get_suspicious_reads(bam_file, fasta_file):
    """
    提取疑似包含新的可变剪接位点或包含融合断点的reads，可疑等级较低。
    大概过程是：根据cigar字符串找到和参考基因组完全匹配的read，然后判断前n个碱基或后n个碱基的错配情况，
    如果错配超过一定数目，则认为他们可能包含新的剪接位点或融合断点
    :param bam_file:
    :param fasta_file:
    :return:
    """
    ref = pysam.FastaFile(fasta_file)
    bam = pysam.AlignmentFile(bam_file, "rb")
    suspicious_left = list()
    suspicious_right = list()
    pattern = re.compile(r'[\dMNDI]+')
    n = 0
    for r in bam:
        n += 1
        if not r.is_paired:
            continue
        cigar = r.cigarstring
        if pattern.fullmatch(cigar):
            r_seq = r.query_sequence
            blocks = r.get_blocks()

            # 判断左边前7个碱基的错配情况
            left_seq = ref.fetch(r.reference_name, *blocks[0])
            b_id = 0
            while len(left_seq) < 7:
                b_id += 1
                left_seq += ref.fetch(r.reference_name, *blocks[b_id])
            for i in range(len(left_seq)):
                if i > 5:
                    break
                if len(set(left_seq[:i+3]) - set(r_seq[:i+3])) >= i + 1:
                    print(cigar, ':', r.query_name)
                    print('refl: ', left_seq[:i+3])
                    print('read: ', r_seq[:i+3])
                    # print(ref.fetch(r.reference_name, r.reference_start, r.reference_start + r.reference_length))
                    suspicious_left.append(r.to_string())
                    break

            # 判断右边前7个碱基的错配情况
            right_seq = ref.fetch(r.reference_name, *blocks[-1])
            b_id = -1
            while len(right_seq) < 7:
                b_id += -1
                right_seq += ref.fetch(r.reference_name, *blocks[b_id])
            for i in range(len(right_seq)):
                if i > 5:
                    break
                if len(set(right_seq[-i-3:]) - set(r_seq[-i-3:])) >= i + 1:
                    print(cigar, ':', r.query_name)
                    print('refr: ', right_seq[-i - 3:])
                    print('read: ', r_seq[-i - 3:])
                    suspicious_right.append(r.to_string())
                    break

        if len(suspicious_right)+len(suspicious_left) > 50:
            print(f'parsed {n} reads')
            break

    ref.close()
    bam.close()
    print(len(suspicious_right), len(suspicious_left))
    return suspicious_left, suspicious_right


def reverse_complement(seq):
    """
    :param seq: 输入序列
    :return: 返回reversed的序列
    """
    seq = seq.upper()
    complement = dict(zip(list("ATCG"), list("TAGC")))
    return ''.join(complement[base] if base in complement else base for base in seq[::-1])


@timer
def splitted_reads(bam_file, min_clip_size=3):
    bam = pysam.AlignmentFile(bam_file, "rb")
    splits = dict()
    for n, r in enumerate(bam):
        if not r.is_paired:
            continue
        cigar = r.cigarstring
        if r.is_duplicate or (not cigar) or ('H' in cigar or 'P' in cigar):
            continue
        cigar_tuples = r.cigartuples
        left_clip = cigar_tuples[0][0] == 4 and cigar_tuples[0][1] > min_clip_size
        right_clip = cigar_tuples[-1][0] == 4 and cigar_tuples[-1][1] > min_clip_size
        if left_clip or right_clip:
            record = splits.setdefault(r.query_name + '_' + str(int(r.is_read1)), [])
            # record[r.get_tag('HI')] = r  # 不能作为多重比对的唯一索引，尤其chimeric alignment
            record.append(r)
    else:
        bam.close()
    return splits


@timer
def paired_splitted_reads(bam_file, splits):
    bam = pysam.AlignmentFile(bam_file, "rb")
    targets = dict()
    read_names = set(x[:-2] for x in splits.keys())
    for n, r in enumerate(bam):
        if not r.is_paired:
            continue
        if r.query_name in read_names:
            record = targets.setdefault(r.query_name + '_' + str(int(r.is_read2)), [])
            record.append(r)
    bam.close()
    # targets = {k: v for k, v in targets.items() if k[:-2]+'_0' in targets and k[:-2]+'_1' in targets}
    return targets


def determine_fusion_type(a, a2, libtype='ISR'):
    break_type = 'fusion'
    if a.reference_name != a2.reference_name:
        break_type = 'fusion'
    elif libtype == 'ISR':
        if a2.is_reverse and not a.is_reverse:
            if a2.reference_start < a.reference_start:
                break_type = 'duplication or translocation'
        elif not a2.is_reverse and a.is_reverse:
            if a2.reference_start > a.reference_start:
                break_type = 'duplication or translocation'
        elif not a2.is_reverse and not a.is_reverse:
            if a.reference_start - a2.reference_end > a.query_length + 1:
                break_type = 'inversion'
        elif a2.is_reverse and a.is_reverse:
            if a.reference_start - a2.reference_end > a2.query_length + 1:
                break_type = 'inversion'
    return break_type


def filter_by_seq_similarity(ref_obj, chr1, break1, chr2, break2, reverse=False,
                             cutoff=0.8, extend_up=80, extend_down=80):
    try:
        seq = ref_obj.fetch(chr1, break1-extend_up, break1+extend_down)
        seq2 = ref_obj.fetch(chr2, break2-extend_up, break2+extend_down)
    except KeyError:
        chr1 = chr1[3:]
        chr2 = chr2[3:]
        try:
            seq = ref_obj.fetch(chr1, break1 - extend_up, break1 + extend_down)
            seq2 = ref_obj.fetch(chr2, break2 - extend_up, break2 + extend_down)
        except:
            return True

    # 通过大小写判断是否在重复区域附近
    seq_in_repeat = sum(x.islower() for x in seq) / len(seq) >= 0.2
    seq2_in_repeat = sum(x.islower() for x in seq2) / len(seq2) >= 0.2
    if seq_in_repeat or seq2_in_repeat:
        return True

    # 判断两个断点附近的序列是否非常相似或者包含的参考序列N太多
    if reverse:
        seq = reverse_complement(seq)
    ratio = difflib.SequenceMatcher(None, seq.lower(), seq2.lower()).ratio()
    if ratio > cutoff or seq.upper().count('N') > len(seq)*0.1 or seq2.upper().count('N') > len(seq2)*0.1:
        return True


def filter_by_coverage(bam_obj, chr_name, start, direction='keep_up', extend=23):
    # 判断断点上游或下游是否有reads覆盖，如果几乎没有read覆盖，则极有可能是一个假的断点
    try:
        if direction == 'keep_up':
            cov = len(list(bam_obj.fetch(chr_name, min([start-extend, 0]), start-3)))
        else:
            cov = len(list(bam_obj.fetch(chr_name, start+3, max([start+extend]))))
        if cov < 4:
            return True
    except KeyError:
        return True
    return False


@timer
def junction_from_single_end(splits, ref_fasta, masked_ref_fasta=None, min_clip_size=15, min_intron_size=30):
    """
    仅仅根据单条read的多重比对结果寻找断点，如一条read有两个比对结果:30M120S和30S120M
    这个断点可以是由于新的可变剪切产生的，也可以是融合基因带来的比对效果
    :param splits: 所有可能的包含断点的reads，是来自split_reads函数的结果
    :param min_intron_size: 最小的内含子大小，方便识别新的可变剪切
    :return:
    """
    sj = dict()
    ref = pysam.FastaFile(ref_fasta)
    if masked_ref_fasta:
        masked_ref = pysam.FastaFile(masked_ref_fasta)
    else:
        masked_ref = ref
    for name, aligns in splits.items():
        # 一条read可能有多个比对结果, 排列组合查看
        for a, a2 in itertools.permutations(aligns, 2):
            if (abs(a.reference_end - a2.reference_start) >= min_intron_size
                    or a.reference_name != a2.reference_name):
                if a.is_reverse == a2.is_reverse:
                    if (a.cigartuples[-1][0] == 4
                            and a2.cigartuples[0][0] == 4
                            and len(a.seq[a.query_alignment_end:]) > min_clip_size
                            and len(a2.seq[:a2.query_alignment_start]) > min_clip_size
                            and a.seq[a.query_alignment_end:] == a2.seq[a2.query_alignment_start:]
                            and a.seq[:a.query_alignment_end] == a2.seq[:a2.query_alignment_start]
                            and not filter_by_seq_similarity(masked_ref, a.reference_name, a.reference_end,
                                                             a2.reference_name, a2.reference_start)
                    ):
                        break1 = a.reference_name + ':' + str(a.reference_end)
                        break2 = a2.reference_name + ':' + str(a2.reference_start)
                        break_pair = break1 + ' ' + break2
                        sj.setdefault(break_pair, set())
                        sj[break_pair].add(a.seq)

                # else:
                #     if (a.cigartuples[-1][0] == 4
                #             and a2.cigartuples[-1][0] == 4 and
                #             len(a.seq[a.query_alignment_end:]) > min_clip_size
                #             and len(a2.seq[a2.query_alignment_end:]) > min_clip_size
                #             and a.seq[:a.query_alignment_end] == reverse_complement(a2.seq[a2.query_alignment_end:])
                #             and a.seq[a.query_alignment_end:] == reverse_complement(a2.seq[:a2.query_alignment_end])
                #             and not filter_by_seq_similarity(ref, a.reference_name, a.reference_end,
                #                                              a2.reference_name, a2.reference_end, reverse=True)
                #     ):
                #         break1 = a.reference_name + ':' + str(a.reference_end)
                #         break2 = a2.reference_name + ':' + str(a2.reference_end)
                #         break_pair = break1 + ' ' + break2
                #         sj.setdefault(break_pair, set())
                #         sj[break_pair].add(a.seq)
    ref.close()
    result = dict()
    for break_pair, reads in sj.items():
        result[break_pair] = len(reads)
        # result[break_pair] = reads
    return dict(sorted(result.items(), key=lambda x: x[1], reverse=True))


@timer
def junction_from_one_pair(split_pair_reads, ref_fasta, max_inner_dist=250, libtype="ISR"):
    ref = pysam.FastaFile(ref_fasta)
    sj = dict()
    fusion_type = dict()
    names = [x[:-2] for x in split_pair_reads.keys() if x.endswith('0')]
    single_support = dict()
    for name in names:
        r1_name, r2_name = name + '_0', name + '_1'
        if (r1_name not in split_pair_reads) or (r2_name not in split_pair_reads):
            continue
        for a, a2 in itertools.product(split_pair_reads[r1_name], split_pair_reads[r2_name]):
            if a.reference_start != a2.next_reference_start:
                continue
            if not a.cigartuples or (not a2.cigartuples):
                continue

            # 如果r1和r2有overlap序列，而且overlap包含断点
            if a.cigartuples[-1][0] == a2.cigartuples[0][0] == 4:
                evidence1 = a.seq[:a.query_alignment_end].endswith(a2.seq[:a2.query_alignment_start])
                evidence2 = a2.seq[a2.query_alignment_start:].startswith(a.seq[a.query_alignment_end:])
                if evidence1 and evidence2:
                    if not filter_by_seq_similarity(
                        ref, a.reference_name, a.reference_end,
                        a2.reference_name, a2.reference_start
                    ):
                        break1 = a.reference_name + ':' + str(a.reference_end)
                        break2 = a2.reference_name + ':' + str(a2.reference_start)
                        break_pair = break1 + ' ' + break2
                        sj.setdefault(break_pair, set())
                        sj[break_pair].add(a.seq+a2.seq)
                        fusion_type[break_pair] = determine_fusion_type(a, a2, libtype)

            # 仅r1包含断点，而r2不包含断点
            elif a.cigartuples[-1][0] == 4 and a2.cigartuples[0][0] == 0:
                clip_seq = a.seq[a.query_alignment_end:]
                # 如果 r1和r2有overlap序列，而且overlap不包含断点，但是r1包含断点
                # 假设overlap至少包含6个碱基，那么r2前6个碱基一定包含在r1右端被clip掉的reads中
                if len(clip_seq) >= 6:
                    # 假设是r1和r2存在overlap,但overlap不包含断点
                    ind = clip_seq.rfind(a2.seq[:6])
                    if ind >= 0:
                        b2_pos = a2.reference_start - ind
                        if abs(b2_pos - a.reference_end) < 3 and a.reference_name == a2.reference_name:
                            logger.info(f'read1: {a.to_string()}')
                            logger.info(f'read2: {a2.to_string()}')
                            logger.info(f"以上两条read可能由于多了几个错配而发生clip, 最后反推回来并没有发现断点")
                        else:
                            if not filter_by_seq_similarity(
                                    ref, a.reference_name, a.reference_end,
                                    a2.reference_name, a2.reference_start
                            ):
                                break1 = a.reference_name + ':' + str(a.reference_end)
                                break2 = a2.reference_name + ':' + str(b2_pos)
                                break_pair = break1 + ' ' + break2
                                sj.setdefault(break_pair, set())
                                sj[break_pair].add(a.seq+a2.seq)
                                fusion_type[break_pair] = determine_fusion_type(a, a2, libtype)

                if len(clip_seq) >= 8:
                    # 如果r1和r2没有overlap，但是r1包含断点，且假设r1右端clip掉的seq和r2之间的真实距离在一定范围内
                    # 根据r1右端clip掉的seq一定在r2的上游区域，定位其在上游区域的位置得到的就是第二个断点
                    # print(f"在r2上游的{max_inner_dist}范围内搜索 {clip_seq}")
                    r2_up_seq = ref.fetch(a2.reference_name,
                                          max([a2.reference_start-max_inner_dist, 0]),
                                          a2.reference_start)
                    ind = r2_up_seq.rfind(clip_seq)
                    if ind >= 0:
                        b2_pos = a2.reference_start-max_inner_dist+ind
                        if abs(b2_pos - a.reference_end) < 3 and a.reference_name == a2.reference_name:
                            logger.info(f'read1: {a.to_string()}')
                            logger.info(f'read2: {a2.to_string()}')
                            logger.info(f"以上两条read可能由于多了几个错配而发生clip, 最后反推回来并没有发现断点")
                        else:
                            if not filter_by_seq_similarity(
                                    ref, a.reference_name, a.reference_end,
                                    a2.reference_name, a2.reference_start
                            ):
                                break1 = a.reference_name + ':' + str(a.reference_end)
                                break2 = a2.reference_name + ':' + str(b2_pos)
                                break_pair = break1 + ' ' + break2
                                sj.setdefault(break_pair, set())
                                sj[break_pair].add(a.seq+a2.seq)
                                fusion_type[break_pair] = determine_fusion_type(a, a2, libtype)

            # 下面假设r1不包含断点，而是r2包含断点
            elif a.cigartuples[-1][0] == 0 and a2.cigartuples[0][0] == 4:
                clip_seq = a2.seq[:a.query_alignment_start]
                # 假设overlap至少包含6个碱基，那么clip_seq中前6个碱基一定在r1中
                if len(clip_seq) >= 6:
                    # 如果 r1和r2有overlap序列，而且overlap不包含断点，但是r2包含断点
                    ind = a.seq.rfind(clip_seq[:6])
                    if ind >= 0:
                        b1_pos = a.reference_end - (a.query_length - ind) + len(clip_seq)
                        if abs(b1_pos - a2.reference_start) < 3 and a.reference_name == a2.reference_name:
                            logger.info(f'read1: {a.to_string()}')
                            logger.info(f'read2: {a2.to_string()}')
                            logger.info(f"以上两条read可能由于多了几个错配而发生clip, 最后反推回来并没有发现断点")
                        else:
                            if not filter_by_seq_similarity(
                                    ref, a.reference_name, a.reference_end,
                                    a2.reference_name, a2.reference_start
                            ):
                                break1 = a.reference_name + ':' + str(b1_pos)
                                break2 = a2.reference_name + ':' + str(a2.reference_start)
                                break_pair = break1 + ' ' + break2
                                sj.setdefault(break_pair, set())
                                sj[break_pair].add(a.seq+a2.seq)
                                fusion_type[break_pair] = determine_fusion_type(a, a2, libtype)

                if len(clip_seq) >= 8:
                    # 如果r1和r2没有overlap，但是r2包含断点，且假设r2左端clip掉的seq和r1之间的真实距离在一定范围内
                    # 根据r2左端clip掉的seq一定在r1的下游区域，定位其在下游区域的位置得到的就是第二个断点
                    # print(f"在r1上游的{max_inner_dist}范围内搜索 {clip_seq}")
                    r1_down_seq = ref.fetch(a.reference_name,
                                            a.reference_end,
                                            a.reference_end+max_inner_dist)
                    ind = r1_down_seq.find(clip_seq)
                    if ind >= 0:
                        b1_pos = a.reference_end + ind
                        if abs(b1_pos - a2.reference_start) < 3 and a.reference_name == a2.reference_name:
                            logger.info(f'read1: {a.to_string()}')
                            logger.info(f'read2: {a2.to_string()}')
                            logger.info(f"以上两条read可能由于多了几个错配而发生clip, 最后反推回来并没有发现断点")
                        else:
                            if not filter_by_seq_similarity(
                                    ref, a.reference_name, a.reference_end,
                                    a2.reference_name, a2.reference_start
                            ):
                                break1 = a.reference_name + ':' + str(b1_pos)
                                break2 = a2.reference_name + ':' + str(a2.reference_start)
                                break_pair = break1 + ' ' + break2
                                sj.setdefault(break_pair, set())
                                sj[break_pair].add(a.seq+a2.seq)
                                fusion_type[break_pair] = determine_fusion_type(a, a2, libtype)

            elif a.cigartuples[0][0] == 4 or a2.cigartuples[-1][0] == 4:
                if a.cigartuples[0][0] == 4:
                    key = a.reference_name+':'+str(a.reference_start)
                    single_support.setdefault(key, set()).add(a.seq)
                if a2.cigartuples[-1][0] == 4:
                    key = a2.reference_name+':'+str(a2.reference_end)
                    single_support.setdefault(key, set()).add(a2.seq)
    ref.close()
    result = {k:len(v) for k, v in sj.items()}
    fusion_paired = dict(sorted(result.items(), key=lambda x: x[1], reverse=True))

    single_support = {k: len(v) for k, v in single_support.items()}
    single_support = dict(sorted(single_support.items(), key=lambda x: x[1], reverse=True))
    return fusion_paired, fusion_type, single_support


@timer
def junction_from_two_pairs(split_pairs, ref_fasta, bam, min_intro_size=30, r1_on_forward=True):
    ref_obj = pysam.FastaFile(ref_fasta)
    read_names = [x[:-2] for x in split_pairs.keys() if x.endswith('0')]
    jac = dict()
    hang_len = 7
    for name in read_names:
        r1_name, r2_name = name + '_0', name + '_1'
        if (r1_name not in split_pairs) or (r2_name not in split_pairs):
            continue
        for r1, r2 in itertools.product(split_pairs[r1_name], split_pairs[r2_name]):
            # 判读是不是r1和r2是不是一对
            if r1.reference_start != r2.next_reference_start:
                continue
            if not r1.cigartuples or (not r2.cigartuples):
                continue
            # 寻找同链上的融合
            # 判断r1是不是直接比对到正链，如果不涉及不同链融合，此时要求r2必定是需要反向互补后才能比对到正链
            if r1.is_reverse != r2.is_reverse == r1_on_forward:
                if r1.cigartuples[0][0] == 4:
                    if len(r1.seq[:r1.query_alignment_start]) >= hang_len:
                        junction_around = r1.seq[r1.query_alignment_start-hang_len:r1.query_alignment_start+hang_len]
                        jac.setdefault(junction_around, set()).add(
                            (r1.reference_name, r1.reference_start, r1.seq, 'keep_down', r1.reference_end)
                        )
                if r1.cigartuples[-1][0] == 4:
                    if len(r1.seq[r1.query_alignment_end:]) >= hang_len:
                        junction_around = r1.seq[r1.query_alignment_end-hang_len:r1.query_alignment_end+hang_len]
                        jac.setdefault(junction_around, set()).add(
                            (r1.reference_name, r1.reference_end, r1.seq, 'keep_up', r1.reference_start)
                        )
                if r2.cigartuples[0][0] == 4:
                    if len(r2.seq[:r2.query_alignment_start]) >= hang_len:
                        junction_around = r2.seq[r2.query_alignment_start-hang_len:r2.query_alignment_start+hang_len]
                        jac.setdefault(junction_around, set()).add(
                            (r2.reference_name, r2.reference_start, r2.seq, 'keep_down', r2.reference_end)
                        )
                if r2.cigartuples[-1][0] == 4:
                    if len(r2.seq[r2.query_alignment_end:]) >= hang_len:
                        junction_around = r2.seq[r2.query_alignment_end - hang_len:r2.query_alignment_end + hang_len]
                        jac.setdefault(junction_around, set()).add(
                            (r2.reference_name, r2.reference_end, r2.seq, 'keep_up', r2.reference_start)
                        )
    fusion_dict = dict()
    bam_obj = pysam.AlignmentFile(bam)
    for k, v in jac.items():
        if len(v) < 3:
            continue
        for p, p2 in itertools.combinations(v, 2):
            if (p[0] != p2[0] or abs(p[1] - p2[1]) > min_intro_size) and p[3] != p2[3]:
                if filter_by_seq_similarity(ref_obj, p[0], p[1], p2[0], p2[1]):
                    continue
                # if filter_by_coverage(bam_obj, p[0], p[4], direction=p[3]):
                #     continue
                # if filter_by_coverage(bam_obj, p2[0], p2[4], direction=p2[3]):
                #     continue
                if p[3] == 'keep_up':
                    key = p[0] + ':' + str(p[1]) + ' ' + p2[0] + ':' + str(p2[1])
                else:
                    key = p2[0] + ':' + str(p2[1]) + ' ' + p[0] + ':' + str(p[1])
                fusion_dict.setdefault(key, 0)
                fusion_dict[key] += 1
    bam_obj.close()
    return dict(sorted(fusion_dict.items(), key=lambda x: x[1], reverse=True))



def fusion_pipe(bam_file, ref_fasta, max_inner_dist=300, libtype='ISR'):
    all_split_reads = splitted_reads(bam_file, min_clip_size=3)
    fusion_fse = junction_from_single_end(all_split_reads, ref_fasta)
    paired_split_reads = paired_splitted_reads(bam_file, all_split_reads)
    fusion_fpe, fusion_type, single_support = junction_from_one_pair(
        paired_split_reads, ref_fasta, max_inner_dist=max_inner_dist, libtype=libtype
    )
    # merge fusion
    merged_fusion = dict()
    for k, v in fusion_fpe.items():
        merged_fusion[k] = [v, 0, 0, 0]
        if k in fusion_fse:
            merged_fusion[k][1] = fusion_fse.pop(k)
    for k, v in fusion_fse.items():
        merged_fusion[k] = [0, v, 0, 0]
    for k, v in merged_fusion.items():
        break1, break2 = k.split()
        if break1 in single_support:
            v[2] = single_support.pop(break1)
        if break2 in single_support:
            v[3] = single_support.pop(break2)
    merged_fusion = dict(sorted(merged_fusion.items(), key=lambda x:sum(x[1]), reverse=True))
    for k, v in merged_fusion.items():
        if k in fusion_type:
            v.append(fusion_type[k])
        else:
            v.append('unknown')
    return merged_fusion, all_split_reads


# if __name__ == '__main__':
#     from xcmds import xcmds
#     xcmds.xcmds(locals(), include=['get_clipped_reads','get_suspicious_reads', 'single_junction_pos'])
#
#
