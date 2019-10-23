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
import os
import pysam
import re
import itertools
import logging
from difflib import SequenceMatcher as similarity
import statistics
import time
from subprocess import check_call
import json

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
def splitted_reads(bam_file, min_clip_size=3, max_hit_num=10, logger_name='soft_clipped.pos.list'):
    break_points = dict()
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
        if left_clip:
            point = r.reference_name + ':' + str(r.reference_start)
            break_points.setdefault(point, 0)
            break_points[point] += 1
        right_clip = cigar_tuples[-1][0] == 4 and cigar_tuples[-1][1] > min_clip_size
        if right_clip:
            point = r.reference_name + ':' + str(r.reference_end)
            break_points.setdefault(point, 0)
            break_points[point] += 1
        if left_clip or right_clip:
            record = splits.setdefault(r.query_name + '_' + str(int(r.is_read1)), [])
            # record[r.get_tag('HI')] = r  # 不能作为多重比对的唯一索引，尤其chimeric alignment
            record.append(r)
    bam.close()
    with open(logger_name, 'w') as f:
        _ = [f.write(k+'\t'+str(v)+'\n') for k, v in break_points.items()]
    return {k: v for k, v in splits.items() if len(v) <= max_hit_num}


@timer
def hard_clipped_reads(splits, bam_file, max_hit_num=10, logger_name='hard_clipped.pos.list'):
    break_points = dict()
    bam = pysam.AlignmentFile(bam_file, "rb")
    hard_splits = dict()
    for n, r in enumerate(bam):
        if not r.is_paired or r.is_duplicate or (not r.cigarstring):
            continue
        cigar_tuples = r.cigartuples
        if cigar_tuples[0][0] == 5 or cigar_tuples[-1][0] == 5:
            record = r.query_name + '_' + str(int(r.is_read1))
            if record in splits:
                soft_aln = splits[record][0]
                if soft_aln.is_reverse == r.is_reverse:
                    r.seq = soft_aln.seq
                    r.query_sequence = soft_aln.query_sequence
                    r.query_qualities = soft_aln.query_qualities
                else:
                    seq = reverse_complement(soft_aln.query_sequence)
                    r.seq = seq
                    r.query_sequence = seq
                    r.query_qualities = soft_aln.query_qualities[::-1]
                if cigar_tuples[0][0] == 5:
                    point = r.reference_name + ':' + str(r.reference_start)
                    break_points.setdefault(point, 0)
                    break_points[point] += 1
                else:
                    point = r.reference_name + ':' + str(r.reference_end)
                    break_points.setdefault(point, 0)
                    break_points[point] += 1
                hard_splits.setdefault(record, []).append(r)
    bam.close()
    with open(logger_name, 'w') as f:
        _ = [f.write(k+'\t'+v+'\n') for k, v in break_points.items()]
    return {k: v for k, v in hard_splits.items() if len(v) <= max_hit_num}


@timer
def junction_around_dict(splits, min_clip=6):
    # 2. get {junction—around-seq: [support-alignments]}
    jad = dict()
    for read, alignments in splits.items():
        for ind, a in enumerate(alignments):
            # start index of the aligned query portion of the sequence (0-based, inclusive)
            # end index of the aligned query portion of the sequence (0-based, exclusive)
            if a.cigartuples[-1][0] == 4 and a.cigartuples[-1][1] >= min_clip:
                if a.query_alignment_end > min_clip:
                    jun_seq = a.seq[a.query_alignment_end - min_clip:a.query_alignment_end + min_clip]
                    if len(set(jun_seq)) > 2:  # 如果jun_seq是高度重复序列则去掉
                        jad.setdefault(jun_seq, list()).append((read, ind, a.reference_name, a.reference_end))
            if a.cigartuples[0][0] == 4 and a.cigartuples[0][1] >= min_clip:
                if a.query_alignment_start+min_clip < a.query_length:
                    jun_seq = a.seq[a.query_alignment_start - min_clip:a.query_alignment_start + min_clip]
                    if len(set(jun_seq)) > 2:
                        jad.setdefault(jun_seq, list()).append((read, ind, a.reference_name, a.reference_start))

    # 3. 如果jun_seq之间反向互补则认为他们包含同一个断点
    for jun_seq in list(jad.keys()):
        if jun_seq in jad:
            jun_seq_reverse = reverse_complement(jun_seq)
            if jun_seq_reverse in jad:
                # 回文序列反向互补为其本尊
                if jun_seq_reverse != jun_seq:
                    info = jad.pop(jun_seq_reverse)
                    jad[jun_seq].extend(info)
    return jad


@timer
def paired_splitted_reads(bam_file, splits, max_hit_num=20):
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
    # 过滤掉没有配对的，过滤掉比对的位置超过10的reads
    targets = {k: v for k, v in targets.items() if len(v) < max_hit_num and k[:-2]+'_0' in targets and k[:-2]+'_1' in targets}
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


def filter_by_seq_similarity(ref_obj, chr1, break1, chr2, break2, ref_masked=True, cutoff=0.85, extend=50):
    up_seq = max(0, break1-extend)
    up_seq2 = max(0, break2-extend)
    try:
        seq = ref_obj.fetch(chr1, up_seq, break1+extend)
        seq2 = ref_obj.fetch(chr2, up_seq2, break2+extend)
    except KeyError:
        chr1 = chr1[3:]
        chr2 = chr2[3:]
        try:
            seq = ref_obj.fetch(chr1, up_seq, break1 + extend)
            seq2 = ref_obj.fetch(chr2, up_seq2, break2 + extend)
        except KeyError:
            return True
        except ValueError:
            return False

    # 通过大小写判断是否在重复区域附近
    if ref_masked:
        seq_in_repeat = sum(x.islower() for x in seq) / len(seq) >= 0.35
        seq2_in_repeat = sum(x.islower() for x in seq2) / len(seq2) >= 0.35
        if seq_in_repeat or seq2_in_repeat:
            return True

    # 判断两个断点附近的序列是否非常相似或者包含的参考序列N太多
    seq = seq.upper()
    seq2 = seq2.upper()
    if seq.count('N') > len(seq)*0.1 or seq2.count('N') > len(seq2)*0.1:
        return True
    # 可以过滤掉假基因或相似基因的融合候选
    for seq in [seq, reverse_complement(seq)]:
        if similarity(None, seq[:extend], seq2[:extend]).ratio() > cutoff:
            return True
        elif similarity(None, seq[:extend], seq2[extend:]).ratio() > cutoff:
            return True
        elif similarity(None, seq[extend:], seq2[:extend]).ratio() > cutoff:
            return True
        elif similarity(None, seq[extend:], seq2[extend:]).ratio() > cutoff:
            return True


def get_overlap(seq, seq2, min_overlap=6, maximum=True):
    """寻找第一条序列尾部和第二条序列的头部的重合序列，默认求最大overlap"""
    match = re.match(seq2[:min_overlap-2]+'.*'+seq[-2:], seq2)
    if match:
        m_seq = match.group()
        indexes = range(min_overlap, len(m_seq)+1)
        if maximum:
            indexes = indexes[::-1]
        for i in indexes:
            if seq.endswith(m_seq[:i]):
                return m_seq[:i]
    return ''


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
                elif a.is_reverse != a2.is_reverse == True:
                    if (a.cigartuples[-1][0] == 4
                            and a2.cigartuples[-1][0] == 4 and
                            len(a.seq[a.query_alignment_end:]) > min_clip_size
                            and len(a2.seq[a2.query_alignment_end:]) > min_clip_size
                            and a.seq[:a.query_alignment_end] == reverse_complement(a2.seq[a2.query_alignment_end:])
                            and a.seq[a.query_alignment_end:] == reverse_complement(a2.seq[:a2.query_alignment_end])
                            and not filter_by_seq_similarity(masked_ref, a.reference_name, a.reference_end,
                                                             a2.reference_name, a2.reference_end, reverse=True)
                    ):
                        break1 = a.reference_name + ':' + str(a.reference_end)
                        break2 = a2.reference_name + ':' + str(a2.reference_end)
                        break_pair = break1 + ' ' + break2
                        sj.setdefault(break_pair, set())
                        sj[break_pair].add(a.seq)
    ref.close()
    result = dict()
    for break_pair, reads in sj.items():
        result[break_pair] = len(reads)
        # result[break_pair] = reads
    return dict(sorted(result.items(), key=lambda x: x[1], reverse=True))


@timer
def junction_from_one_pair(split_pair_reads, ref_fasta, masked_ref_fasta=None, min_clip_size=7,
                           max_inner_dist=250, libtype="ISR", logger=None):
    logger = set_logger('jc_from_one_pair.log', logger_id='jfop') if logger is None else logger
    ref = pysam.FastaFile(ref_fasta)
    if masked_ref_fasta:
        masked_ref = pysam.FastaFile(masked_ref_fasta)
    else:
        masked_ref = ref
    sj = dict()
    fusion_type = dict()
    names = [x[:-2] for x in split_pair_reads.keys() if x.endswith('0')]
    for name in names:
        r1_name, r2_name = name + '_0', name + '_1'
        if (r1_name not in split_pair_reads) or (r2_name not in split_pair_reads):
            continue
        r1r2_pair = itertools.product(split_pair_reads[r1_name], split_pair_reads[r2_name])
        r2r1_pair = itertools.product(split_pair_reads[r2_name], split_pair_reads[r1_name])
        for a, a2 in itertools.chain(r1r2_pair, r2r1_pair):
            if a.reference_start != a2.next_reference_start:
                continue
            if not a.cigartuples or (not a2.cigartuples):
                continue

            # 如果r1和r2有overlap序列，而且overlap包含断点
            if (a.cigartuples[-1][0] == a2.cigartuples[0][0] == 4
                    and a.is_reverse != a2.is_reverse
                    and (a.reference_end < a2.reference_start or a.reference_name != a2.reference_name)
                    and get_overlap(a.seq, a2.seq, min_clip_size*2)
                    and a.seq[:a.query_alignment_end].endswith(a2.seq[:a2.query_alignment_start])
                    and a2.seq[a2.query_alignment_start:].startswith(a.seq[a.query_alignment_end:])
                    and not filter_by_seq_similarity(
                        masked_ref, a.reference_name, a.reference_end,
                        a2.reference_name, a2.reference_start)
            ) or (a.cigartuples[-1][0] == a2.cigartuples[0][0] == 4
                    and a.is_reverse == a2.is_reverse
                    # and (a.reference_end < a2.reference_start or a.reference_name != a2.reference_name)
                    and get_overlap(a.seq, reverse_complement(a2.seq), min_clip_size*2)
                    and a.seq[:a.query_alignment_end].endswith(reverse_complement(a2.seq[:a2.query_alignment_start]))
                    and a2.seq[a2.query_alignment_start:].startswith(reverse_complement(a.seq[a.query_alignment_end:]))
                    and not filter_by_seq_similarity(
                        masked_ref, a.reference_name, a.reference_end,
                        a2.reference_name, a2.reference_start, reverse=True)
            ):
                break1 = a.reference_name + ':' + str(a.reference_end)
                break2 = a2.reference_name + ':' + str(a2.reference_start)
                break_pair = break1 + ' ' + break2
                sj.setdefault(break_pair, set())
                sj[break_pair].add(a.seq+a2.seq)
                fusion_type[break_pair] = determine_fusion_type(a, a2, libtype)
                if a.is_reverse == a2.is_reverse:
                    logger.info(break_pair + ' diff_strand '+' Overlap contain break points '+ a.cigarstring+' '+a2.cigarstring)
                    logger.info(a.seq)
                    logger.info(a2.seq)
                else:
                    logger.info(break_pair + ' same_strand '+' Overlap contain break points '+ a.cigarstring+' '+a2.cigarstring)
                    logger.info(a.seq)
                    logger.info(a2.seq)

            # (1) 仅r1包含断点，而r2不包含相应断点
            elif a.cigartuples[-1][0] == 4 and a2.cigartuples[0][0] == 0:
                clip_seq = a.seq[a.query_alignment_end:]
                if len(clip_seq) >= min_clip_size and get_overlap(a.seq, a2.seq, min_clip_size):
                    # 假设是r1和r2存在overlap,但overlap不包含断点
                    # 假设overlap至少包含6个碱基，那么r2前6个碱基一定包含在r1右端被clip掉的reads中
                    ind = clip_seq.rfind(reverse_complement(a2.seq[:min_clip_size]))
                    if ind >= 0:
                        b2_pos = a2.reference_start - ind
                        if abs(b2_pos - a.reference_end) < 3 and a.reference_name == a2.reference_name:
                            logger.info(f'read1: {a.to_string()}')
                            logger.info(f'read2: {a2.to_string()}')
                            logger.info(f"以上两条read最后反推回来并没有发现断点")
                        else:
                            if not filter_by_seq_similarity(
                                    masked_ref, a.reference_name, a.reference_end,
                                    a2.reference_name, a2.reference_start
                            ):
                                break1 = a.reference_name + ':' + str(a.reference_end)
                                break2 = a2.reference_name + ':' + str(b2_pos)
                                break_pair = break1 + ' ' + break2
                                sj.setdefault(break_pair, set())
                                sj[break_pair].add(a.seq+a2.seq)
                                fusion_type[break_pair] = determine_fusion_type(a, a2, libtype)

                elif len(clip_seq) >= min_clip_size and not get_overlap(a.seq, a2.seq, 5):
                    # 如果r1和r2没有overlap，但是r1包含断点，且假设r1右端clip掉的seq和r2之间的真实距离在一定范围内
                    # 根据r1右端clip掉的seq一定在r2的上游区域，定位其在上游区域的位置得到的就是第二个断点
                    # print(f"在r2上游的{max_inner_dist}范围内搜索 {clip_seq}")
                    r2_up_seq = ref.fetch(a2.reference_name,
                                          max([a2.reference_start-max_inner_dist, 0]),
                                          a2.reference_start)
                    ind = r2_up_seq.rfind(reverse_complement(clip_seq))
                    if ind >= 0:
                        b2_pos = a2.reference_start-max_inner_dist+ind
                        if abs(b2_pos - a.reference_end) < 3 and a.reference_name == a2.reference_name:
                            logger.info(f'read1: {a.to_string()}')
                            logger.info(f'read2: {a2.to_string()}')
                            logger.info("以上两条read最后反推回来并没有发现断点")
                        else:
                            if not filter_by_seq_similarity(
                                    masked_ref, a.reference_name, a.reference_end,
                                    a2.reference_name, a2.reference_start
                            ):
                                break1 = a.reference_name + ':' + str(a.reference_end)
                                break2 = a2.reference_name + ':' + str(b2_pos)
                                break_pair = break1 + ' ' + break2
                                sj.setdefault(break_pair, set())
                                sj[break_pair].add(a.seq+a2.seq)
                                fusion_type[break_pair] = determine_fusion_type(a, a2, libtype)

            # (2) 下面假设r1不包含断点，而是r2包含断点，同链融合, 分是否有overlap讨论
            elif a.cigartuples[-1][0] == 0 and a2.cigartuples[0][0] == 4:
                clip_seq = a2.seq[:a.query_alignment_start]
                if len(clip_seq) >= min_clip_size and get_overlap(a.seq, a2.seq, min_clip_size):
                    # 如果 r1和r2有overlap序列，而且overlap不包含断点，但是r2包含断点
                    # 假设overlap至少包含6个碱基，那么clip_seq中前6个碱基一定在r1中
                    ind = a.seq.rfind(reverse_complement(clip_seq[:min_clip_size]))
                    if ind >= 0:
                        b1_pos = a.reference_end - (a.query_length - ind) + len(clip_seq)
                        if abs(b1_pos - a2.reference_start) < 3 and a.reference_name == a2.reference_name:
                            logger.info(f'read1: {a.to_string()}')
                            logger.info(f'read2: {a2.to_string()}')
                            logger.info("以上两条read最后反推回来并没有发现断点")
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

                elif len(clip_seq) >= min_clip_size and not get_overlap(a.seq, a2.seq, 5):
                    # 如果r1和r2没有overlap，但是r2包含断点，且假设r2左端clip掉的seq和r1之间的真实距离在一定范围内
                    # 根据r2左端clip掉的seq一定在r1的下游区域，定位其在下游区域的位置得到的就是第二个断点
                    # print(f"在r1上游的{max_inner_dist}范围内搜索 {clip_seq}")
                    r1_down_seq = ref.fetch(a.reference_name,
                                            a.reference_end,
                                            a.reference_end+max_inner_dist)
                    ind = r1_down_seq.find(reverse_complement(clip_seq))
                    if ind >= 0:
                        b1_pos = a.reference_end + ind
                        if abs(b1_pos - a2.reference_start) < 3 and a.reference_name == a2.reference_name:
                            logger.info(f'read1: {a.to_string()}')
                            logger.info(f'read2: {a2.to_string()}')
                            logger.info("以上两条read最后反推回来并没有发现断点")
                        else:
                            if not filter_by_seq_similarity(
                                    masked_ref, a.reference_name, a.reference_end,
                                    a2.reference_name, a2.reference_start
                            ):
                                break1 = a.reference_name + ':' + str(b1_pos)
                                break2 = a2.reference_name + ':' + str(a2.reference_start)
                                break_pair = break1 + ' ' + break2
                                sj.setdefault(break_pair, set())
                                sj[break_pair].add(a.seq+a2.seq)
                                fusion_type[break_pair] = determine_fusion_type(a, a2, libtype)

    ref.close()
    result = {k:len(v) for k, v in sj.items()}
    fusion_paired = dict(sorted(result.items(), key=lambda x: x[1], reverse=True))

    return fusion_paired, fusion_type


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
                if p[3] == 'keep_up' and p2[3] == 'keep_down':
                    key = p[0] + ':' + str(p[1]) + ' ' + p2[0] + ':' + str(p2[1])
                else:
                    key = p2[0] + ':' + str(p2[1]) + ' ' + p[0] + ':' + str(p[1])
                fusion_dict.setdefault(key, 0)
                fusion_dict[key] += 1
    bam_obj.close()
    return dict(sorted(fusion_dict.items(), key=lambda x: x[1], reverse=True))


def parse_gtf_line(file_obj, sep='\t') -> dict:
    tmp_dict = dict()
    for line in file_obj:
        if line.startswith("#"):
            continue
        tmp_list = line.rstrip().split(sep)
        tmp_dict['chr'] = tmp_list[0]
        tmp_dict['feature'] = tmp_list[2]
        tmp_dict['start'] = tmp_list[3]
        tmp_dict['end'] = tmp_list[4]
        tmp_dict['strand'] = tmp_list[6]
        # parse the column 9
        col_9 = tmp_list[8].strip().split(";")
        for each in col_9[:-1]:
            name = each.split()[0].strip()
            value = each.split()[1].strip().strip('"')
            tmp_dict[name] = value
        yield tmp_dict


def junction_dict_from_gtf(gtf, exon='exon'):
    result = dict()
    with open(gtf) as f:
        for line in parse_gtf_line(f):
            if line['feature'] == 'exon':
                pos = line['chr']+':'+line['start']

@timer
def junction_from_clipped_reads(splits, ref_fasta, masked_ref_fasta=None, min_clip=12, prefix='', max_bp_num=5):
    """good！"""
    # 1. setting logger
    log_name = prefix+'warn.log'
    if os.path.exists(log_name):
        os.remove(log_name)
    logger = set_logger(log_name, logger_id="11")

    log_name = prefix + 'all.junction.points.log'
    if os.path.exists(log_name):
        os.remove(log_name)
    all_break_log = set_logger(log_name, logger_id="22")

    log_name = prefix + 'break_point_pair.tsv'
    if os.path.exists(log_name):
        os.remove(log_name)
    break_pair_txt = set_logger(log_name, logger_id='33')

    log_name = prefix + 'break_nearby.bed'
    if os.path.exists(log_name):
        os.remove(log_name)
    break_nearby_bed = set_logger(log_name, logger_id='55')

    log_name = prefix + 'single_break.log'
    if os.path.exists(log_name):
        os.remove(log_name)
    single_break = set_logger(log_name, logger_id='66')

    # 2. get {junction—around-seq: [support-alignments]}
    jad = junction_around_dict(splits, min_clip=min_clip)

    # 3. 提取断点信息
    single_break_info = dict()
    single_break_points = dict()
    break_pairs = dict()
    ref = pysam.FastaFile(ref_fasta)
    if masked_ref_fasta:
        masked_ref = pysam.FastaFile(masked_ref_fasta)
    else:
        masked_ref = ref
    for jun_seq, reads in jad.items():
        if len(reads) <= 2:
            continue
        break_points = list()
        read_info = list()
        for name, ind, chr_name, point in reads:
            break_points.append((chr_name, point))
            read_info.append((name, ind))
        break_point_set = set(break_points)
        all_break_log.info(jun_seq+":"+str(break_point_set))
        if len(break_point_set) < 2:
            single_break_info[jun_seq] = reads
            bp = list(break_point_set)[0]
            single_break_points[bp] = len(
                {splits[r[0]][r[1]].seq for p, r in zip(break_points, read_info) if p == bp}
            )
            single_break.info(jun_seq+':'+str(break_point_set))
            continue
        elif len(break_point_set) <= max_bp_num:
            if len(break_point_set) > 1:
                logger.info('junction_around_seq: '+jun_seq)
                logger.info('\n'.join(str(x) for x in jad[jun_seq]))

            candidates = dict()
            for p1, p2 in itertools.combinations(break_point_set, 2):
                p1_support_aln = [splits[r[0]][r[1]] for p, r in zip(break_points, read_info) if p == p1]
                p2_support_aln = [splits[r[0]][r[1]] for p, r in zip(break_points, read_info) if p == p2]
                p1_support_seq_count = len({x.seq for x in p1_support_aln})
                p2_support_seq_count = len({x.seq for x in p2_support_aln})
                p1_support_max_aln_len = max(x.query_alignment_length for x in p1_support_aln)
                p2_support_max_aln_len = max(x.query_alignment_length for x in p2_support_aln)
                p1_primary_count = sum(not x.is_secondary for x in p1_support_aln)
                p2_primary_count = sum(not x.is_secondary for x in p2_support_aln)
                supports = [p1_support_seq_count, p2_support_seq_count]
                if (p1_support_seq_count > 1
                        and p2_support_seq_count > 1
                        and p1_primary_count + p2_primary_count > 1
                        and max(supports)/min(supports) < 20
                        and p1_support_max_aln_len > 35
                        and p2_support_max_aln_len > 35
                        and (not filter_by_seq_similarity(masked_ref, p1[0], p1[1], p2[0], p2[1]))
                ):
                    # 方向的确定可能仍然不够精确
                    b1_up_seq = ref.fetch(p1[0], max([p1[1] - min_clip, 0]), p1[1]).upper()
                    b2_up_seq = ref.fetch(p2[0], max([p2[1] - min_clip, 0]), p2[1]).upper()
                    jun_seq_half = jun_seq[:len(jun_seq)//2]
                    score = similarity(None, b1_up_seq, jun_seq_half).ratio()
                    score2 = similarity(None, b2_up_seq, jun_seq_half).ratio()
                    if score > score2 and score > 0.8:
                        fusion = p1[0] + ':' + str(p1[1]) + ' ' + p2[0] + ':' + str(p2[1])
                        candidates[fusion] = (p1_support_seq_count, p2_support_seq_count)
                    else:
                        score = similarity(None, reverse_complement(b1_up_seq), jun_seq_half).ratio()
                        score2 = similarity(None, reverse_complement(b2_up_seq), jun_seq_half).ratio()
                        if score > score2 and score > 0.8:
                            fusion = p1[0] + ':' + str(p1[1]) + ' ' + p2[0] + ':' + str(p2[1])
                            candidates[fusion] = (p1_support_seq_count, p2_support_seq_count)
                        else:
                            fusion = p2[0] + ':' + str(p2[1]) + ' ' + p1[0] + ':' + str(p1[1])
                            candidates[fusion] = (p2_support_seq_count, p1_support_seq_count)
            if candidates:
                fusion, counts = sorted(candidates.items(), key=lambda x: sum(x[1]), reverse=True)[0]
                if fusion not in break_pairs:
                    break_pairs[fusion] = counts
                else:
                    logger.info(f'不同jun_seq:推出了同一个断点{fusion}')
                    break_pairs[fusion] = (break_pairs[fusion][0]+counts[0], break_pairs[fusion][1]+counts[1])
        else:
            logger.info(f'{jun_seq}对应{len(jad[jun_seq])}断点, 而阈值是10, 因此放弃它！具体如下:')
            # logger.info('\n'.join(str(x) for x in jad[jun_seq]))

    # 过滤掉两个断点所支持的reads数目严重不平衡的
    after_filter_pairs = dict()
    for k,v in break_pairs.items():
        ratio = max(v)/min(v)
        if ratio >= 20:
            logger.info(f'{k} was discarded: support ratio={max(v)}/{min(v)}>=20')
            continue
        after_filter_pairs[k] = v

    break_pairs = sorted(after_filter_pairs.items(), key=lambda x: sum(x[1]), reverse=True)

    log_name = 'for_annotate.bed'
    if os.path.exists(log_name):
        os.remove(log_name)
    annotate_bed = set_logger(log_name)
    for fusion, counts in break_pairs:
        chr_name, pos = fusion.split()[0].split(':')
        chr_name2, pos2 = fusion.split()[1].split(':')
        break_pair_txt.info(fusion + '\t' + str(counts[0]) + '\t' + str(counts[1]))
        break_nearby_bed.info(chr_name + '\t' + str(int(pos) - 50) + '\t' + pos)
        break_nearby_bed.info(chr_name2 + '\t' + pos2 + '\t' + str(int(pos2) + 50))
        # 允许3个碱基的误差
        err = 5
        annotate_bed.info(chr_name+'\t'+str(int(pos)-err)+'\t'+str(int(pos)+err)+'\t'+fusion+'\t'+str(counts[0]))
        annotate_bed.info(chr_name2+'\t'+str(int(pos2)-err)+'\t'+str(int(pos2)+err)+'\t'+fusion+'\t'+str(counts[1]))

    ref.close()
    masked_ref.close()
    return break_pairs, single_break_points


def parse_gtf_annotated_bed(file_obj, sep='\t'):
    tmp_dict = dict()
    for line in file_obj:
        if line.startswith("#"):
            continue
        tmp_list = line.rstrip().split(sep)
        tmp_dict['fusion'] = tmp_list[3]
        tmp_list[1] = str(int((int(tmp_list[1]) + int(tmp_list[2]))/2))
        tmp_dict['break_point'] = tmp_list[:2]
        partner = 'break1'
        if tmp_dict['fusion'].split()[1] == ":".join(tmp_list[:2]):
            partner = 'break2'
        tmp_dict['fusion_partner'] = partner
        tmp_dict['support'] = tmp_list[4]
        tmp_dict['chr'] = tmp_list[0+5]
        tmp_dict['feature'] = tmp_list[2+5]
        tmp_dict['start'] = tmp_list[3+5]
        tmp_dict['end'] = tmp_list[4+5]
        tmp_dict['strand'] = tmp_list[6+5]
        # parse the column 9
        col_9 = tmp_list[8+5].strip().split(";")
        for each in col_9[:-1]:
            name = each.split()[0].strip()
            value = each.split()[1].strip().strip('"')
            tmp_dict[name] = value
        yield tmp_dict


@timer
def annotate_break_position(bed, gtf, bam=None, gene_id='gene_id', exon_number='exon_number',
                            gene_name='gene_name', transcript_id='transcript_id',
                            bedtools="bedtools"):
    log_name = 'annotate.log'
    if os.path.exists(log_name):
        os.remove(log_name)
    logger = set_logger(log_name, logger_id='annotate')
    # bedtools intersect
    cmd = f"{bedtools} intersect -loj "
    cmd += f"-a {bed} "
    cmd += f"-b {gtf} "
    cmd += "> {} ".format(bed+'.gtf.annotated')
    check_call(cmd, shell=True)

    annotated_dict = dict()
    # {fusion:{ break1:[support, match_boundary, {gene1:{gene_name, trans:{t1: info, t2:info}}}], break2:[...]}
    with open(bed+'.gtf.annotated', 'r') as fr:
        for line in parse_gtf_annotated_bed(fr):
            chr_name, break_pos = line['break_point']
            break_pos = int(break_pos)
            match_boundary = False
            if chr_name == line['chr']:
                # 判断断点都是已知的外显子边界的
                if abs(break_pos - int(line['start'])) <= 5:
                    match_boundary = chr_name+':'+line['start']
                elif abs(break_pos - int(line['end'])) <= 5:
                    match_boundary = chr_name+':'+line['end']
            fusion_info = annotated_dict.setdefault(line['fusion'], dict())
            detail = fusion_info.setdefault(line['fusion_partner'], [line['support'], dict()])
            if gene_id in line:
                tmp = detail[1].setdefault(line[gene_id], {'match_boundary': False})
                tmp['strand'] = line['strand']
                if (not tmp['match_boundary']) and match_boundary:
                    # 只要有一次匹配则该值永远为真
                    tmp['match_boundary'] = match_boundary
                if gene_name in line:
                    tmp['gene_name'] = line[gene_name]
                if transcript_id in line:
                    trans_info = tmp.setdefault('transcript', dict())
                    trans_info[line[transcript_id]] = line[transcript_id]
                    if line['feature'] == "exon":
                        # 当出现在外显子上，则把intron改成外显子number
                        trans_info[line[transcript_id]] = line[transcript_id]+'|'+line[exon_number]
                    if line['feature'] == "CDS":
                        trans_info[line[transcript_id]] += '|CDS'
                    if line['feature'] == 'UTR':
                        trans_info[line[transcript_id]] += '|UTR'

    with open('annotate.dict', 'w') as f:
        json.dump(annotated_dict, f, indent=2)
    # 有时候断点对应两个基因，无法确定真实对应哪一个基因
    # 如果一个两个断点刚好对应两个已知的junction，则以这两个junction对应的基因为准
    # 如果两个断点都没有对应的基因注释，那么也过滤掉
    discarded_fusion = set()
    if bam:
        bam_obj = pysam.AlignmentFile(bam)
    else:
        bam_obj = None
    paired_gene_dict = dict()
    for fusion in annotated_dict:
        if fusion in discarded_fusion:
            continue
        b1_support, break1_info = annotated_dict[fusion]['break1']
        b2_support, break2_info = annotated_dict[fusion]['break2']
        if not break1_info:  # 没有基因注释, 可以认为是在基因间区
            break1_info['None'] = dict(gene_name=fusion.split()[0], transcript={'None': 'None'})
        if not break2_info:
            break2_info['None'] = dict(gene_name=fusion.split()[1], transcript={'None': 'None'})
        # 如果两个断点可以注释到的基因是一样的, 而且超过两个, 则过滤掉这个融合
        if len(break1_info) >= 2 and not(set(break1_info.keys()) - set(break2_info.keys())):
            discarded_fusion.add(fusion)
            logger.info(f'{fusion} was discarded: 两个断点可以注释到的基因是一样的, 而且超过两个')
            continue
        # 根据注释信息过滤
        for g1, g2 in itertools.product(break1_info.keys(), break2_info.keys()):
            paired_gene_dict.setdefault(g1, set()).add(g2)
            paired_gene_dict.setdefault(g2, set()).add(g1)
            if break1_info[g1]['match_boundary'] and break2_info[g2]['match_boundary'] and g1 == g2:
                # 符合该条件则不考虑其他的基因组合
                discarded_fusion.add(fusion)
                logger.info(f'{fusion} was discarded: 断点均为同一个基因的已知外显子边界')
                break
            elif g1 == g2 == 'None':
                discarded_fusion.add(fusion)
                logger.info(f'{fusion} was discarded: 两个断点均无基因注释')
                break
        if fusion in discarded_fusion:
            continue

        # 根据覆盖度分析和过滤
        if bam:
            # 获得断点1的信息
            chr_name, pos = fusion.split()[0].split(':')
            pos = int(pos)
            break1_around_cov = 0
            break1_around_primary = 0
            break1_around_uniq_seqs = set()
            break1_max_aln_len = 0
            for aln in bam_obj.fetch(chr_name, pos - 50, pos + 3):
                break1_around_cov += 1
                if not aln.is_secondary:
                    break1_around_primary += 1
                break1_around_uniq_seqs.add(aln.seq)
                if aln.query_alignment_length > break1_max_aln_len:
                    break1_max_aln_len = aln.query_alignment_length
            break1_around_uniq_seq_num = len(break1_around_uniq_seqs)

            # 获得断点2的信息
            chr_name, pos = fusion.split()[1].split(':')
            pos = int(pos)
            break2_around_cov = 0
            break2_around_primary = 0
            break2_around_uniq_seqs = set()
            break2_max_aln_len = 0
            for aln in bam_obj.fetch(chr_name, pos - 50, pos + 3):
                break2_around_cov += 1
                if not aln.is_secondary:
                    break2_around_primary += 1
                break2_around_uniq_seqs.add(aln.seq)
                if aln.query_alignment_length > break2_max_aln_len:
                    break2_max_aln_len = aln.query_alignment_length
            break2_around_uniq_seq_num = len(break2_around_uniq_seqs)

            # 整合断点1和断点2的信息进行过滤
            uniq_covs = [break1_around_uniq_seq_num, break2_around_uniq_seq_num]
            if break1_around_primary / break1_around_cov < 0.1:
                discarded_fusion.add(fusion)
                logger.info(f'{fusion} was discarded: break1_around_primary/break1_around_cov < 0.1')
            elif break2_around_primary / break2_around_cov < 0.1:
                discarded_fusion.add(fusion)
                logger.info(f'{fusion} was discarded: break2_around_primary/break2_around_cov < 0.1')
            elif break1_max_aln_len < 50:
                discarded_fusion.add(fusion)
                logger.info(f'{fusion} was discarded: break1_max_aln_len < 50')
            elif break2_max_aln_len < 50:
                discarded_fusion.add(fusion)
                logger.info(f'{fusion} was discarded: break2_max_aln_len < 50')
            elif break1_around_uniq_seq_num - int(b1_support) < 3:
                discarded_fusion.add(fusion)
                logger.info(f'{fusion} was discarded: break1_around_uniq_seq_num-int(b1_support) < 3')
            elif break2_around_uniq_seq_num - int(b2_support) < 3:
                discarded_fusion.add(fusion)
                logger.info(f'{fusion} was discarded: break2_around_uniq_seq_num-int(b2_support) < 3')
            elif break1_around_cov + break2_around_cov <= 10:
                logger.info(f'{fusion} was discarded: break1_around_cov + break2_around_cov <= 10')
                discarded_fusion.add(fusion)
            elif break1_around_cov + break1_around_cov <= 30 and max(uniq_covs) / min(uniq_covs) > 5:
                discarded_fusion.add(fusion)
                logger.info(f'{fusion} was discarded: break1_around_cov+break1_around_cov<=30 and max(covs)/min(covs)>5')
            elif max(uniq_covs) / min(uniq_covs) > 100:
                discarded_fusion.add(fusion)
                logger.info(f'{fusion} was discarded: max(covs)/min(covs)={uniq_covs[0]}/{uniq_covs[1]} > 100')
            elif int(b1_support)/break1_around_uniq_seq_num < 0.0001:
                discarded_fusion.add(fusion)
                logger.info(f'{fusion} was discarded: int(b1_support)/break1_around_uniq_seq_num < 0.0001')
            elif int(b2_support) / break2_around_uniq_seq_num < 0.0001:
                discarded_fusion.add(fusion)
                logger.info(f'{fusion} was discarded: int(b2_support)/break2_around_uniq_seq_num < 0.0001')
            else:
                pass
    if bam:
        bam_obj.close()

    for fusion in list(discarded_fusion):
        discarded_fusion.add(' '.join(fusion.split()[::-1]))

    # write out result
    out_file = os.path.join(os.path.dirname(bed), 'InterGene.fusions.csv')
    out_file2 = os.path.join(os.path.dirname(bed), 'IntraGene.fusions.csv')
    with open(out_file, 'w') as f1, open(out_file2, 'w') as f2:
        header = ['symbol-pair', 'pos-pair', 'b1_support', 'b2_support', 'id-pair', 'b1_annotate', 'b2_annotate']
        f1.write(','.join(header)+'\n')
        f2.write(','.join(header)+'\n')
        for fusion in annotated_dict:
            if fusion in discarded_fusion:
                continue
            b1_support, break1_info = annotated_dict[fusion]['break1']
            b2_support, break2_info = annotated_dict[fusion]['break2']

            if not break1_info:
                break1_info['None'] = dict(gene_name=fusion.split()[0], transcript={'None': 'None'})

            if not break2_info:
                break2_info['None'] = dict(gene_name=fusion.split()[0], transcript={'None': 'None'})

            for g1, g2 in itertools.product(break1_info.keys(), break2_info.keys()):
                if len(paired_gene_dict[g1]) >= 5 or len(paired_gene_dict[g1]) >= 5:
                    logger.info(f'{fusion} was discarded: 一个基因对>=5个融合')
                    continue
                g1d = break1_info[g1]
                g2d = break2_info[g2]
                line = f'{g1d["gene_name"]}--{g2d["gene_name"]},'
                line += f'{fusion},'
                line += f'{b1_support},'
                line += f'{b2_support},'
                line += f'{g1}--{g2},'
                line += f'{";".join(g1d["transcript"].values())},'
                line += f'{";".join(g2d["transcript"].values())}\n'
                if g1 != g2:
                    f1.write(line)
                else:
                    f2.write(line)


# if __name__ == '__main__':
#     from xcmds import xcmds
#     xcmds.xcmds(locals(), include=['get_clipped_reads','get_suspicious_reads', 'single_junction_pos'])
#
#
