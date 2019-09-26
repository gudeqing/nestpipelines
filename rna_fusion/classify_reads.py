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

spanned_between_genes = dict()
spanned_within_genes = dict()
splitted = dict()
possible_splitted = dict()


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


def single_junction_pos(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    junction = dict()
    for n, r in enumerate(bam):
        if not r.is_paired:
            continue
        cigar = r.cigarstring
        if (not cigar) or ('H' in cigar or 'P' in cigar):
            continue
        cigar_tuples = r.cigartuples
        if cigar_tuples[0][0] == 4 and cigar_tuples[0][1] > 10:
            # print(cigar)
            left_break = r.reference_start
            key = r.reference_name + ':' + str(left_break)
            junction.setdefault(key, 0)
            junction[key] += 1
        if cigar_tuples[-1][0] == 4 and cigar_tuples[-1][1] > 10:
            right_break = r.reference_end
            key = r.reference_name + ':' + str(right_break)
            junction.setdefault(key, 0)
            junction[key] += 1
    else:
        bam.close()

    # sort
    print(dict(sorted(junction.items(), key=lambda x: x[1], reverse=True)[:50]))
    print(len(junction))
    if 'chr1:154170400' in junction:
        print(junction['chr1:154170400'])
    return dict(sorted(junction.items(), key=lambda x: x[1], reverse=True))


def split_reads(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    splits = dict()
    for n, r in enumerate(bam):
        if not r.is_paired:
            continue
        cigar = r.cigarstring
        if (not cigar) or ('H' in cigar or 'P' in cigar):
            continue
        cigar_tuples = r.cigartuples
        left_clip = cigar_tuples[0][0] == 4 and cigar_tuples[0][1] > 6
        right_clip = cigar_tuples[-1][0] == 4 and cigar_tuples[-1][1] > 6
        if left_clip or right_clip:
            record = splits.setdefault(r.query_name + '_' + str(int(r.is_read1)), {})
            record[r.get_tag('HI')] = r
    else:
        bam.close()
    return splits


def paired_splits(bam_file, splits):
    bam = pysam.AlignmentFile(bam_file, "rb")
    target_reads = dict()
    split_names = set(x[:-2] for x in splits.keys())
    for n, r in enumerate(bam):
        if not r.is_paired:
            continue
        if r.query_name in split_names:
            record = target_reads.setdefault(r.query_name + '_' + str(int(r.is_read1)), {})
            record[r.get_tag('HI')] = r
    bam.close()
    return target_reads


def junction_from_one_end(splits):
    sj = set()
    for name, aligns in splits.items():
        # 一条read可能有多个比对结果, 排列组合查看
        for a, a2 in itertools.permutations(aligns.values(), 2):
            if abs(a.reference_start - a2.reference_start) > 200 \
                    or a.reference_name != a2.reference_name:
                if a.cigartuples[-1][0] == 4 and a2.cigartuples[0][0] == 4:
                    m = a.seq[a.query_alignment_end:] == a2.seq[a2.query_alignment_start:]
                    m2 = a.seq[:a.query_alignment_end] == a2.seq[:a2.query_alignment_start]
                    if m and m2:
                        # print(a.cigarstring, a2.cigarstring)
                        break1 = a.reference_name + ':' + str(a.reference_end)
                        break2 = a2.reference_name + ':' + str(a2.reference_start)
                        if abs(a2.reference_start - a.reference_end) > 50 or\
                                (a.reference_name != a2.reference_name):
                            sj.add(break1+'  '+break2)
    print(len(sj))
    sj = sorted(list(sj))
    return sj


def junction_from_pair_end(paired_splits, ref_fasta, max_inner_dist=300):
    ref = pysam.FastaFile(ref_fasta)
    sj = set()
    names = [x[:-1] for x in paired_splits.keys() if x.endswith('0')]
    for name in names:
        r1_name, r2_name = name + '0', name + '1'
        if (r1_name not in paired_splits) or (r2_name not in paired_splits):
            continue
        r1 = paired_splits.pop(r1_name)
        r2 = paired_splits.pop(r2_name)
        for a, a2 in itertools.product(r1.values(), r2.values()):
            distance = abs(a.reference_start - a2.reference_start)
            if not a.cigartuples or (not a2.cigartuples):
                continue
            if distance > 200 or (a.reference_name != a.reference_name):
                # 如果r1和r2有overlap序列，而且overlap包含断点
                if a.cigartuples[-1][0] == a2.cigartuples[0][0] == 4:
                    evidence1 = a.seq[:a.query_alignment_end].endswith(a2.seq[:a2.query_alignment_start])
                    evidence2 = a2.seq[a2.query_alignment_start:].startswith(a.seq[a.query_alignment_end:])
                    if evidence1 and evidence2:
                        break1 = a.reference_name + ':' + str(a.reference_end)
                        break2 = a2.reference_name + ':' + str(a2.reference_start)
                        sj.add(break1+'  '+break2)

                # 仅r1包含断点，而r2不包含断点
                elif a.cigartuples[-1][0] == 4 and a2.cigartuples[0][0] == 0:
                    clip_seq = a.seq[a.query_alignment_end:]
                    # 如果 r1和r2有overlap序列，而且overlap不包含断点，但是r1包含断点
                    # 假设overlap至少包含6个碱基，那么r2前6个碱基一定包含在r1右端被clip掉的reads中
                    if len(clip_seq) >= 6:
                        # 假设是r1和r2存在overlap,但overlap不完全包含断点
                        ind = clip_seq.rfind(a2.seq[:6])
                        if ind >= 0:
                            b2_pos = a2.reference_start - ind
                            if abs(b2_pos - a.reference_end) < 3 and a.reference_name == a2.reference_name:
                                print("有时会由于多了几个错配而发生clip, 最后反推回来并没有发现断点")
                            else:
                                break1 = a.reference_name + ':' + str(a.reference_end)
                                break2 = a2.reference_name + ':' + str(b2_pos)
                                sj.add(break1+'  '+break2)

                    if len(clip_seq) >= 8:
                        # 如果r1和r2没有overlap，但是r1包含断点，且假设r1右端clip掉的seq和r2之间的真实距离在一定范围内
                        # 根据r1右端clip掉的seq一定在r2的上游区域，定位其在上游区域的位置得到的就是第二个断点
                        # print(f"在r2上游的{max_inner_dist}范围内搜索 {clip_seq}")
                        r2_up_seq = ref.fetch(a2.reference_name,
                                              a2.reference_start-max_inner_dist,
                                              a2.reference_start)
                        ind = r2_up_seq.rfind(clip_seq)
                        if ind >= 0:
                            b2_pos = a2.reference_start-max_inner_dist+ind
                            if abs(b2_pos - a.reference_end) < 3 and a.reference_name == a2.reference_name:
                                print("有时会由于多了几个错配而发生clip, 最后反推回来并没有发现断点")
                            else:
                                break1 = a.reference_name + ':' + str(a.reference_end)
                                break2 = a2.reference_name + ':' + str(b2_pos)
                                sj.add(break1+'  '+break2)

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
                                print("有时会由于多了几个错配而发生clip, 最后反推回来并没有发现断点")
                            else:
                                break1 = a.reference_name + ':' + str(b1_pos)
                                break2 = a2.reference_name + ':' + str(a2.reference_start)
                                sj.add(break1+'  '+break2)

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
                                print("有时会由于多了几个错配而发生clip, 最后反推回来并没有发现断点")
                            else:
                                break1 = a.reference_name + ':' + str(b1_pos)
                                break2 = a2.reference_name + ':' + str(a2.reference_start)
                                sj.add(break1+'  '+break2)
    ref.close()
    return sj


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['get_clipped_reads','get_suspicious_reads', 'single_junction_pos'])


