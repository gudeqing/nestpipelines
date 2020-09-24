import pysam


def is_contained(a, b):
    """
    a      s1---e1
    b   s2----------e2
    a s1---------------e1
    :param a:(s1, e1)
    :param b: (s2, e2)
    :return:
    """
    overlap = 0
    s_dif = a[0]-b[0]
    e_dif = a[1]-b[1]
    if s_dif*e_dif <= 0:
        if s_dif > 0:
            overlap = a[1] - a[0]
        elif s_dif == 0 and e_dif < 0:
            overlap = a[1] - a[0]
        else:
            overlap = b[1] - b[0]
    return overlap


def is_part_overlap(a, b):
    """
    不包含完全包含的关系
    s2-----e2
        s1-----e1
            s2-------e2
    :param a:
    :param b:
    :return: overlap
    """
    overlap = 0
    if (b[0] < a[0]) and (b[1] < a[1]) and (a[0] <= b[1]):
        overlap = b[1] - a[0]
    elif (b[0] > a[0]) and (b[1] > a[1]) and (a[1] >= b[0]):
        overlap = a[1] - b[0]
    return overlap


def annot_bed(bed, out, gff='gencode.v19.annotation.gff3.gz', transcript=None,
              trans_col=1, level='exon', strip_version=False, prior_genes=None):
    """
    使用gff文件对bed文件进行注释
    info order: [gene_name, strand, gene_id, trans_id, exon_number, phase]
    :param bed: 至少三列
    :param out: 输出文件名
    :param gff: gff3文件路径，目前仅测试了gencode的gff3
    :param transcript: 如果提供该文件，该文件某一列包含转录本信息
    :param trans_col: 指示transcript文件哪一列包含转录本信息，默认第二列
    :param level: one of ['gene', 'transcript', 'exon'(default), 'CDS']
     for CDS:
     A phase of "0" indicates that a codon begins on the first nucleotide of the CDS feature (i.e. 0 bases forward),
     a phase of "1" indicates that the codon begins at the second nucleotide of this CDS feature and
     a phase of "2" indicates that the codon begins at the third nucleotide of this region.
    :param strip_version: 如果设置该参数，则剔除基因或转录本的version信息
    :return:
    """
    features = ['gene', 'transcript', 'exon', 'CDS']
    if level not in features:
        raise Exception('level must be one of', features)
    gff = pysam.TabixFile(gff, parser=pysam.asGFF3())
    if transcript:
        transcripts = {x.strip().split()[trans_col].rsplit('.', 1)[0] for x in open(transcript)}
    else:
        transcripts = set()
    if prior_genes:
        prior_genes = {x.strip() for x in open(prior_genes)}
    else:
        prior_genes = set()
    # fetch one or more rows in a region using 0-based indexing.
    # The region is specified by reference, start and end.
    # Alternatively, a samtools region string can be supplied.
    # 每一个record的内容虽然可以用字典的方式访问，但是其本身不是字典，不可以使用in来判断key是否存在
    hit_genes = set()
    with open(bed) as fr, open(out, 'w') as fw:
        for line in fr:
            lst = line.strip().split('\t')
            contig, start, end = lst[0], int(lst[1]), int(lst[2])
            info = []
            gene_info = []
            gene_to_discard = set()
            for record in gff.fetch(contig, start, end):
                gene_name = record['gene_name']
                gene_id = record['gene_id']
                strand = record['strand']
                trans_id = record['transcript_id']
                if strip_version:
                    gene_id = gene_id.rsplit('.', 1)[0]
                    trans_id = trans_id.rsplit('.', 1)[0]

                # 单独提取基因信息，当区域注释不到其他指定的feature时，则使用这个基因注释信息
                if record['feature'] == 'gene':
                    if 'gene_type' in record.keys() and \
                            (record['gene_type'] == 'antisense'
                             or record['gene_type'] == 'processed_transcript'
                             or 'pseudogene' in record['gene_type']
                             or record['gene_type'] == 'lincRNA'
                             or record['gene_type'] == 'sense_intronic'
                             or record['gene_type'] == 'sense_overlapping'
                            ):
                        # 根据注释到的基因类型考虑是否添加该基因注释
                        if gene_name not in prior_genes:
                            gene_to_discard.add(gene_name)
                    else:
                        # 计算基因feature与当前区域的交集情况，根据交集情况进行过滤
                        a = sorted([start, end])
                        b = sorted([record['start'], record['end']])
                        ovelap = is_contained(a, b) or is_part_overlap(a, b)
                        if ovelap > (a[1] - a[0]) * 0.9:
                            gene_info.append([gene_name, strand, gene_id])
                        else:
                            # 要求查询区域的90%包含于gene区域, 否则丢弃该基因相关注释
                            gene_to_discard.add(gene_name)

                if gene_name in gene_to_discard:
                    continue

                if record['feature'] == level:
                    if level == 'exon':
                        exon_number = record['exon_number']
                        info.append([gene_name, strand, gene_id, trans_id, str(exon_number)])
                    elif level == 'gene':
                        info.append([gene_name, strand, gene_id])
                    elif level == 'transcript':
                        info.append([gene_name, strand, gene_id, trans_id])
                    elif level == 'CDS':
                        exon_number = record['exon_number']
                        phase = record['phase']
                        info.append([gene_name, strand, gene_id, trans_id, str(exon_number), str(phase)])

            if gene_to_discard:
                info = [x for x in info if x[0] not in gene_to_discard]
                print(f'基因{gene_to_discard}和查询的{start, end}只有<=90%的交集或其为antisense或lincRNA或假基因，故丢弃其相关注释')

            # 如果还有多个基因的注释，仅保留优先基因集中的基因的注释，如果不在优先集，则保留所有
            includes = [x for x in info if x[0] in prior_genes]
            if includes:
                discards = {x[0] for x in info if x[0] not in prior_genes} # 不在优先集的基因
                # print(discards, '不在优先集')
                gene_to_discard |= discards # 更新要提取的基因
                gene_info = [x for x in gene_info if x[0] not in prior_genes]
                info = includes

            if len(set([x[0] for x in gene_info])) > 1:
                # print(info)
                print(f'{lst[:3]}在排除反义基因等等过滤后仍然可以注释到多个基因，它们是{gene_info}')

            if not info:
                # 没有指定的feature注释
                info = gene_info
                if not info:
                    # 也没有基因注释
                    info.append(['NoGeneInfo'])
                    if len(lst) >= 4:
                        # 当 原始bed文件有4列时，在没有注释的情况下，使用原有的注释
                        info[-1].append(lst[3])
                else:
                    # 有基因注释，但是没有指定feature注释
                    if level != 'gene':
                        info[-1].append(f'No{level.capitalize()}Info')
                    else:
                        raise Exception('info没有基因注释，但是gene_info却有注释，矛盾！')

            if transcripts:
                new_info = []
                discard_genes = set()
                for each in info:
                    if len(each) >= 4:
                        # 如果有转录本注释，则看转录本是否为经典转录本
                        if each[3] in transcripts:
                            new_info.append(each)
                        else:
                            discard_genes.add(each[0])
                if new_info:
                    # 如果有经典转录本注释，则只使用提取出来的经典转录本注释
                    info = new_info

            hit_genes |= set([x[0] for x in info if not x[0].startswith('No')])
            info = ';'.join(':'.join(x) for x in info)
            new_line = [contig, str(start), str(end), info]
            fw.write('\t'.join(new_line)+'\n')

    with open('hit.gene.list', 'w') as f:
        for name in hit_genes:
            f.write(name+'\n')


def get_canonical_transcript_exon(gff, canonical_trans, out):
    trans = {x.strip().split()[1].rsplit('.', 1)[0] for x in open(canonical_trans)}
    with open(out, 'w') as fw, open(gff) as fr:
        for line in fr:
            if line.startswith('#'):
                continue
            lst = line.strip().split()
            if lst[2] == 'exon':
                col9 = dict(x.split('=', 1) for x in lst[8].split(';'))
                if col9['transcript_id'].rsplit('.', 1)[0] in trans:
                    fw.write(line)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())

