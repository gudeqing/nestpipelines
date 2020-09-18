import pysam


def annot_bed(bed, out, gff='gencode.v19.annotation.gff3.gz', transcript=None, trans_col=1, level='exon', strip_version=False):
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
    # fetch one or more rows in a region using 0-based indexing.
    # The region is specified by reference, start and end.
    # Alternatively, a samtools region string can be supplied.
    # 每一个record的内容虽然可以用字典的方式访问，但是其本身不是字典，不可以使用in来判断key是否存在
    with open(bed) as fr, open(out, 'w') as fw:
        for line in fr:
            lst = line.strip().split('\t')
            contig, start, end = lst[0], lst[1], lst[2]
            info = []
            gene_info = []
            for record in gff.fetch(contig, int(start)-1, int(end)):
                gene_name = record['gene_name']
                gene_id = record['gene_id']
                trans_id = record['transcript_id']
                if strip_version:
                    gene_id = gene_id.rsplit('.', 1)[0]
                    trans_id = trans_id.rsplit('.', 1)[0]
                strand = record['strand']
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
                if record['feature'] == 'gene':
                    gene_info.append([gene_name, strand, gene_id])

            if not info:
                info = gene_info
                if not info:
                    # print('NoInformation for', line)
                    info.append([f'No{level.capitalize()}OrGeneInformation'])
                else:
                    info[-1].append(f'NotIn{level.capitalize()}')
                if len(lst) >=4:
                    info[-1].append(lst[3])

            if transcripts:
                new_info = []
                for each in info:
                    if len(each) >= 4:
                        if each[3] in transcripts:
                            new_info.append(each)
                if new_info:
                    # 如果没有常用转录本，则使用注释到的所有转录本信息
                    info = new_info
            info = ';'.join(':'.join(x) for x in info)
            new_line = [contig, start, end, info]
            fw.write('\t'.join(new_line)+'\n')


def stat_gene_and_bed(gene_bed, cov_bed, exon_bed):
    cov_bed = pysam.TabixFile(cov_bed, parser=pysam.asBed())
    exon_bed = pysam.TabixFile(exon_bed, parser=pysam.asBed())
    result = dict()
    with open(gene_bed) as f:
        for line in f:
            chr_, start, end, gene = line.strip().split()
            regions = cov_bed.fetch(chr_, start, end)
            for region in regions:
                pass


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

