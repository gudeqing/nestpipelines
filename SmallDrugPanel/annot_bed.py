import pysam


def annot_bed(bed, out, gff='gencode.v19.annotation.gff3.gz', transcript=None, trans_col=1, level='exon', strip_version=False):
    """
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
    transcripts = {x.strip().split()[trans_col].rsplit('.', 1)[0] for x in open(transcript)}
    # fetch one or more rows in a region using 0-based indexing.
    # The region is specified by reference, start and end.
    # Alternatively, a samtools region string can be supplied.
    # 每一个record的内容虽然可以用字典的方式访问，但是其本身不是字典，不可以使用in来判断key是否存在
    with open(bed) as fr, open(out, 'w') as fw:
        for line in fr:
            lst = line.strip().split('\t')
            contig, start, end = lst[0], lst[1], lst[2]
            info = []
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

            if transcripts:
                new_info = []
                for each in info:
                    if each[3] in transcripts:
                        new_info.append(each)
                if not new_info:
                    # 如果没有常用转录本，则使用注释到的所有转录本信息
                    new_info = info
            info = new_info
            if not info:
                info = [[f'No{level.capitalize()}Information']]
                if len(lst) >=4:
                    info[0][0] = info[0][0]+ ';'+lst[3]
            info = ';'.join(':'.join(x) for x in info)
            new_line = [contig, start, end, info]
            fw.write('\t'.join(new_line)+'\n')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())

