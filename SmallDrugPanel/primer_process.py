import pysam
genome = '/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta'
gn = pysam.FastaFile(genome)
gtf = pysam.TabixFile('gencode.v19.annotation.gff3.gz', parser=pysam.asGFF3())


def reverse_complement(seq):
    """
    :param seq: 输入序列
    :return: 返回reversed的序列
    """
    seq = seq.upper()
    complement = dict(zip(list("ATCG"), list("TAGC")))
    return ''.join(complement[base] if base in complement else base for base in seq[::-1])


with open('primer.txt') as f, open('primer_for_cutadapt.fasta', 'w') as fw, open('primer.bed', 'w') as fw2:
    for line in f:
        lst = line.strip().split()
        chr_, pos, length = lst[0], int(lst[1]), len(lst[3])
        # guess coordinate
        start = pos - length + 1
        end = pos + 1
        ref = gn.fetch(chr_, start, end).upper()
        if ref == lst[3]:
            strand = '+'
        else:
            # guess again
            start = pos
            end = pos + length
            ref = gn.fetch(chr_, start, end).upper()
            if reverse_complement(ref) == lst[3]:
                ref = reverse_complement(ref)
                strand = '-'
            else:
                raise Exception(f'参考基因组和当前序列信息不一致: {line}')
        # get gene annotation
        info = dict()
        for each in gtf.fetch(chr_, start, end):
            gene_name = each['gene_name']
            gene_id = each['gene_id']
            trans = each['transcript_id']
            direct = list(each)[6]
            info.setdefault(gene_id, [])
            if direct not in info[gene_id]:
                info[gene_id].append(direct)
            if gene_name not in info[gene_id]:
                info[gene_id].append(gene_name)
            if trans not in info[gene_id]:
                info[gene_id].append(trans)
        info = str(info).replace(' ', '')
        # to fasta file
        # name = f'{chr_}:{start + 1}-{end}:len={length}:strand={strand}:annot={info}'
        # name = f'{chr_}:{start+1}-{end}:len={length}:strand={strand}:gene={gene_name}'
        name = f'{chr_}:{start+1}-{end}:{length}:{strand}:{gene_name}'
        fw.write(f'>{name}\n^{ref}\n')
        # to bed file
        fw2.write(f'{chr_}\t{start}\t{end}\t{info}\t.\t{strand}\n')

gn.close()
gtf.close()
