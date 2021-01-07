import pysam


def add_contig_header(vcf, ref='hg19'):
    contig_info = [
        f'##assembly={ref}',
        "##contig=<ID=chr1,length=249250621>",
        "##contig=<ID=chr2,length=243199373>",
        "##contig=<ID=chr3,length=198022430>",
        "##contig=<ID=chr4,length=191154276>",
        "##contig=<ID=chr5,length=180915260>",
        "##contig=<ID=chr6,length=171115067>",
        "##contig=<ID=chr7,length=159138663>",
        "##contig=<ID=chrX,length=155270560>",
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
        "##contig=<ID=chr20,length=63025520>",
        "##contig=<ID=chr19,length=59128983>",
        "##contig=<ID=chr22,length=51304566>",
        "##contig=<ID=chr21,length=48129895>",
        "##contig=<ID=chrMT,length=16569>",
        "##contig=<ID=chrY,length=59373566>",
        "##contig=<ID=chrX,length=156040895>",
    ]
    for line in contig_info:
        vcf.header.add_line(line)


def filter_vcf(vcf, out):
    vcf = pysam.VariantFile(vcf)
    if vcf.header.contigs.__len__() < 1:
        add_contig_header(vcf)
    vcf_out = pysam.VariantFile(out, "w", header=vcf.header)
    samples = list(vcf.header.samples)
    for r in vcf:
        info = r.info
        af = r.samples[samples[-1]]['AF'][0]
        f0 = af >= 0.04
        f1 = True if info['esp6500siv2_all'] is None else float(info['esp6500siv2_all']) <= 0.01
        f2 = True if info['1000g2015aug_all'] is None else float(info['1000g2015aug_all']) <= 0.01
        f3 = True if info['ExAC_ALL'] is None else float(info['ExAC_ALL']) <= 0.01
        f4 = True if info['gnomAD_exome_ALL'] is None else float(info['gnomAD_exome_ALL']) <= 0.01
        f5 = True if info['gnomAD_exome_EAS'] is None else float(info['gnomAD_exome_EAS']) <= 0.01
        f6 = True if info['ExAC_EAS'] is None else float(info['ExAC_EAS']) <= 0.01
        for j, d in zip(
                [f0, f1, f2, f3, f4, f5, f6],
                ['af', 'esp6500siv2_all', '1000g2015aug_all', 'ExAC_ALL', 'gnomAD_exome_ALL', 'gnomAD_exome_EAS', 'ExAC_EAS']
        ):
            if not j:
                # print(samples[-1], r.pos, r.info['Gene_refGene'], d)
                pass
        # f7 = a['SBF'] >= 0.05
        if all([f0, f1, f2, f3, f4, f5, f6]):
            vcf_out.write(r)
    else:
        vcf_out.close()

def filter_nirvan_json(js):
    import json
    import gzip
    # a['positions'][0]['variants'][0]['transcripts'][0]['hgnc']


def reorder_vcf(vcf, out):
    vcf = pysam.VariantFile(vcf)
    vcf_out = pysam.VariantFile(out, "w", header=vcf.header)
    newlines = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype: 0/1 if AF < 0.9 else 1/1 ">',
        '##FORMAT=<ID=AD,Number=2,Type=Integer,Description="(ref depth, Variant Depth)">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        # '##FORMAT=<ID=VF,Number=A,Type=Float,Description="Allele Frequency">'
    ]
    for line in newlines:
        vcf_out.header.add_line(line)

    for r in vcf:
        record = vcf_out.new_record()
        record.contig = r.contig
        record.start = r.start
        record.ref = r.ref
        record.alts = r.alts
        record.qual = r.qual
        for each in r.filter:
            record.filter.add(each)
        record.info.update(r.info)
        sample = list(r.samples)[0]
        record.samples[sample]['GT'] = r.samples[sample]['GT']
        record.samples[sample]['AD'] = (r.samples[sample]['DP']-r.samples[sample]['AO'][0], r.samples[sample]['AO'][0])
        record.samples[sample]['DP'] = r.samples[sample]['DP']
        record.samples[sample]['AF'] = r.samples[sample]['AF']
        vcf_out.write(record)
    vcf.close()
    vcf_out.close()


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())


