import os
import gzip


def main(target_genes='target.genes', target_mutations='target.mutations',
         cosmic='/nfs2/database/COSMIC/V91/CosmicCodingMuts.normal.vcf.gz'):
    """
    根据cosmic的数据CosmicCodingMuts.normal.vcf.gz提取目标基因的突变信息
    由于hgvs格式的不确定性，提取结果不一定正确，请自己核查
    :param target_genes: 目标基因文件，每行一个gene symbol
    :param target_mutations: 目标突变，每行一个突变，格式如 AKT1:NM_001014431:exon3:c.49G>A:p.E17K
    :param cosmic: cosmic的突变文件
    :return:
    """
    targets = set(x.strip() for x in open(target_genes))

    mut_c_dict = dict()
    mut_p_dict = dict()

    if os.path.exists(target_mutations):
        target_mutations = set(x.strip() for x in open(target_mutations))
    else:
        target_mutations = set()

    for x in target_mutations:
        g,t,e,c,p = x.split(':')
        mut_c_dict.setdefault(g, set())
        mut_p_dict.setdefault(g, set())
        mut_c_dict[g].add(c)
        mut_p_dict[g].add(p)

    with gzip.open(cosmic, 'rt') as f, \
            open('target.genes.vcf', 'w') as fw, \
            open('target.genes.target_mutation.vcf', 'w') as fw2:
        for line in f:
            if line.startswith('#'):
                fw.write(line)
                fw2.write(line)
                continue

            lst = line.strip().split('\t')
            info = dict(x.split('=') for x in lst[7].split(';') if len(x.split('='))==2)
            if 'CDS' not in info:
                info['CDS'] = ''
            if 'HGVSC' not in info:
                info['HGVSC'] = ''
            if 'AA' not in info:
                info['AA'] = ''
            if 'HGVSP' not in info:
                info['HGVSP'] = ''
            gene = info['GENE']
            if gene in targets:
                fw.write(line)
                if gene in mut_c_dict:
                    if info['CDS'] in mut_c_dict[gene] \
                            or info['HGVSC'] in mut_c_dict[gene]\
                            or info['AA'] in mut_p_dict[gene] \
                            or info['HGVSP'] in mut_p_dict[gene]:
                        fw2.write(line)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['main'])
