import gzip
targets = set(x.strip() for x in open('target.genes'))
target_mutations = set(x.strip() for x in open('target.mutation.list'))
# AKT1:NM_001014431:exon3:c.49G>A:p.E17K
mut_c_dict = dict()
mut_p_dict = dict()
for x in target_mutations:
    g,t,e,c,p = x.split(':')
    mut_c_dict.setdefault(g, set())
    mut_p_dict.setdefault(g, set())
    mut_c_dict[g].add(c)
    mut_p_dict[g].add(p)

with gzip.open('/nfs2/database/COSMIC/V91/CosmicCodingMuts.normal.vcf.gz', 'rt') as f, \
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
            new_line = 'chr'+line.rstrip('\n')+';'+'is_target=0'+'\n'
            if gene in mut_c_dict:
                if info['CDS'] in mut_c_dict[gene] \
                        or info['HGVSC'] in mut_c_dict[gene]\
                        or info['AA'] in mut_p_dict[gene] \
                        or info['HGVSP'] in mut_p_dict[gene]:
                    new_line = 'chr'+line.rstrip('\n')+';'+'is_target=1'+'\n'
                    fw2.write(new_line)
            fw.write(new_line)



