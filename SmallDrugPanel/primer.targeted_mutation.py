# bedtools intersect -a primer_targeted_250bp.bed -b ../TargetMutations/target_genes.cosmic_coding_mutation.in_ROI.vcf -wao > primer.target_mutation.txt


def primer_target_mutation_stat(file, mut_format=('GENE', 'CDS', 'AA')):
    """
    主要用来统计出每个gene primer对应哪些已知突变信息，针对从cosmic提取的突变信息
    :param file: 使用bedtools intersect产生的文件，如:
    bedtools intersect -a primer_targeted_250bp.bed -b ../TargetMutations/target_genes.cosmic_coding_mutation.in_ROI.vcf -wao > primer.target_mutation.txt
    :param mut_format: 表示突变的主要字段信息
    :return:
    """
    result = dict()
    # cosmic = dict()
    with open(file) as f, open('primer_mutations.detail.txt', 'w') as fw:
        for line in f:
            lst = line.strip().split('\t')
            target_region = lst[0]+':'+lst[1]+'-'+lst[2]
            primer = lst[3]+'|'+target_region
            cosmic_id = lst[8]
            mutation_info = lst[13]
            fw.write(primer+'\t'+mutation_info+'\n')
            result.setdefault(primer, list())
            # cosmic.setdefault(primer, list())
            # cosmic[primer].append(cosmic_id)
            if len(mutation_info) <= 1:
                if len(cosmic_id) > 1:
                    result[primer].append(mutation_info+':'+cosmic_id)
                else:
                    result[primer].append(mutation_info)
            else:
                try:
                    info_dict = dict(x.split('=', 1) for x in mutation_info.split(';') if len(x.split('=', 1))>=2)
                except:
                    print(mutation_info)
                simple_mut = [info_dict['is_target']]
                for each in mut_format:
                    if each in info_dict:
                        simple_mut.append(info_dict[each])
                simple_mut.append(cosmic_id)
                mut = ':'.join(simple_mut)
                if mut not in result[primer]:
                    result[primer].append(mut)

    with open('primer_mutations.simple.txt', 'w') as f:
        f.write('primer\ttarget_region\tcosmic_mutations\ttarget_mutation\n')
        for k, v in result.items():
            muts = '|'.join(sorted(v, key=lambda x:x.split(':')[0], reverse=True))
            primer, target_region = k.split('|')
            target_muts = '|'.join([x for x in v if x.startswith('1:')])
            if not target_muts:
                target_muts = '.'
            f.write(f'{primer}\t{target_region}\t{muts}\t{target_muts}\n')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())

