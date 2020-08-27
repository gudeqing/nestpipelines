import gzip
import re
import pandas as pd
import numpy as np

"""
该脚本是针对panel设计初步阶段进行的开发的，主要是为了获取：cosmic和FusionGDB的融合，PharmGKB的drug位点信息
"""

def get_cosmic_fusion(infile='target.fusion_genes.list', break_extend=100, out_off_target_break=False,
                      database='fusionDatabase/CosmicFusionExport.tsv.gz'):
    cols = [
        'Sample ID',
        'Sample name',
        'Primary site',
        'Site subtype 1',
        'Site subtype 2',
        'Site subtype 3',
        'Primary histology',
        'Histology subtype 1',
        'Histology subtype 2',
        'Histology subtype 3',
        'Fusion ID',
        'Translocation Name',
        "5'_CHROMOSOME",
        "5'_GENOME_START_FROM",
        "5'_GENOME_START_TO",
        "5'_GENOME_STOP_FROM",
        "5'_GENOME_STOP_TO",
        "5'_STRAND",
        "3'_CHROMOSOME",
        "3'_GENOME_START_FROM",
        "3'_GENOME_START_TO",
        "3'_GENOME_STOP_FROM",
        "3'_GENOME_STOP_TO",
        "3'_STRAND",
        'Fusion type',
        'Pubmed_PMID'
    ]

    fusion_id_tuple = [
        'Fusion ID',
        'Translocation Name',
        "5'_CHROMOSOME",
        "5'_GENOME_START_FROM",
        "5'_GENOME_START_TO",
        "5'_GENOME_STOP_FROM",
        "5'_GENOME_STOP_TO",
        "5'_STRAND",
        "3'_CHROMOSOME",
        "3'_GENOME_START_FROM",
        "3'_GENOME_START_TO",
        "3'_GENOME_STOP_FROM",
        "3'_GENOME_STOP_TO",
        "3'_STRAND",
    ]

    # target_genes = {x.strip() for x in open('gene.list')}
    # print(target_genes)

    target_fusion_genes = {x.strip() for x in open(infile)}
    print(target_fusion_genes)
    data = pd.read_table(database)
    def to_int(x):
        try:
            return int(x)
        except:
            return x
    data = data.fillna('nan')
    data = data.applymap(to_int)
    primary_sites = set(data['Primary site'])
    print('primary_sites number:', len(primary_sites))
    result = dict()
    break_lst = []
    genes = set()
    fusions = []
    for index, row in data.iterrows():
        # fid = row['Fusion ID']
        fid = [row[x] for x in fusion_id_tuple]
        # fid = '|'.join([str(x) for x in fid])
        name = str(row['Translocation Name'])
        if not fid or (not name):
            continue
        symbols = [x.strip('(').strip(')') for x in re.findall(r'\(.*?\)', name)]
        # if len(set(symbols) - target_genes) == 0 and len(set(symbols) & target_fusion_genes)>=1:
        if len(set(symbols) & target_fusion_genes)>=1:
            genes |= set(symbols) & target_fusion_genes
            symbols = [x for x in symbols if not x.isnumeric()]
            if len(symbols) > 2:
                symbol_lst = []
                for each in symbols:
                    if each not in symbol_lst:
                        symbol_lst.append(each)
                symbols = symbol_lst
                if len(symbols) > 2:
                    print('found more than 2 genes in fusion', name)
            symbol_pair = '--'.join(symbols)
            fid.append(symbol_pair)
            chr_name = 'chr' + str(row["5'_CHROMOSOME"])
            chr_name2 = 'chr' + str(row["3'_CHROMOSOME"])
            if row["5'_STRAND"] == '+':
                break1 = row["5'_GENOME_STOP_TO"]
                break1_str = chr_name + ':' + str(break1)
            else:
                break1 = row["5'_GENOME_START_FROM"]
                break1_str =  chr_name + ':' + str(break1)

            if row["3'_STRAND"] == '+':
                break2 = row["3'_GENOME_START_FROM"]
                break2_str = chr_name2 + ':' + str(break2)
            else:
                break2 = row["3'_GENOME_STOP_FROM"]
                break2_str = chr_name2 + ':' + str(break2)

            fid.append(break1_str)
            fid.append(break2_str)
            fid = tuple(fid)
            strand5 = row["5'_STRAND"]
            strand3 = row["3'_STRAND"]
            mutation = f'{symbols[0]}:{strand5}:{chr_name}:{break1}--{symbols[1]}:{strand3}:{chr_name2}:{break2}'
            # to bed
            if type(break1) == int and type(break2) == int:
                if out_off_target_break:
                    break_point = [chr_name, break1-break_extend-1, break1+break_extend, mutation]
                    if break_point not in break_lst:
                        break_lst.append(break_point)
                    break_point2 = [chr_name2, break2-break_extend-1, break2+break_extend, mutation]
                    if break_point2 not in break_lst:
                        break_lst.append(break_point2)
                else:
                    if symbols[0] in target_fusion_genes:
                        break_point = [chr_name, break1-break_extend-1, break1+break_extend, mutation]
                        if break_point not in break_lst:
                            break_lst.append(break_point)
                    if symbols[1] in target_fusion_genes:
                        break_point2 = [chr_name2, break2-break_extend-1, break2+break_extend, mutation]
                        if break_point2 not in break_lst:
                            break_lst.append(break_point2)

                bedpe = [
                    chr_name, str(break1-1), str(break1),
                    chr_name2, str(break2-1), str(break2),
                    f'{symbols[0]}--{symbols[1]}', '.',
                    strand5, strand3, '.'
                ]
                if bedpe not in fusions:
                    fusions.append(bedpe)

            # print('--'.join(symbols))
            result.setdefault(fid, dict())
            result[fid].setdefault(row['Primary site'], set())
            result[fid][row['Primary site']].add(row['Sample ID'])
    if out_off_target_break:
        out_name = f'CosmicFusion.target.breaks.{break_extend}_extended.bed'
    else:
        out_name = f'CosmicFusion.OnTarget.breaks.{break_extend}_extended.bed'
    with open(out_name, 'w') as f:
        for each in break_lst:
            f.write('\t'.join([str(x) for x in each])+'\n')

    # 计算样本个数
    for fusion, site_dict in result.items():
        for k, v in site_dict.items():
            result[fusion][k] = len(v)

    not_found = target_fusion_genes - genes
    if not_found:
        print('No fusion found for genes:', not_found)
    else:
        print('found fusion for all target genes!')

    result = pd.DataFrame(result)
    result = result.transpose().fillna(0)
    result.index.names = fusion_id_tuple + ['SymbolPair', 'Break1', 'Break2']
    result.to_csv('CosmicFusion.target.txt', sep='\t')
    result.to_excel('CosmicFusion.target.xlsx', merge_cells=False)
    # print(all_fusion_pairs)
    with open('CosmicFusion.target.bedpe', 'w') as f:
        for each in fusions:
            f.write('\t'.join(each)+'\n')
    return fusions


def get_PharmGKB_mutation(infile='target_drug_genes.list',
                          drug_ann='drugDatabase/clinical_ann_metadata.tsv',
                          dbsnp='/nfs2/database/annovar/humandb/hg19_avsnp150.txt'):
    targets = set([x.strip() for x in open(infile)])
    dbsnp_ids = dict()
    # Clinical Annotation ID  Location        Gene    Level of Evidence
    # Clinical Annotation Types       Genotype-Phenotype IDs
    # Annotation Text Variant Annotations IDs Variant Annotations
    # PMIDs   Evidence Count  Related Chemicals       Related Diseases
    # Biogeographical Groups  Chromosome
    with open(drug_ann) as f, open('target.PharmGKB.annot.txt', 'w') as fw:
        fw.write(f.readline())
        for line in f:
            id_, loc, gene, *_ = line.split('\t', 4)
            if gene.strip():
                gene = gene.split()[0]
                if gene in targets:
                    fw.write(line)
                    if loc.startswith('rs'):
                        dbsnp_ids[loc] = gene

    with open(dbsnp) as f, open('target.dbsnp.mutation.txt', 'w') as fw:
        # 1       10019   10020   TA      T       rs775809821
        for line in f:
            rs_id = line.strip().rsplit('\t', 1)[1]
            if rs_id in dbsnp_ids:
                fw.write('chr'+line[:-1]+'\t'+dbsnp_ids[rs_id]+'\n')


def get_FusionGDB_fusion(infile='target.fusion_genes.list', break_extend=100, out_off_target_break=False,
                         database='fusionDatabase/TCGA_ChiTaRS_combined_fusion_information_on_hg19.txt'):
    targets = set([x.strip() for x in open(infile)])
    break_lst = []
    fusions = []
    genes = set()
    with open(database) as f, open('FusionGDB.target.txt', 'w') as fw:
        # Data-Source-Ctype-Sample-Hgene-Hchr-Hbp-Hstrand-Tgene-Tchr-Tbp-Tstrand
        # TCGA	LD	UVM	TCGA-VD-A8KM-01A	TMEM259	chr19	1020771	-	PRSS57	chr19	694970	-
        fw.write('Data-Source-Ctype-Sample-Hgene-Hchr-Hbp-Hstrand-Tgene-Tchr-Tbp-Tstrand\n'.replace('-', '\t'))
        for line in f:
            lst = line.rstrip().split('\t')
            Data, Source, Ctype, Sample, Hgene, Hchr, Hbp, Hstrand, Tgene, Tchr, Tbp, Tstrand = lst
            if Hgene in targets or Tgene in targets:
                genes |= targets & {Hgene, Tgene}
                fw.write(line)
                mutation = f'{Hgene}:{Hstrand}:{Hchr}:{Hbp}--{Tgene}:{Tstrand}:{Tchr}:{Tbp}'
                break1 = f'{Hchr}\t{int(Hbp) - break_extend-1}\t{int(Hbp) + break_extend}\t{mutation}\n'
                break2 = f'{Tchr}\t{int(Tbp) - break_extend-1}\t{int(Tbp) + break_extend}\t{mutation}\n'
                if not out_off_target_break:
                    if break1 not in break_lst and (Hgene in targets):
                        break_lst.append(break1)
                    if break2 not in break_lst and (Tgene in targets):
                        break_lst.append(break2)
                else:
                    if break1 not in break_lst:
                        break_lst.append(break1)
                    if break2 not in break_lst:
                        break_lst.append(break2)

                bedpe = [
                    Hchr, str(int(Hbp) - 1), Hbp,
                    Tchr, str(int(Tbp) - 1), Tbp,
                    f'{Hgene}--{Tgene}', '.',
                    Hstrand, Tstrand, '.'
                ]
                if bedpe not in fusions:
                    fusions.append(bedpe)
    not_found = targets - genes
    if not_found:
        print('No fusion found for genes:', not_found)
    else:
        print('found fusion for all target genes!')

    if out_off_target_break:
        out_name = f'FusionGDB.target.breaks.{break_extend}_extended.bed'
    else:
        out_name = f'FusionGDB.OnTarget.breaks.{break_extend}_extended.bed'
    with open(out_name, 'w') as fw:
        for each in break_lst:
            fw.write(each)
    with open('FusionGDB.target.bedpe', 'w') as f:
        for each in fusions:
            f.write('\t'.join(each)+'\n')
    return fusions


def get_hotspot(infile='target.gene.list', hots='simplified.hotspot.txt', out='target.gene.hotspot.txt'):
    targets = {x.strip() for x in open(infile)}
    found = set()
    with open(hots) as fr, open(out, 'w') as fw:
        fw.write(fr.readline())
        for line in fr:
            lst = line.strip().split('\t')
            if lst[9] in targets or lst[11] in targets:
                found.add(lst[9])
                found.add(lst[11])
                fw.write(line)
    not_found = targets - found
    if found:
        print("Found no hotspot for genes:", not_found)
    else:
        pass


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
