import itertools

""""
"""


def parse_gtf_line(file_obj, sep='\t') -> dict:
    tmp_dict = dict()
    for line in file_obj:
        if line.startswith("#"):
            continue
        tmp_list = line.rstrip().split(sep)
        tmp_dict['chr'] = tmp_list[0]
        tmp_dict['feature'] = tmp_list[2]
        tmp_dict['start'] = tmp_list[3]
        tmp_dict['end'] = tmp_list[4]
        tmp_dict['strand'] = tmp_list[6]
        # parse the column 9
        col_9 = tmp_list[8].strip().split(";")
        for each in col_9[:-1]:
            name = each.split()[0].strip()
            value = each.split()[1].strip().strip('"')
            tmp_dict[name] = value
        yield tmp_dict


def annotate_dict_from_gtf(gtf, gene_id='gene_id', exon_id='exon_number',
                            gene_name='gene_name', transcript_id='transcript_id'):
    """
    先梳理出每个转录本的外显子信息，再把外显子的起始坐标排序，最后都到外显子和外显子之间的junction坐标
    exon1: s1 e1
    exon2: s2 e2
    exon3: s3 e3
    => junction (e1, s2), (e2,s3)
    {'tid': {
        gene_id: gene_id,
        gene_name: gene_name,
        strand: '+',
        exon: { exon1: (s1, e1), exon2:(s2, e2), exon3:(s3, e3)}
        }
    }
    :param gtf:
    :param exon_id:
    :return:
    """
    trans_info = dict()
    with open(gtf) as f:
        for line in parse_gtf_line(f):
            if line['feature'] == 'exon':
                info = trans_info.setdefault(line[transcript_id], dict())
                info['chr'] = line['chr']
                info['gene_id'] = line[gene_id]
                info['gene_name'] = line[gene_name]
                info['strand'] = line['strand']
                exon_dict = info.setdefault('exon', dict())
                exon_dict[line[exon_id]] = (line['start'], line['end'])

    junctions = dict()
    exon_boundary = dict()
    for tid in trans_info:
        exon_coords = trans_info[tid]['exon'].values()
        exon_coords = sorted([int(x) for y in exon_coords for x in y])
        gene_info = dict(
            gene_id=trans_info[tid]['gene_id'],
            gene_name=trans_info[tid]['gene_name'],
            strand=trans_info[tid]['strand'],
            transcript_id=tid,
        )
        # junction pair as key
        for i in range(len(exon_coords)-1):
            if i % 2 == 1:
                key = '-'.join([trans_info[tid]['chr'], str(exon_coords[i]), str(exon_coords[i+1])])
                if key not in junctions:
                    junctions[key] = gene_info
                else:
                    if trans_info[tid]['gene_id'] == junctions[key]['gene_id']:
                        if tid not in junctions[key][transcript_id].split('|'):
                            junctions[key][transcript_id] += '|' + tid
                    else:
                        print(f"发现极端情况: 一个Junction {key}属于不同的基因"
                              f"{junctions[key]['gene_id']}和{trans_info[tid]['gene_id']}")

        # exon boundary as key
        for coord in exon_coords:
            key = '-'.join([trans_info[tid]['chr'], str(coord)])
            if key not in exon_boundary:
                exon_boundary[key] = gene_info
            else:
                if trans_info[tid]['gene_id'] == exon_boundary[key]['gene_id']:
                    if tid not in exon_boundary[key][transcript_id].split('|'):
                        exon_boundary[key][transcript_id] += '|' + tid
                else:
                    print(f"发现极端情况: 一个Exon boundary {key}属于不同的基因"
                          f"{exon_boundary[key]['gene_id']}和{trans_info[tid]['gene_id']}")

    return junctions, exon_boundary


def extract_novel_mapsplice_junction(junction_file, gtf, miss=3, min_cov=2, min_uniq_read_count=1):
    """
    依据GTF文件从mapsplice2软件的输出结果过滤出新的junction位点，并提供基因注释。
    novel junction位点定义：两个位点必须是已知外显子的边界
    :param junction_file:
    :param gtf:
    :param miss:
    :param min_cov:
    :param min_uniq_read_count:
    :return:
    """
    known_junctions, exon_boundary = annotate_dict_from_gtf(gtf)
    with open(junction_file) as fr,\
            open('known.junction.txt', 'w') as fwk,\
            open('novel.junction.txt', 'w') as fwn:
        for line in fr:
            lst = line.split()
            chr_name, donor, acceptor, name, coverage, strand = lst[:6]
            unique_read = lst[20]
            multiple_read = lst[21]
            if int(coverage) >= min_cov and int(unique_read) >= min_uniq_read_count:
                sites = sorted([int(donor), int(acceptor)])
                site1_range = range(sites[0]-miss, sites[0]+miss+1)
                site2_range = range(sites[1]-miss, sites[1]+miss+1)
                is_known_junction = False
                is_new_junction = False
                key = None
                key1, key2 = None, None
                for s1, s2 in itertools.product(site1_range, site2_range):
                    key = '-'.join([chr_name, str(s1), str(s2)])
                    if key in known_junctions:
                        is_known_junction = True
                        break
                    key1 = '-'.join([chr_name, str(s1)])
                    key2 = '-'.join([chr_name, str(s2)])
                    if (key1 in exon_boundary) and (key2 in exon_boundary):
                        # 要求两个位点都是已知外显子边界或边界附近, 且要求其在是同一个基因上
                        if exon_boundary[key1]['gene_id'] == exon_boundary[key2]['gene_id']:
                            is_new_junction = True
                            break

                line_info = lst[:6] + [unique_read, multiple_read]
                if is_known_junction:
                    line_info += list(known_junctions[key].values())
                    fwk.write('\t'.join(line_info) + '\n')
                elif is_new_junction:
                    info = exon_boundary[key1].copy()
                    info['transcript_id'] += '|' + exon_boundary[key2]['transcript_id']
                    line_info += list(info.values())
                    fwn.write('\t'.join(line_info) + '\n')
                else:
                    pass


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['extract_novel_mapsplice_junction'])

