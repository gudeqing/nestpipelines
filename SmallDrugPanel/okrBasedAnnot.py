import re
import pandas as pd
from dataclasses import dataclass, make_dataclass
import pysam
from typing import List
import sqlite3

"""
临床类证据：
1,TRIAL,trial,"Clinical Trials",Trials,,

药监局证据：
2,EMA,label,EMA,EMA,www.ema.europa.eu/ema,"European Medicine Agency"
3,FDA,label,FDA,FDA,www.fda.gov,"United States-Food and Drug Administration"

指南证据：
4,ESMO,guideline,ESMO,ESMO,www.esmo.org,"European Society for Medical Oncology"
5,NCCN,guideline,NCCN,NCCN,www.nccn.org,"United States-National Comprehensive Cancer Network"

注释思路：
目前仅考虑small indel
1.variant 分类，形成marker
2.根据分类marker找到所有三类证据
3.汇总三类证据
4.给出基因的注释信息（narrative.csv）
5.结合实际疾病类型信息，对证据进一步分类，找到推荐治疗指导信息（疾病相关和其他疾病相关）
6.tier分级：(1)指南为一级证据 (2)临床为二级证据 （3）其他证据 （学习dataaccess/GenerateDataAccess.py）

"""


@dataclass
class Variant(object):
    """Class for tracking item in inventory."""
    contig: str
    start: int
    end: int
    ref: str
    alt = str
    chgvs: str = '.'
    phgvs: str = '.'
    vtype: str = 'SNV'
    # marker_id 即variant class，临床证据搜索将依赖该id
    marker_id: str = '.'


def variant_to_marker_id(variant):
    pass


def get_gene_narrative(genes, narrative='narrative.csv'):
    df = pd.read_csv(narrative, index_col=0, header=0)
    df = df.loc[[x in genes for x in df['gene_symbol']]]
    df = df.sort_values(by=['sort_order'])
    pattern = r'\[PMID:.*?\]'
    def get_pmids(section_text):
        pmid_lst = re.findall(pattern, section_text)
        pmid_lst = [x[1: -1].split(';') for x in pmid_lst]
        pmid_lst = [y for x in pmid_lst for y in x]
        return ';'.join(pmid_lst)
    df['PMIDs'] = df['section_text'].apply(get_pmids)
    return df


def marker_to_evidence_id(table='label_evidence_variant_class.csv'):
    """
    该函数适用于3种evidence信息的提取，既可以处理*evidence_variant_class.csv三张表。
    ==> guideline_evidence_variant_class.csv（指南解读） <==
    guideline_evidence_variant_class_id,guideline_evidence_id,variant_class
    1,1,t(12;21)(p13.2;q22.1)
    2,2,"BCR-ABL1 fusion"
    3,3,"ABL1 fusion"
    4,4,"ETV6-RUNX1 fusion"

    ==> label_evidence_variant_class.csv (通常对应的是已知药物信息）<==
    label_evidence_variant_class_id,label_evidence_id,variant_class
    1,1,"ERBB2 overexpression"
    2,2,"ERBB2 amplification"
    3,3,"EGFR S768I mutation"
    4,4,"EGFR L861Q mutation"

    ==> trial_evidence_variant_class.csv （临床解读）<==
    trial_evidence_variant_class_id,open_trial_evidence_id,variant_class
    1,1,"PML-RARA fusion"
    2,2,t(15;17)
    3,3,t(16;16)
    4,4,t(16;16)

    下面以label为为例说明：
    同一个id可能对应不同的marker，如下，但这种情况目前仅发现2个：
    506,529,"RET D898_E901del mutation"
    507,529,"RET D903_S904delinsEP mutation"

    不同的marker可能共享一个id，如下
    477,500,"RET fusion"
    484,507,"RET fusion"
    486,509,"RET fusion"
    487,510,"RET fusion"
    488,511,"RET fusion"
    489,512,"RET fusion"

    :param table:
    :return: a dict likes {'EGFR S768I mutation': {3, 74, 270}, ...}
    """
    data = pd.read_csv(table, index_col=0, header=0)
    marker2label = dict()
    for ind, (label_evidence_id, marker) in data.iterrows():
        marker2label.setdefault(marker, set()).add(label_evidence_id)
    return marker2label


def load_evidences(table='label_evidence.csv'):
    """
    读取具体的evidence信息，该函数可以读取*_evidence.csv三张表:
    ==> guideline_evidence.csv <==
    guideline_evidence_id,guideline_id,evidence_number,marker_type,indication_name,risk_status,gene_symbol,other_criteria,indication_status,therapy_label,simplified_therapy_label,summary,category,population_segment,diagnostic_class
    1,1,19,PROGNOSTIC,"Acute Lymphoblastic Leukemia","ESMO: Favorable","",,INDICATED,,,"",,"",""
    2,1,49,PROGNOSTIC,"Acute Lymphoblastic Leukemia","ESMO: High","",,INDICATED,,,"",,"",""
    3,1,22,PROGNOSTIC,"Acute Lymphoblastic Leukemia","ESMO: High",ABL1,,INDICATED,,,"",,"",""
    4,1,20,PROGNOSTIC,"Acute Lymphoblastic Leukemia","ESMO: Favorable",ETV6,,INDICATED,,,"",,"",""

    ==> label_evidence.csv <==
    label_evidence_id,label_id,evidence_number,indication_name,risk_status,gene_symbol,other_criteria,indication_status,therapy_label,simplified_therapy_label,has_variant_classes
    1,1,1,"Breast Cancer",,ERBB2,,INDICATED,"ado-trastuzumab emtansine","ado-trastuzumab emtansine",t
    2,1,2,"Breast Cancer",,ERBB2,,INDICATED,"ado-trastuzumab emtansine","ado-trastuzumab emtansine",t
    3,2,6,"Non-Small Cell Lung Cancer",,EGFR,,INDICATED,afatinib,afatinib,t
    4,2,5,"Non-Small Cell Lung Cancer",,EGFR,,INDICATED,afatinib,afatinib,t

    ==> open_trial_evidence.csv <==
    open_trial_evidence_id,open_clinical_trial_id,evidence_number,indication_name,risk_status,gene_symbol,other_criteria,exclusion_criteria_flag,regimen,therapy_label
    1,1,8,"Acute Promyelocytic Leukemia",,RARA,,f,"1, 2","supplement, chemotherapy"
    2,1,7,"Acute Promyelocytic Leukemia",,RARA,,f,"1, 2","supplement, chemotherapy"
    3,2,6,"Acute Myeloid Leukemia",,CBFB,,f,"1, 2, 3","gemtuzumab ozogamicin, cytokine, chemotherapy"
    4,2,6,"Myelodysplastic Syndrome",,CBFB,,f,"1, 2, 3","gemtuzumab ozogamicin, cytokine, chemotherapy"

    以label_evidence为例，不同的label_id可能对应多个label_evidence_id, 他们对应不同的用药说明如FDA，EMA
    当“other_criteria”列不为空时，需同时满足“other_criteria”中的变异条件，才能使用该条目中的用药指导信息。
    :param table:
    :param indication: 适应症，如"Breast Cancer"
    :return:返回以label_evidence_id为key的字典，类似：
    {1: {'label_id': 1,
     'evidence_number': 1,
     'indication_name': 'Breast Cancer',
     'risk_status': nan,
     'gene_symbol': 'ERBB2',
     'other_criteria': nan,
     'indication_status': 'INDICATED',
     'therapy_label': 'ado-trastuzumab emtansine',
     'simplified_therapy_label': 'ado-trastuzumab emtansine',
     'has_variant_classes': 't'}, ...}
    """
    data = pd.read_csv(table, index_col=0, header=0)
    return data.to_dict(orient='index')


def evidence_source_detail(table='label.csv'):
    """
    用于获取每个证据的来源信息，如根据guidelines_id从guideline.txt中查询当前evidence的来源
    :param table:
    :return: label id as key, dict likes:
    {1: {'label_type': 'EMA',
     'content_name': 'EMA-ado-trastuzumab emtansine',
     'therapy_label': 'ado-trastuzumab emtansine',
     'country': 'European Union',
     'indication_and_usage': nan,
     'label_url': 'https://www.ema.europa.eu/en/documents/product-information/kadcyla-epar-product-information_en.pdf',
     'label_approved_date': '2020-01-20'}, ...}
    """
    data = pd.read_csv(table, index_col=0, header=0)
    counts = data.iloc[:, 0].value_counts().to_dict()
    print(f'{data.columns[0]} stat: {counts}')
    return data.to_dict(orient='index')


def open_clinical_trial_id_to_local_segment(table='trial_location.csv'):
    """
    该函数用来读表trial_location.csv trial_population_segment.csv
    理出open_clinical_trial_id和研究地点及人群信息的关系
    :param table:
    :return:
    """
    data = pd.read_csv(table, index_col=0, header=0)
    trial2location = dict()
    for ind, (open_clinical_trial_id, location) in data.iterrows():
        trial2location.setdefault(open_clinical_trial_id, set()).add(location)
    return trial2location


def location_detail(table='location.csv'):
    data = pd.read_csv(table, index_col=0, header=0)
    return data.to_dict(orient='index')


def therapy_summary(markers:tuple=("EGFR L861Q mutation", )):
    label_evidences_dict = marker_to_evidence_id('label_evidence_variant_class.csv')
    guideline_evidences_dict = marker_to_evidence_id('guideline_evidence_variant_class.csv')
    trial_evidences_dict = marker_to_evidence_id('trial_evidence_variant_class.csv')

    print(f'there are {len(label_evidences_dict)} markers associated with label evidences.')
    print(f'there are {len(guideline_evidences_dict)} markers associated with guideline evidences.')
    print(f'there are {len(trial_evidences_dict)} markers associated with clinical trial evidences.')

    label_evidence_detail = load_evidences('label_evidence.csv')
    guideline_evidence_detail = load_evidences('guideline_evidence.csv')
    trial_evidence_detail = load_evidences('open_trial_evidence.csv')

    for marker in markers:
        label_evidences = label_evidences_dict[marker]
        guideline_evidences = guideline_evidences_dict[marker]
        trial_evidences = trial_evidences_dict[marker]

        # get detail
        label_detail = [label_evidence_detail[x] for x in label_evidences]
        guideline_detail = [guideline_evidence_detail[x] for x in guideline_evidences]
        trial_detail = [trial_evidence_detail[x] for x in trial_evidences]

        # print
        print(marker, ':')
        label_df = pd.DataFrame(label_detail)
        guideline_df = pd.DataFrame(guideline_detail)
        trial_df = pd.DataFrame(trial_detail)
        label_df.to_csv('E1.csv')
        guideline_df.to_csv('E2.csv')
        trial_df.to_csv('E3.csv')
        gene_df = get_gene_narrative(gene, )



if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
