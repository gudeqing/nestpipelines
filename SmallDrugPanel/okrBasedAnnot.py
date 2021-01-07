import pandas as pd
from dataclasses import dataclass, make_dataclass
from typing import List
import sqlite3
# 目前仅考虑small indel
print(f'there are {len(marker2label)} markers associated with guidelines')
print(f'there are {data.shape[0]} label evidences for indication {indication}')

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
1.variant 分类，形成marker
2.根据分类marker找到所有三类证据
3.汇总三类证据
4.给出基因的注释信息
5.结合实际疾病类型信息，对证据进一步分类，找到推荐治疗指导信息
6.tier分级：(1)指南为一级证据 (2)临床为二级证据 （3）其他证据

"""


@dataclass
class Variant(object):
    """Class for keeping track of an item in inventory."""
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


def variant_to_marker_id():
    pass


def get_gene_narrative(table=''):
    pass


def marker_to_evidence_id(table='label_evidence_variant_class.csv'):
    """
    该函数适用于3种evidence信息的提取，下面以label为为例说明

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


def load_evidences(table='label_evidence.csv', indication=None):
    """
    不同的label_id可能对应多个label_evidence_id, 他们对应不同的用药说明如FDA，EMA
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
    if indication:
        data = data[data['indication_name'] == indication]
    return data.to_dict(orient='index')


def evidence_source_detail(table='label.csv'):
    """
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


def therapy_summary(vcf):
    pass


def report(summary):
    pass