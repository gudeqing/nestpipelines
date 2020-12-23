from dataclasses import dataclass, make_dataclass
from typing import List
import sqlite3


class Variant(object):
    # 构建一个容器存储突变信息，可以从vcf的pysam.vcf.record开始
    pass

    def variant_to_variant_class(self):
        # 把variant转换为variant_class，方便后续兵分三路找治疗证据
        pass

def get_gene_narrative():
    pass


def variant_class_to_label():
    pass


def variant_class_to_guideline():
    pass


def variant_class_to_trail():
    pass


def therapy_summary(vcf):
    pass


def report(summary):
    pass