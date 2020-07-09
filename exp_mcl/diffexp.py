import numpy as np
import pandas as pd
from scipy.stats import f_oneway


def oneway_anova_test(data, group, out='result.xls'):
    # 单因素方差分析, 设各总体服从正态分布，且方差相同
    data = pd.read_csv(data, header=0, index_col=0, sep=None, engine='python')
    gd = pd.read_csv(group, header=0, index_col=0, sep=None, engine='python')
    gd = gd.iloc[:, [0]]
    gd.columns = ['group']
    data = data.T.join(gd)
    grouped = data.groupby('group')
    stat_func = lambda col:f_oneway(*grouped[col].apply(lambda x:list(x)))[1]
    pvals = map(stat_func, data.columns)
    data.loc['annova_pvalue'] = list(pvals)
    data.T.to_csv(out, sep='\t')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())

