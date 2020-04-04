import math
from math import log10
import numpy as np
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import pandas as pd


def get_18gene_score_weights(table=None):
    if table:
        wt = pd.read_csv(table, header=None, sep=None, engine='python')
        wt_pair = zip(wt.iloc[:, 0], wt.iloc[:, 1])
        return dict(wt_pair)

    weights = {
        'CCL5': 0.008346,
        'CD27': 0.072293,
        'CD274': 0.042853,
        'CD276': -0.0239,
        'CD8A': 0.031021,
        'CMKLR1': 0.151253,
        'CXCL9': 0.074135,
        'CXCR6': 0.004313,
        # 'HLA.DQA1': 0.020091,
        'HLA-DQA1': 0.020091,
        # 'HLA.DRB1': 0.058806,
        'HLA-DRB1': 0.058806,
        # 'HLA.E': 0.07175,
        'HLA-E': 0.07175,
        'IDO1': 0.060679,
        'LAG3': 0.123895,
        'NKG7': 0.075524,
        'PDCD1LG2': 0.003734,
        'PSMB10': 0.032999,
        'STAT1': 0.250229,
        'TIGIT': 0.084767
    }
    return weights


def get_house_keeping_genes():
    return [
        'STK11IP',
        'ZBTB34',
        'TBC1D10B',
        'OAZ1',
        'POLR2A',
        'G6PD',
        'ABCF1',
        # 'C14orf102', replace with alia NRDE2
        'NRDE2',
        'UBB',
        'TBP',
        'SDHA'
    ]


def scoring(exp_table, out=None):
    table = pd.read_csv(exp_table, header=0, index_col=0, sep=None, engine='python')
    hkg_include = []
    hkg_missing = []
    for x in get_house_keeping_genes():
        if x in table.index:
            hkg_include.append(x)
        else:
            hkg_missing.append(x)
    if hkg_missing:
        print(f'These house keeping genes were not found: {hkg_missing}')

    predictor_include = []
    predictor_missing = []
    weight_dict = get_18gene_score_weights()
    for each in weight_dict:
        if each in table.index:
            predictor_include.append(each)
        else:
            predictor_missing.append(each)
    if predictor_missing:
        print(f'These predictor genes were not found: {predictor_missing}')

    hkg_exp = table.loc[hkg_include]
    hkg_exp = hkg_exp.applymap(lambda x: np.log10(x))
    hkg_exp_mean = np.sum(hkg_exp, axis=0) / hkg_exp.shape[0]

    predictor_exp = table.loc[predictor_include]
    norm_exp = predictor_exp.applymap(lambda x: np.log10(x)).apply(lambda x: x - hkg_exp_mean, axis=1)
    w = [weight_dict[x] for x in predictor_exp.index]
    score = norm_exp.apply(lambda x: sum(x * w), axis=0)
    score.to_csv(out or 'gep_score.csv', header=False)


def compare_score(score, score2):
    parameter = np.polyfit(score, score2, 1)  # 1次幂拟合
    func = np.poly1d(parameter)  # 方程
    # plt.plot(score, score, 'g--')
    plt.plot(score, func(score), 'r--')
    plt.scatter(score, score2)
    plt.annotate(f'y={round(func[1], 3)}x{round(func[0], 3)}', xy=(0, 0.5))
    plt.xlabel("old score")
    plt.ylabel('new score')
    plt.savefig('score.png', dpi=300)
    # plt.grid(True)
    plt.close()


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['scoring'])
