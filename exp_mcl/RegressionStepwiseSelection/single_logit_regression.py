import matplotlib
matplotlib.use('agg')

import pandas as pd
import statsmodels.api as sm
from matplotlib import pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve, auc
from bokeh.plotting import figure, save, output_file
from bokeh.models import ColumnDataSource, CDSView, GroupFilter, HoverTool, LabelSet, Legend
from bokeh.layouts import gridplot
import sys


def single_lgr(data, y_col='y',  x_cols:list=None, drop_cols:list=None, factorize:tuple=None, prefix='Result'):
    """

    :param data:
    :param y_col:
    :param x_cols:
    :param drop_cols:
    :param prefix:
    :param factorize:
    :return:
    """
    df = pd.read_csv(data, header=0, sep=None, index_col=0, engine='python')
    # Dependent and Independent Variables
    if factorize is not None:
        fac = {x:int(y) for x, y in zip(factorize[::2], factorize[1:][::2])}
        y_data = [fac[x] for x in df[y_col]]
    else:
        if df[y_col].dtype != int:
            y_data, rule = df[y_col].factorize()
            print(dict(zip(rule, range(len(rule)))))
        else:
            y_data = df[y_col]

    data = df.drop(columns=y_col)
    if drop_cols:
        for each in drop_cols:
            data = data.drop(columns=each)
    if x_cols:
        data = df[x_cols]
    target_cols = data.columns

    data[y_col] = y_data
    # manually add the intercept
    data['intercept'] = 1.0

    # # 一起做回归分析
    # final_data = data.dropna()
    # unfit_model = sm.Logit(final_data[y_col], final_data.drop(columns=y_col))
    # model = unfit_model.fit()
    # print(model.summary())

    # 逐一进行回归分析
    res_data = list()
    row_num = len(target_cols)//3
    if row_num == 0:
        row_num = 1
    elif len(target_cols)%3 != 0:
        row_num += 1
    plot_options = dict(
        width=250,
        plot_height=250,
        tools='pan,wheel_zoom,box_select, reset,save',
    )
    plots = list()
    for col in target_cols:
        print('>>>Analysis variable', col)
        # fit the model
        my_data = data[[y_col, col, 'intercept']].dropna()
        y_data = my_data[y_col]
        train_data = my_data[[col, 'intercept']]
        print(y_data)
        unfit_model = sm.Logit(y_data, train_data)
        model = unfit_model.fit()
        print(model.summary())
        stat = model.conf_int()
        stat.columns = ['conf_lower', 'conf_upper']
        stat['pvalues'] = model.pvalues
        stat['coef'] = model.params
        target = stat.loc[[col]]
        intercept = stat.loc['intercept', 'coef']
        predict_data = model.predict(train_data)
        fpr, tpr, thresholds = roc_curve(y_data, predict_data)
        roc_auc = auc(fpr, tpr)
        train_data[y_col] = y_data
        train_data['intercept'] = intercept
        train_data['predict'] = predict_data

        ####################################
        # The optimal cut off would be where tpr is high and fpr is low
        # tpr - (1-fpr) is zero or near to zero is the optimal cut off point
        ####################################
        optimal_idx = np.argmax(tpr - fpr)
        optimal_threshold = thresholds[optimal_idx]
        # prob = 1/(1+np.exp(-(b*x+c)))
        cutoff = (np.log(1/optimal_threshold - 1)*(-1) - intercept)/target['coef']
        target['cutoff'] = cutoff
        print('# predict result')
        print(train_data)
        print("Area under the ROC curve : %f" % roc_auc)
        print('Optimal probability threshold is', optimal_threshold)
        print(f'Conclude that optimal variable cutoff of', cutoff)
        res_data.append(target[['coef', 'cutoff', 'pvalues', 'conf_lower', 'conf_upper']])
        i = np.arange(len(tpr)) # index for df
        roc = pd.DataFrame({'fpr' : pd.Series(fpr, index=i),'tpr' : pd.Series(tpr, index = i), '1-fpr' : pd.Series(1-fpr, index = i), 'tf' : pd.Series(tpr - (1-fpr), index = i), 'thresholds' : pd.Series(thresholds, index = i)})
        # optimal_threshold = roc.iloc[(roc.tf-0).abs().argsort()[:1]]['thresholds']
        # print(roc)

        # Plot tpr vs 1-fpr
        # plt.plot(roc['tpr'], color='green', label='TPR')
        # plt.plot(roc['1-fpr'], color = 'red', label='1-FPR')
        s = figure(**plot_options)
        s.line(roc['fpr'], roc['tpr'], color='blue', legend_label=f'{col}  AUC={roc_auc:.2f}')
        s.xaxis.axis_label = 'FPR'
        s.yaxis.axis_label = 'TPR'
        s.legend.location = 'bottom_right'
        plots.append((s, roc_auc))

    plots = [x[0] for x in sorted(plots, key=lambda x:x[1], reverse=True)]
    p = gridplot(plots, sizing_mode='stretch_{}'.format('width'), ncols=3)
    output_file(f'{prefix}.ROC.html', title="ROC")
    save(p)
    pd.concat(res_data).sort_values(by='pvalues').to_csv(f'{prefix}.xls', sep='\t')




if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['single_lgr'])

