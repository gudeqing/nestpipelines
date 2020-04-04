import matplotlib
matplotlib.use('agg')

import pandas as pd
import statsmodels.api as sm
import pylab as pl
import numpy as np
from sklearn.metrics import roc_curve, auc
import sys


def single_lgr(data, y_col='y',  x_cols:list=None, drop_cols:list=None, prefix='Result'):
    df = pd.read_csv(data, header=0, sep=None, engine='python')
    # Dependent and Independent Variables
    y = df[y_col].copy()
    X = df.drop(columns= y_col)
    if drop_cols:
        for each in drop_cols:
    	    X = X.drop(columns=each)
    if x_cols:
        X = df[x_cols]


    data = X
    train_cols = data.columns

    # manually add the intercept
    data['intercept'] = 1.0
    data[y_col] = y
    res_data = list()
    for col in train_cols:
        print('>>>Analysis variable', col)
        # fit the model
        input_data = data[[y_col, col, 'intercept']].dropna()
        y_data = input_data[y_col]
        train_data = input_data[[col, 'intercept']]
        result = sm.Logit(y_data, train_data).fit()
        print(result.summary())
        model = result
        stat = model.conf_int()
        stat.columns = ['conf_lower', 'conf_upper']
        stat['pvalues'] = model.pvalues
        stat['coef'] = model.params
        target = stat.loc[[col]]
        intercept = stat.loc['intercept', 'coef']
        predict_data = result.predict(train_data)
        fpr, tpr, thresholds =roc_curve(y_data, predict_data)
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
        print(roc)
        
        # Plot tpr vs 1-fpr
        fig, ax = pl.subplots()
        pl.plot(roc['tpr'])
        pl.plot(roc['1-fpr'], color = 'red')
        pl.xlabel('1-False Positive Rate')
        pl.ylabel('True Positive Rate')
        pl.title('Receiver operating characteristic')
        ax.set_xticklabels([])
        pl.savefig(f'{col}.ROC.png')
        pl.close()

    pd.concat(res_data).to_csv(f'{prefix}.xls', sep='\t')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['single_lgr'])

