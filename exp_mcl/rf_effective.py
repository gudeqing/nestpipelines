import matplotlib
matplotlib.use('agg')
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from numpy import random
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn import decomposition, preprocessing
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.metrics import roc_curve, auc
import pickle


def rf_result(exp_matrix, group_info, standard_scale=False, target_rows=None, target_cols=None,
              train_samples=None, n_times=1000):
    data = read_data(exp_matrix, group_info, standard_scale=standard_scale,
                     target_rows=target_rows, target_cols=target_cols)
    if train_samples is None:
        x_train, x_test, y_train, y_test = train_test_split(data.iloc[:, :-1], data.iloc[:, -1], test_size=0.2)
    else:
        train_sample_list = [x.strip().split()[0] for x in open(train_samples)]
        x_train = data.loc[train_sample_list].iloc[:, :-1]
        y_train = data.loc[train_sample_list].iloc[:, -1]
        x_test = data.loc[[x for x in data.index if x not in train_sample_list]].iloc[:, :-1]
        y_test = data.loc[[x for x in data.index if x not in train_sample_list]].iloc[:, -1]

    test_scores = list()
    estimator_list = list()
    f_importance = pd.DataFrame()
    for i in range(n_times):
        estimator = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
                                           max_depth=None, max_features=16, max_leaf_nodes=None,
                                           min_impurity_decrease=0.0, min_impurity_split=None,
                                           min_samples_leaf=1, min_samples_split=2,
                                           min_weight_fraction_leaf=0.0, n_estimators=110, n_jobs=1,
                                           oob_score=False, random_state=None, verbose=0,
                                           warm_start=False)
        estimator.fit(x_train, y_train)
        estimator_list.append(estimator)
        y_predict = estimator.predict(x_test)
        # print(classification_report(y_test, y_predict))
        test_scores.append(accuracy_score(y_test, y_predict))
        f_importance[i] = estimator.feature_importances_
    else:
        f_importance.index = x_train.columns
        feature_score = f_importance.transpose().describe().transpose().sort_values(by=['mean', 'std'], ascending=False)
        print(feature_score.head(50))
        feature_score.to_csv('feature_score.csv')
        test_score_pd = pd.Series(test_scores).describe()
        print(test_score_pd)
        test_score_pd.to_csv('test_data_scores.csv')

    # plot
    models = [estimator_list[test_scores.index(max(test_scores))], estimator_list[test_scores.index(test_score_pd['50%'])]]
    plt.title('Receiver Operating Characteristic')
    for name, model in zip(['max', 'median'], models):
        pickle.dump(model, open(f"{name}_test_score.model", 'wb'))
        probs = model.predict_proba(x_test)
        # print(probs)
        preds = probs[:, 1]
        fpr, tpr, threshold = roc_curve(y_test, preds)
        roc_auc = auc(fpr, tpr)
        color = 'r' if name == 'max' else 'b'
        plt.plot(fpr, tpr, color, label=f'model_with_{name}_test_score: AUC = %0.2f' % roc_auc)
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], 'g--')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.savefig('ROC.png', dpi=300)
    plt.close()


def read_data(exp_matrix, group_info, standard_scale=False, target_rows=None, target_cols=None):
    data = pd.read_csv(exp_matrix, header=0, index_col=0, sep=None, engine='python')
    if target_rows is not None:
        data = data.loc[[x.strip().split()[0] for x in open(target_rows)]]
    if target_cols is not None:
        data = data[[x.strip().split()[0] for x in open(target_cols)]]
    # print(((data.std(axis=1) / data.mean(axis=1)).abs()).describe())
    # data = data[(data.std(axis=1) / data.mean(axis=1)).abs() > 0.01]
    data = data.transpose()

    if standard_scale:
        data = data.apply(preprocessing.scale)
        print(data.sum(axis=0))
    group_info = pd.read_csv(group_info, header=0, index_col=0, sep=None, engine='python')
    group_info = group_info.loc[[x for x in data.index if x in group_info.index]].iloc[:, 0]
    y_dict = dict()
    for ind, group in enumerate(set(group_info)):
        y_dict[group] = ind
    target = [y_dict[x] for x in group_info]
    print(data.shape)
    print(y_dict)
    data['ClassID'] = target
    return data


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['rf_result'])
