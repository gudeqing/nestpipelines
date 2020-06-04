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
from sklearn.inspection import permutation_importance
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.preprocessing import StandardScaler
from sklearn import decomposition, preprocessing
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
# We'll use this library to make the display pretty
from tabulate import tabulate
from sklearn.metrics import auc
from sklearn.metrics import plot_roc_curve
from sklearn.model_selection import StratifiedKFold
import pickle


def read_data(exp_matrix, group_info, scale_on_row=False, scale_on_col=False,  target_rows=None, target_cols=None):
    data = pd.read_csv(exp_matrix, header=0, index_col=0, sep=None, engine='python')
    data = data.fillna(0)
    if target_rows is not None:
        data = data.loc[[x.strip().split()[0] for x in open(target_rows)]]
        data.to_csv(exp_matrix+'.target.csv')
    if target_cols is not None:
        data = data[[x.strip().split()[0] for x in open(target_cols)]]
        data.to_csv(exp_matrix + '.target.csv')
    # print(((data.std(axis=1) / data.mean(axis=1)).abs()).describe())
    # data = data[(data.std(axis=1) / data.mean(axis=1)).abs() > 0.01]
    if scale_on_col:
        data = data.apply(preprocessing.scale)

    data = data.transpose()
    if scale_on_row:
        data = data.apply(preprocessing.scale)
    group_info = pd.read_csv(group_info, header=0, index_col=0, sep=None, engine='python')
    group_info = group_info.loc[[x for x in data.index if x in group_info.index]].iloc[:, 0]
    y_dict = dict()
    for ind, group in enumerate(set(group_info)):
        y_dict[group] = ind
    target = [y_dict[x] for x in group_info]
    print(data.shape)
    print(y_dict)
    return data, np.array(target)


def spearman_corr(X):
    from scipy.stats import spearmanr
    from scipy.cluster import hierarchy

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
    corr = spearmanr(X).correlation
    corr_linkage = hierarchy.ward(corr)
    dendro = hierarchy.dendrogram(corr_linkage, labels=X.columns, ax=ax1, leaf_rotation=90)
    dendro_idx = np.arange(0, len(dendro['ivl']))

    ax2.imshow(corr[dendro['leaves'], :][:, dendro['leaves']])
    ax2.set_xticks(dendro_idx)
    ax2.set_yticks(dendro_idx)
    ax2.set_xticklabels(dendro['ivl'], rotation='vertical')
    ax2.set_yticklabels(dendro['ivl'])
    fig.tight_layout()
    plt.savefig('spearman_corr.pdf')


def pca(data, explained_ratio=0.95):
    pca = decomposition.PCA()
    pca.fit(data)
    _ratio = list(enumerate(pca.explained_variance_ratio_, start=1))
    total_ratio, n_components = 0, 0
    for ind, each in _ratio:
        total_ratio += each
        if total_ratio >= explained_ratio:
            n_components = ind
            break
    if n_components <= 1:
        n_components = 2
    _ratio = _ratio[:n_components]
    result = pd.DataFrame(pca.transform(data), index=data.index)
    result = result.iloc[:, :n_components]
    return result


def train_and_predict(data, target, classifier="rf", max_iter:int=None, test_size=0.15, cv=10):
    """
    先把数据集拆分为两份，一份为训练集，一份为测试集，比例由test_size决定
    用上面的训练集进行grid search，grid search本身还需要一个交叉验证的参数cv
    这意味着每次训练时，还要对上面的训练集进一步划分，每次取其中的如90%进行训练，剩下的10%进行验证
    :param data:
    :param target:
    :param classifier:
    :param max_iter:
    :param test_size:
    :param cv:
    :return:
    """
    best_scores = [0]
    k = 0
    model_id = 0
    feature_num = data.shape[1]
    train_f = open(f'{classifier}.train_data_records.txt', 'a')
    test_f = open(f'{classifier}.test_data_record.txt', 'a')
    estimator_f = open(f'{classifier}.best_estimator.txt', 'a')
    my_scores = []
    feature_importance = dict()
    while True:
        # 随机分割数据，然后grid_search最佳模型
        # 对每次得到的最佳模型，都会画基于cross-validation的ROC曲线
        k += 1
        x_train, x_test, y_train, y_test = train_test_split(data, target, test_size=test_size)
        if classifier == 'rf':
            estimator = RandomForestClassifier()
            # 下面的字典给出参数的可选组合
            param_dict = dict(
                n_estimators=range(50, 200, 20),
                # criterion=('gini', 'entropy'),
                max_features=range(2, int(np.log2(feature_num)+np.sqrt(feature_num)), 2),
                oob_score=(True, False)
            )
        elif classifier == 'svm':
            estimator = SVC()
            param_dict = dict(
                kernel=('linear', 'poly', 'rbf', 'sigmoid'),
                degree=range(2, 7),
                gamma=(0.01, 0.1, 0.5, 1, 2),
                C=(1, 5, 10)
            )
        else:
            estimator = LogisticRegression()
            param_dict = dict(penalty=('l1', 'l2'))

        estimator = GridSearchCV(estimator, param_grid=param_dict, cv=cv, n_jobs=5)
        estimator.fit(x_train, y_train)
        # if estimator.best_score_ > best_scores[-1]:
        model_id += 1
        best_scores.append(estimator.best_score_)
        best_estimator = estimator.best_estimator_
        with open(f'best_model_{model_id}.pickle', 'wb') as f:
            pickle.dump(best_estimator, f)
        feature_importance[f'model_{model_id}'] = best_estimator.feature_importances_
        print(f'{model_id}.best estimator:', best_estimator, flush=True)
        estimator_f.write(f'{model_id}: '+str(best_estimator)+'\n')
        print(f'{model_id}.Mean cross-validated score of the best_estimator:', best_scores[-1], flush=True)
        train_f.write(f'{model_id}.mean_cv_score:{best_scores[-1]}\t'+'\t'.join(x_train.index)+'\n')
        y_predict = best_estimator.predict(x_test)
        test_accuracy = accuracy_score(y_test, y_predict)
        my_scores.append([model_id, best_estimator, best_scores[-1], test_accuracy])
        print(f'{model_id}.test_accuracy:{test_accuracy}', flush=True)
        test_f.write(f'{model_id}.test_accuracy:{test_accuracy}\t'+'\t'.join(x_test.index)+'\n')
        # final_x_train = x_train
        # final_y_train = y_train
        # final_x_test = x_test
        # final_y_test = y_test
        print(classification_report(y_test, y_predict), flush=True)
        # roc plot
        roc_cross_validation(best_estimator, data, target, out=f'model{model_id}.roc.pdf')

        if max_iter is None:
            # 倒数5个结果相差之和小于0.1
            if sum(abs(best_scores[-1] - x) for x in best_scores[-6:-1]) < 0.1:
                break
        else:
            if k == max_iter:
                break
        train_f.flush()
        test_f.flush()
        estimator_f.flush()

    train_f.close()
    test_f.close()
    estimator_f.close()
    model_id, best_estimator, score, test_accuracy = sorted(my_scores, key=lambda x:x[2]*x[3])[-1]
    print('performed grid searching {} times'.format(k))
    # print('Final Mean cross-validated score of the best_estimator:', best_scores[-1])
    print(f'Final best estimator is the {model_id}th:', best_estimator)
    print(f"It's score is {score}")
    print("It's accuracy on test data is:", test_accuracy)
    print("--Following is its performance on all input data---")
    y_predict = best_estimator.predict(data)
    print(classification_report(target, y_predict))
    model = best_estimator
    headers = ["name", "score"]
    values = sorted(zip(data.columns, model.feature_importances_), key=lambda x: x[1] * -1)
    print(tabulate(values, headers, tablefmt="plain"))
    # 输出所有最佳模型的feature importance，最佳模型有好有坏，均值不一定有意义
    df = pd.DataFrame(feature_importance, index=data.columns)
    df['mean_importance'] = df.mean(axis=1)
    df.sort_values(by='mean_importance', ascending=False).to_csv('feature_importance.csv')
    # 输出基于最优模型的permutation importance，对于共线性feature较多的数据，极有可能结果都是0
    result = permutation_importance(best_estimator, data, target, n_repeats=10, random_state=42, n_jobs=2)
    df = pd.DataFrame(result.importances, index=data.columns)
    df['mean_importance'] = df.mean(axis=1)
    df.sort_values(by='mean_importance', ascending=False).to_csv('permutation_feature_importance.csv')


def run(exp_matrix, group_info, classifier='rf', scale_on_row=False, scale_on_col=False,
        target_rows=None, target_cols=None, pca_first=False,
        max_iter:int=None, test_size=0.15, cv=10):
    data, target = read_data(exp_matrix, group_info, scale_on_row=scale_on_row, scale_on_col=scale_on_col,
                             target_rows=target_rows, target_cols=target_cols)
    if pca_first:
        data = pca(data, explained_ratio=0.95)
        print('after pca:', data.shape)
    train_and_predict(data, target, classifier=classifier, max_iter=max_iter, test_size=test_size, cv=cv)


def roc_cross_validation(classifier, X, y, out='roc.pdf', n_splits=5):
    """
    把数据集随机分成n_splits份，取n-1份作为训练集，剩下的一份作为测试集，共有n种取法。
    计算每次
    :param classifier:
    :param X:
    :param y:
    :param out:
    :return:
    """
    cv = StratifiedKFold(n_splits=n_splits)
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    fig, ax = plt.subplots()
    for i, (train, test) in enumerate(cv.split(X, y)):
        classifier.fit(X.iloc[train], y[train])
        viz = plot_roc_curve(classifier, X.iloc[test], y[test],
                             name='ROC fold {}'.format(i),
                             alpha=0.3, lw=1, ax=ax)
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)

    ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(mean_fpr, mean_tpr, color='b',
            label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
            lw=2, alpha=.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                    label=r'$\pm$ 1 std. dev.')

    ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="Receiver operating characteristic")
    ax.legend(loc="lower right")
    plt.savefig(out)
    plt.close()
    return mean_auc


def plot_decision_regions(X, y, classifier, test_idx=None, resolution=0.02, out='svm.decision.pdf'):
    # setup marker generator and color map
    markers = ('s', 'x', 'o', '^', 'v')
    colors = ('red', 'blue', 'lightgreen', 'gray', 'cyan')
    cmap = ListedColormap(colors[:len(np.unique(y))])

    # plot the decision surface
    x1_min, x1_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    x2_min, x2_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    xx1, xx2 = np.meshgrid(np.arange(x1_min, x1_max, resolution),
                           np.arange(x2_min, x2_max, resolution))
    Z = classifier.predict(np.array([xx1.ravel(), xx2.ravel()]).T)
    Z = Z.reshape(xx1.shape)
    plt.contourf(xx1, xx2, Z, alpha=0.4, cmap=cmap)
    plt.xlim(xx1.min(), xx1.max())
    plt.ylim(xx2.min(), xx2.max())

    for idx, cl in enumerate(np.unique(y)):
        plt.scatter(x=X[y == cl, 0], y=X[y == cl, 1],
                    alpha=0.8, c=cmap(idx),
                    marker=markers[idx], label=cl)

    # highlight test samples
    if test_idx:
        # plot all samples
        X_test, y_test = X[test_idx, :], y[test_idx]
        plt.scatter(X_test[:, 0],
                    X_test[:, 1],
                    c='',
                    alpha=1.0,
                    linewidths=1,
                    marker='o',
                    s=55, label='test set')
    plt.savefig(out)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['run'])
