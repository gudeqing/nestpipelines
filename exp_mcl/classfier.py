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
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.preprocessing import StandardScaler
from sklearn import decomposition, preprocessing
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
# We'll use this library to make the display pretty
from tabulate import tabulate


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
    return data, target


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


def train_and_predict(data, target, classifier="rf", max_iter=None, test_size=0.15, cv=10):
    best_scores = [0]
    k = 0
    model_id = 0
    feature_num = data.shape[1]
    train_f = open(f'{classifier}.train_data_records.txt', 'a')
    test_f = open(f'{classifier}.test_data_record.txt', 'a')
    estimator_f = open(f'{classifier}.best_estimator.txt', 'a')
    while True:
        k += 1
        x_train, x_test, y_train, y_test = train_test_split(data, target, test_size=test_size)
        if classifier == 'rf':
            estimator = RandomForestClassifier()
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
        if estimator.best_score_ > best_scores[-1]:
            model_id += 1
            best_scores.append(estimator.best_score_)
            best_estimator = estimator.best_estimator_
            print(f'{model_id}.best estimator:', best_estimator, flush=True)
            estimator_f.write(f'{model_id}: '+str(best_estimator)+'\n')
            print(f'{model_id}.Mean cross-validated score of the best_estimator:', best_scores[-1], flush=True)
            train_f.write(f'{model_id}.mean_cv_score:{best_scores[-1]}\t'+'\t'.join(x_train.index)+'\n')
            y_predict = best_estimator.predict(x_test)
            test_accuracy = accuracy_score(y_test, y_predict)
            print(f'{model_id}.test_accuracy:{test_accuracy}', flush=True)
            test_f.write(f'{model_id}.test_accuracy:{test_accuracy}\t'+'\t'.join(x_test.index)+'\n')
            final_x_train = x_train
            final_y_train = y_train
            final_x_test = x_test
            final_y_test = y_test
            print(classification_report(final_y_test, y_predict), flush=True)
            if max_iter is None:
                if sum(best_scores[-1] - x for x in best_scores[-6:-1]) < 0.05:
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
    print('performed grid searching {} times'.format(k))
    print('Final Mean cross-validated score of the best_estimator:', best_scores[-1])
    print('best estimator:', best_estimator)
    print('for test data:')
    y_predict = best_estimator.predict(final_x_test)
    print(classification_report(final_y_test, y_predict))
    model = best_estimator
    headers = ["name", "score"]
    values = sorted(zip(final_x_train.columns, model.feature_importances_), key=lambda x: x[1] * -1)
    print(tabulate(values, headers, tablefmt="plain"))


def run(exp_matrix, group_info, classifier='rf', standard_scale=False,
        target_rows=None, target_cols=None, pca_first=False,
        max_iter=None, test_size=0.15, cv=10):
    data, target = read_data(exp_matrix, group_info, standard_scale=standard_scale,
                             target_rows=target_rows, target_cols=target_cols)
    if pca_first:
        data = pca(data, explained_ratio=0.95)
        print('after pca:', data.shape)
    train_and_predict(data, target, classifier=classifier, max_iter=max_iter, test_size=test_size, cv=cv)


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
