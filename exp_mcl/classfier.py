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
from boruta import BorutaPy
import pickle


def read_data(exp_matrix, group_info, scale_on_row=False, scale_on_col=False,  target_rows=None, target_cols=None):
    """
    读取目标数据，用于模型训练和特征选择
    :param exp_matrix:
    :param group_info:
    :param scale_on_row:
    :param scale_on_col:
    :param target_rows:
    :param target_cols:
    :return: 返回X（仍未dataframe）和y（为numpy.array)
    """
    data = pd.read_csv(exp_matrix, header=0, index_col=0, sep=None, engine='python')
    data = data.fillna(0)
    if target_rows is not None:
        targets = [x.strip().split()[0] for x in open(target_rows)]
        data = data.loc[[x for x in targets if x in data.index]]
        data.to_csv(exp_matrix+'.target.csv')
    if target_cols is not None:
        targets = [x.strip().split()[0] for x in open(target_cols)]
        data = data[[x for x in targets if x in data.columns]]
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


def group_collinear_vars(X, corr_cutoff=0.8, target_rows=None, distance_cutoff=0.0, method='spearman'):
    from scipy.cluster import hierarchy
    if type(X) == str:
        X = pd.read_csv(X, header=0, index_col=0, sep=None, engine='python')
    if target_rows is not None:
        targets = [x.strip().split()[0] for x in open(target_rows)]
        X = X.loc[[x for x in targets if x in X.index]]
        # X.to_csv('target_data.csv')
    corr = X.T.corr(method=method)

    # 计算高相关性的基因
    high_corr_group = set()
    for ind, row in corr.abs().iterrows():
        group = set(row[row > corr_cutoff].index)
        if len(group) > 1:
            high_corr_group.add(tuple(sorted(group)))

    # merge high corr group
    merged = set()
    for each in high_corr_group:
        tmp = set(each)
        _ = [tmp.update(set(x)) for x in high_corr_group if set(x)&tmp]
        merged.add(tuple(sorted(tmp)))
    with open('HighCorrGroup.list', 'w') as f:
        for each in merged:
            f.write(each[0]+ '\t' + ';'.join(each)+'\n')

    # 计算VIF，计算出与其他基因有较高共线性的基因
    vif = np.linalg.inv(corr.values).diagonal()
    vif_df = pd.DataFrame(vif, index=corr.index, columns=['VIF']).sort_values(by='VIF', ascending=False)
    vif_df.round(2).to_csv('all.vif.txt', sep='\t')
    vif_dict = dict(zip(vif_df.index, vif_df['VIF']))
    # high_vif_genes = vif_df[vif_df['VIF'] >= 5].index

    # 根据聚类距离计算距离比较近的基因簇
    corr_linkage = hierarchy.ward(corr)
    distance_desc = pd.DataFrame(corr_linkage)[2].describe()
    if distance_cutoff <= 0:
        lq = distance_desc['25%']
        print('25% distance cutoff:', lq)
    else:
        lq = distance_cutoff
    cluster_ids = hierarchy.fcluster(corr_linkage,lq, criterion='distance')
    distance_group = dict()
    ids = corr.columns
    for ind, big_id in enumerate(cluster_ids):
        distance_group.setdefault(ids[big_id], list())
        distance_group[ids[big_id]].append(ids[ind])

    # 找代表
    representing = list()
    for _, v in distance_group.items():
        if len(v) > 1:
            vif_lst = [vif_dict[x] for x in v]
            # 选择vif最大的一个作为代表
            represent = sorted(zip(v, vif_lst), key=lambda x:x[1])[-1][0]
            # 如果v中某个变量与代表的相关系数低于0.7，则该变量不能被代表, 只能自己作为代表
            # print(represent, v)
            for x, corr_value in corr.loc[represent, v].iteritems():
                if abs(corr_value) < 0.7:
                    v.remove(x)
                    representing.append([x, x, 1, vif_dict[x]])
            if len(v) > 1:
                inter_corr = corr.loc[v, v]
                mean_corr = (inter_corr.values.flatten().sum() - len(v)) / (len(v) ** 2 - len(v))
            else:
                mean_corr = 1
            vif_lst = [vif_dict[x] for x in v]
            mean_vif = sum(vif_lst) / len(vif_lst)
            representing.append([represent, ';'.join(sorted(v)), mean_corr, mean_vif])
        else:
            representing.append([v[0], v[0], 1, vif_dict[v[0]]])

    with open('representing_decision.txt', 'w') as f:
        f.write('represent\tmembers\tmean_corr\tVIF\n')
        for each in representing:
            f.write('\t'.join(str(x) for x in each)+'\n')

    # new VIF using representing
    new_corr = X.loc[[x[0] for x in representing]].T.corr(method=method)
    vif = np.linalg.inv(new_corr.values).diagonal()
    vif_df = pd.DataFrame(vif, index=new_corr.index, columns=['VIF']).sort_values(by='VIF', ascending=False)
    vif_df.round(2).to_csv('representing.vif.txt', sep='\t')

    # 可视化聚类结果
    # 如果筛选的基因全是差异基因，在不对相关系数取绝对值时，聚成两大类，分别对应上调基因和下调基因成团
    # 如果对相关系数取绝对值，可以聚成3大类，即增加了一个负相关基因成团的聚类
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
    dendro = hierarchy.dendrogram(corr_linkage, ax=ax1,  leaf_rotation=90)
    ax2.imshow(corr.iloc[dendro['leaves'], :].iloc[:, dendro['leaves']])
    ax1.set_xticklabels(corr.index[dendro['leaves']])
    fig.tight_layout()
    plt.savefig(f'{method}_corr_cluster.pdf')
    return representing


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


def multi_grid_search_model(data, target, classifier="rf", max_iter:int=None, test_size=0.15, cv=10, prefix=''):
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
    :param prefix: 输出结果前缀，目前仅支持所有结果输出到当前目录
    :return:返回确定好参数的最佳模型
    """
    best_scores = [0]
    k = 0
    model_id = 0
    feature_num = data.shape[1]
    train_f = open(f'{prefix}{classifier}.train_data_records.txt', 'a')
    test_f = open(f'{prefix}{classifier}.test_data_record.txt', 'a')
    estimator_f = open(f'{prefix}{classifier}.best_estimator.txt', 'a')
    my_scores = []
    feature_importance = dict()
    while True:
        # 随机分割数据，然后grid_search最佳模型参数
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
        with open(f'{prefix}best_model_{model_id}.pickle', 'wb') as f:
            pickle.dump(best_estimator, f)

        if classifier == 'rf':
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
        if len(set(target)) == 2:
            roc_cross_validation(best_estimator, data, target, out=f'{prefix}model{model_id}.roc.pdf')

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
    if classifier == 'rf':
        model = best_estimator
        headers = ["name", "score"]
        values = sorted(zip(data.columns, model.feature_importances_), key=lambda x: x[1] * -1)
        print(tabulate(values, headers, tablefmt="plain"))
        # 输出所有最佳模型的feature importance，最佳模型有好有坏，均值不一定有意义
        df = pd.DataFrame(feature_importance, index=data.columns)
        df['mean_importance'] = df.mean(axis=1)
        df.sort_values(by='mean_importance', ascending=False).to_csv(f'{prefix}feature_importance.csv')
        # 输出基于最优模型的permutation importance，对于共线性feature较多的数据，极有可能结果都是0
        result = permutation_importance(best_estimator, data, target, n_repeats=10, random_state=42, n_jobs=2)
        df = pd.DataFrame(result.importances, index=data.columns)
        df['mean_importance'] = df.mean(axis=1)
        df.sort_values(by='mean_importance', ascending=False).to_csv(f'{prefix}permutation_feature_importance.csv')
    return best_estimator


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


def run(exp_matrix, group_info, classifier='rf',
        scale_on_row=False, scale_on_col=False,  # data preprocess
        target_rows=None, target_cols=None, pca_first=False,  # data selection
        max_iter:int=10, test_size=0.15, cv=10, # grid search model
        no_feature_selection=False,  # if use boruta select feature, Need a randomforest classfier
        percent=90, alpha=0.05  # for boruta feature selection
        ):
    """
    :param exp_matrix:
    :param group_info:
    :param classifier:
    :param scale_on_row:
    :param scale_on_col:
    :param target_rows:
    :param target_cols:
    :param pca_first:
    :param max_iter:
    :param test_size:
    :param cv:
    :param no_feature_selection:
    :param percent:
    :param alpha:
    :return:
    """
    # step1: preprocess
    X, y = read_data(
        exp_matrix, group_info,
        scale_on_row=scale_on_row,
        scale_on_col=scale_on_col,
        target_rows=target_rows,
        target_cols=target_cols
    )
    if pca_first:
        X = pca(X, explained_ratio=0.95)
        print('after pca:', X.shape)

    # step2: optimize model parameter
    param_optimized_model = multi_grid_search_model(
        X, y, classifier=classifier,
        max_iter=max_iter, test_size=test_size, cv=cv
    )

    # step3: feature selection
    if not no_feature_selection:
        feat_selector = BorutaPy(
            param_optimized_model,
            n_estimators='auto', verbose=2, random_state=1,
            max_iter=100, perc=percent,
        )
        # find all relevant features
        feat_selector.fit(X.values, y)
        # check selected features - first 5 features are selected
        print('selected feature:\n', X.columns[feat_selector.support_])
        # check ranking of features
        rank = pd.DataFrame({'feature':X.columns, 'rank': feat_selector.ranking_})
        rank.sort_values(by='rank').to_csv('feature_ranking.txt', sep='\t')

        # step4: 找共线性并举出代表
        # re-train model by using selected feature
        X = X.iloc[:, feat_selector.support_]
        represents = group_collinear_vars(
            X.T, corr_cutoff=0.8, method='spearman'
        )

        # step5: 使用最终筛选出来的代表性feature进行最终的模型训练和评估
        X = X[[x[0] for x in represents]]
        multi_grid_search_model(
            X, y, classifier=classifier, prefix='final.',
            max_iter=max_iter, test_size=test_size, cv=cv
        )

        # save final data used for building model
        data = X.T
        rep_dict = dict((x[0], x[1]) for x in represents)
        data.insert(0, 'members', [rep_dict[x] for x in data.index])
        data.to_csv('final_model_based_data.txt')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['run', 'group_collinear_vars'])
