import os
import matplotlib
matplotlib.use('agg')
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import pickle
from tabulate import tabulate

from sklearn import preprocessing
from sklearn import decomposition

from sklearn.ensemble import RandomForestClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LogisticRegressionCV

from sklearn.svm import SVC
# feature selection module
from boruta import BorutaPy

from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold

from sklearn.metrics import auc
from sklearn.metrics import roc_curve
from sklearn.metrics import plot_roc_curve
from sklearn.metrics import plot_confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.inspection import permutation_importance


def read_data(exp_matrix, group_info, scale_on_row=False, scale_on_col=False,  target_rows=None, target_cols=None):
    """
    读取目标数据，用于模型训练和特征选择
    :param exp_matrix:
    :param group_info:
    :param scale_on_row:
    :param scale_on_col:
    :param target_rows:
    :param target_cols:
    :return: 返回X（dataframe）和y（numpy.array)
    """
    if type(exp_matrix) == str:
        data = pd.read_csv(exp_matrix, header=0, index_col=0, sep=None, engine='python')
    else:
        data = exp_matrix
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
    for ind, group in enumerate(sorted(set(group_info))):
        y_dict[group] = ind
    target = [y_dict[x] for x in group_info]
    print(data.shape)
    print(y_dict)
    return data, np.array(target)


def group_collinear_vars(X, corr_cutoff=0.75, target_rows=None, distance_cutoff=0.0, method='spearman'):
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
    with open('all.HighCorrGroup.xls', 'w') as f:
        for each in merged:
            f.write(each[0]+ '\t' + ';'.join(each)+'\n')

    # 计算VIF，计算出与其他基因有较高共线性的基因
    vif = np.linalg.inv(corr.values).diagonal()
    vif_df = pd.DataFrame(vif, index=corr.index, columns=['VIF']).sort_values(by='VIF', ascending=False)
    vif_df.round(2).to_csv('all.vif.xls', sep='\t')
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

    # 找代表, 上面已经根据聚类聚类将feature进行了小组划分，选择VIF最大得作为代表
    # 接下来，如果某个成员和代表的相关性不够大，则独立出去
    representing = list()
    for _, v in distance_group.items():
        if len(v) > 1:
            vif_lst = [vif_dict[x] for x in v]
            # 选择vif最大的一个作为代表
            represent = sorted(zip(v, vif_lst), key=lambda x:x[1])[-1][0]
            # 如果v中某个变量与代表的相关系数低于0.7，则该变量不能被代表, 只能自己作为代表
            # print(represent, v)
            for x, corr_value in corr.loc[represent, v].iteritems():
                if abs(corr_value) < corr_cutoff:
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

    representing_df = pd.DataFrame(representing, columns=['represent', 'members', 'mean_corr', 'VIF_before'])
    representing_df.set_index('represent', inplace=True)
    # new VIF using representing
    new_corr = X.loc[[x[0] for x in representing]].T.corr(method=method)
    vif = np.linalg.inv(new_corr.values).diagonal()
    vif_df = pd.DataFrame(vif, index=new_corr.index, columns=['VIF_after'])
    representing_df = representing_df.join(vif_df).sort_values(by='VIF_before', ascending=False)
    representing_df.round(2).to_csv('feature_represent_decision.xls', sep='\t')

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


def multi_grid_search_model(data, target, classifier="rf", max_iter:int=None, ovr=False, test_size=0.15, cv=10, prefix=''):
    """
    先把数据集拆分为两份，一份为训练集，一份为测试集，比例由test_size决定
    用上面的训练集进行grid search，grid search本身还需要一个交叉验证的参数cv
    这意味着每次训练时，还要对上面的训练集进一步划分，每次取其中的如90%进行训练，剩下的10%进行验证
    :param data:
    :param target:
    :param classifier:
    :param max_iter: 随机分割训练集和测试集的次数，如设为10，表示进行10次参数调优，然后从10个最优中再选最优。
    :param ovr: 多分类时，是否采用OneVsRestClassifier策略
    :param test_size:
    :param cv:
    :param prefix: 输出结果前缀，目前仅支持所有结果输出到当前目录
    :return:返回确定好参数的最佳模型
    """
    best_scores = [0]
    k = 0
    model_id = 0
    feature_num = data.shape[1]
    # train_f = open(f'{prefix}{classifier}.train_data_records.txt', 'a')
    # test_f = open(f'{prefix}{classifier}.test_data_record.txt', 'a')
    # estimator_f = open(f'{prefix}{classifier}.best_estimator.txt', 'a')
    report_f = open(f'{prefix}{classifier}.cls_report.txt', 'w')
    my_scores = []
    feature_importance = dict()
    outdir = f'{prefix}models'
    os.makedirs(outdir, exist_ok=True)
    while True:
        # 随机分割数据，然后grid_search最佳模型参数
        # 对每次得到的最佳模型，都会画基于cross-validation的ROC曲线
        k += 1
        x_train, x_test, y_train, y_test = train_test_split(data, target, test_size=test_size, stratify=target)
        if classifier == 'rf':
            estimator = RandomForestClassifier()
            if ovr and len(set(y_test)) > 2:
                estimator = OneVsRestClassifier(estimator)
                # 下面的字典给出参数的可选组合
                param_dict = dict(
                    estimator__n_estimators=range(50, 200, 20),
                    estimator__criterion=('gini', ),
                    estimator__max_features=range(2, int(np.log2(feature_num) + np.sqrt(feature_num)), 2),
                    estimator__oob_score=(True, False),
                    estimator__max_depth=[5, ]
                )
            else:
                # 下面的字典给出参数的可选组合
                param_dict = dict(
                    n_estimators=range(50, 200, 20),
                    # criterion=('gini', 'entropy'),
                    max_features=range(2, int(np.log2(feature_num)+np.sqrt(feature_num)), 2),
                    oob_score=(True, False),
                    max_depth=[5,]
                )
        elif classifier == 'svm':
            estimator = SVC()
            param_dict = dict(
                kernel=('linear', 'poly', 'rbf', 'sigmoid'),
                degree=range(2, 7),
                gamma=(0.01, 0.1, 0.5, 1, 2),
                C=(1, 5, 10)
            )
        elif classifier == 'lgr':
            estimator = LogisticRegression(warm_start=True, max_iter=100)
            param_dict = dict(penalty=('l2',))
        else:
            raise Exception('Only support rf or svm or lgr classifier!')

        estimator = GridSearchCV(estimator, param_grid=param_dict, cv=cv, n_jobs=5)
        estimator.fit(x_train, y_train)
        # if estimator.best_score_ > best_scores[-1]:
        model_id += 1
        best_scores.append(estimator.best_score_)
        best_estimator = estimator.best_estimator_
        with open(f'{outdir}/{prefix}best_model_{model_id}.pickle', 'wb') as f:
            pickle.dump(best_estimator, f)

        if classifier == 'rf' and not ovr:
            feature_importance[f'model_{model_id}'] = best_estimator.feature_importances_
        model_str = best_estimator.__str__().replace("\n", '').replace(' ', '')
        print(f'{model_id}.BestModel<mean_cv_score={best_scores[-1]:.2f}>:', model_str, flush=True)
        # estimator_f.write(f'{model_id}: '+str(best_estimator)+'\n')
        # train_f.write(f'{model_id}.mean_cv_score:{best_scores[-1]}\t'+'\t'.join(x_train.index)+'\n')
        y_predict = best_estimator.predict(x_test)
        missed = x_test.loc[[y1 != y2 for y1, y2 in zip(y_predict, y_test)]].index
        print(f'Miss labeled samples by best estimator are:', missed)
        test_accuracy = accuracy_score(y_test, y_predict)
        # print(f'{model_id}.test_accuracy:{test_accuracy}', flush=True)
        # test_f.write(f'{model_id}.test_accuracy:{test_accuracy}\t'+'\t'.join(x_test.index)+'\n')
        print(
            f'model{model_id}[score:{best_scores[-1]:.2f}] performance on test data:\n',
            classification_report(y_test, y_predict), flush=True, file=report_f
        )
        # plot confusion matrix
        plot_confusion_table(best_estimator, x_test, y_test, out=f'{outdir}/{prefix}model{model_id}.test.confusion.pdf')
        plot_confusion_table(best_estimator, data, target, out=f'{outdir}/{prefix}model{model_id}.all.confusion.pdf')

        # roc plot
        if len(set(target)) == 2:
            # mean_auc = roc_cross_validation(best_estimator, data, target, out=f'{outdir}/{prefix}model{model_id}.roc.pdf')
            mean_auc = binary_roc_plot(best_estimator, x_test, y_test, out=f'{outdir}/{prefix}model{model_id}.roc.pdf')
        else:
            mean_auc = multiclass_roc_plot(best_estimator, x_train, y_train, x_test, y_test, out=f'{outdir}/{prefix}model{model_id}.roc.pdf')
        my_scores.append([model_id, best_estimator, best_scores[-1], test_accuracy, mean_auc])
        if max_iter is None:
            # 倒数5个结果相差之和小于0.1
            if sum(abs(best_scores[-1] - x) for x in best_scores[-6:-1]) < 0.1:
                break
        else:
            if k == max_iter:
                break
        # train_f.flush()
        # test_f.flush()
        # estimator_f.flush()

    # train_f.close()
    # test_f.close()
    # estimator_f.close()
    sort_scores = sorted(my_scores, key=lambda x: x[2] * x[3], reverse=True)
    model_id, best_estimator, score, test_accuracy, mean_auc = sort_scores[0]
    # 整理临时文件
    all_best_model_info = pd.DataFrame(
        sort_scores,
        columns=['model_id', 'model', 'mean_cv_score', 'test_accuracy', 'mean_auc']
    )
    all_best_model_info['model'] = [x.__str__().replace("\n", '').replace(' ', '') for x in all_best_model_info['model']]
    all_best_model_info.to_csv(f'{prefix}best_models.xls', sep='\t', index=False)
    print('performed grid searching {} times'.format(k), file=report_f)
    # print('Final Mean cross-validated score of the best_estimator:', best_scores[-1])
    model_str = best_estimator.__str__().replace("\n", '').replace(' ', '')
    print(f'Final best estimator is the {model_id}th:', model_str, file=report_f)
    print(f"It's mean cv score is {score}", file=report_f)
    print("It's accuracy on test data is:", test_accuracy, file=report_f)
    print("It's AUC on test data is:", mean_auc, file=report_f)
    print("--Following is its performance on all input data---", file=report_f)
    y_predict = best_estimator.predict(data)
    print(classification_report(target, y_predict), file=report_f)
    report_f.close()

    if classifier == 'rf' and not ovr:
        # model = best_estimator
        # headers = ["name", "score"]
        # values = sorted(zip(data.columns, model.feature_importances_), key=lambda x: x[1] * -1)
        # print(tabulate(values, headers, tablefmt="plain"))
        # 输出所有最佳模型的feature importance，最佳模型有好有坏，均值不一定有意义
        df = pd.DataFrame(feature_importance, index=data.columns)
        df['mean_importance'] = df.mean(axis=1)
        df.sort_values(by='mean_importance', ascending=False).to_csv(f'{outdir}/{prefix}feature_importance.csv')
        # 输出基于最优模型的permutation importance，对于共线性feature较多的数据，极有可能结果都是0
        result = permutation_importance(best_estimator, data, target, n_repeats=10, random_state=42, n_jobs=2)
        df = pd.DataFrame(result.importances, index=data.columns)
        df['mean_importance'] = df.mean(axis=1)
        df.sort_values(by='mean_importance', ascending=False).to_csv(f'{outdir}/{prefix}permutation_feature_importance.csv')

    return best_estimator


def binary_roc_plot(clf, X_test, y_test, out='roc.pdf'):
    fig, ax = plt.subplots()
    viz = plot_roc_curve(clf, X_test, y_test, alpha=0.8, lw=1.5, ax=ax, drop_intermediate=False)
    ax.plot([0, 1], [0, 1], linestyle='--', lw=1, label='Chance', alpha=0.3)
    fig.tight_layout()
    fig.savefig(out)
    plt.close()
    return viz.roc_auc


def roc_cross_validation(classifier, X, y, out='roc.pdf', n_splits=5, shuffle=True):
    """
    把数据集随机分成n_splits份，取n-1份作为训练集，剩下的一份作为测试集，共有n种取法。
    计算每次
    :param classifier:
    :param X:
    :param y:
    :param out:
    :return:
    """
    cv = StratifiedKFold(n_splits=n_splits, shuffle=shuffle)
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    fig, ax = plt.subplots()
    for i, (train, test) in enumerate(cv.split(X, y)):
        classifier.fit(X.iloc[train], y[train])
        viz = plot_roc_curve(classifier, X.iloc[test], y[test],
                             name='ROC fold {}'.format(i),
                             alpha=0.3, lw=1, ax=ax)
        pred = classifier.predict(X.iloc[test])
        missed = X.iloc[test].loc[[y1!=y2 for y1,y2 in zip(pred, y[test])]].index
        print(f'Miss labeled sample in {i}th fold:', missed)
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

    ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title=f"{n_splits}-cross validation ROC")
    ax.legend(loc="lower right")
    plt.savefig(out)
    plt.clf()
    return mean_auc


def plot_confusion_table(clf, X_test, y_test, normalize=None, cmap='bwr', out='confusion.pdf'):
    # plt.figure()
    plot_confusion_matrix(clf, X_test, y_test, normalize=normalize, cmap=cmap)
    plt.tight_layout()
    plt.savefig(out)
    plt.clf()


def multiclass_roc_plot(clf, X_train, y_train, X_test, y_test, out='multiLabel.roc.pdf'):
    # https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc.html
    # classifier
    if not clf.__str__().startswith('OneVsRestClassifier'):
        print('当前模型不是OneVsRestClassifier, 为了画ROC图，将会使用OneVsRestClassifier重新训练和预测')
        print('OneVsRestClassifier和原始模型的性能可能不一样，因此ROC图仅供参考')
        clf = OneVsRestClassifier(clf)
        y_score = clf.fit(X_train, y_train).predict_proba(X_test)
    else:
        y_score = clf.predict_proba(X_test)
        # print(y_score, clf.predict(X_test))
        # 下面的代码测试发现，再次将OneVsRestClassifier嵌入OneVsRestClassifier后，结果或有提升或有下降
        # print('before', accuracy_score(y_test, clf.predict(X_test)))
        # clf = OneVsRestClassifier(clf)
        # clf.fit(X_train, y_train)
        # print('after', accuracy_score(y_test, clf.predict(X_test)))

    y_test = preprocessing.label_binarize(y_test, classes=clf.classes_)
    if y_test.shape[1] == 1:
        y_test = np.hstack((y_test, (~y_test.astype(bool)).astype(int)))
    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i, cls in enumerate(clf.classes_):
        fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i], drop_intermediate=False)
        roc_auc[i] = auc(fpr[i], tpr[i])

    # Plot of a ROC curve for a specific class
    # plt.figure()
    for i, cls in enumerate(clf.classes_):
        plt.plot(fpr[i], tpr[i], label=f'{cls}:AUC={roc_auc[i]:.2f}')
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC of OneVsRestClassifier')
        plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(out)
    plt.clf()
    return sum(roc_auc.values())/len(roc_auc)


def get_color_pool(n):
    # https://plot.ly/ipython-notebooks/color-scales/
    import colorlover
    if n <= 8:
        if n <= 3:
            n = 3
        return colorlover.scales[str(n)]['qual']['Set1']
    if n <= 12:
        return colorlover.scales[str(n)]['qual']['Paired']


def per_logit_reg(exp_matrix, group_info, target_rows=None, target_cols=None, min_auc=0.5,
                  link='https://www.proteinatlas.org/search/{}'):
    from bokeh.plotting import figure, save, output_file
    from bokeh.models import ColumnDataSource, CDSView, GroupFilter, \
        HoverTool, LabelSet, Legend, CustomJS, TapTool, OpenURL,Div
    from bokeh.models.annotations import Title
    from bokeh.events import ButtonClick
    from bokeh.layouts import gridplot

    X, y = read_data(
        exp_matrix=exp_matrix,
        group_info=group_info,
        target_cols=target_cols,
        target_rows=target_rows
    )
    plot_options = dict(
        width=250,
        plot_height=250,
        # tools='pan,wheel_zoom,box_select,reset,save,tap',  # tap 支持点击打开超链接
        tools='wheel_zoom,reset,save,tap',  # tap 支持点击打开超链接
    )
    plots = []
    report = []
    if len(set(y)) == 2:
        plot_data = []
        for ind, col in enumerate(X.columns):
            clf = LogisticRegressionCV(cv=5, random_state=0)
            train = X[col].values.reshape(-1, 1)
            clf.fit(train, y)
            y_pred = clf.predict_proba(train)[:, 1]
            fpr, tpr, threshold = roc_curve(y, y_pred, drop_intermediate=False)
            roc_auc = auc(fpr, tpr)
            report.append([col, roc_auc])
            if roc_auc >= min_auc:
                plot_data.append((col, fpr, tpr, roc_auc))

        plot_data = sorted(plot_data, key=lambda x:x[3], reverse=True)
        for col, fpr, tpr, roc_auc in plot_data:
            s = figure(**plot_options, title=f'{col}  AUC={roc_auc:.2f}')
            url = link.format(col)
            taptool = s.select(type=TapTool)
            taptool.callback = OpenURL(url=url)
            source = ColumnDataSource(data=dict(fpr=fpr, tpr=tpr))
            hover = HoverTool(
                tooltips=[
                    ("FPR", "@{}".format('fpr')),
                    ("TPR", "@{}".format('tpr')),
                ]
            )
            s.add_tools(hover)
            s.line('fpr', 'tpr', color='blue', legend_label=f'{col}  AUC={roc_auc:.2f}', source=source)
            s.xaxis.axis_label = 'FPR'
            s.yaxis.axis_label = 'TPR'
            s.legend.location = 'bottom_right'
            plots.append((s, roc_auc))
    else:
        plot_data = []
        for col in X.columns:
            clf = OneVsRestClassifier(LogisticRegressionCV(cv=5, random_state=0))
            train = X[col].values.reshape(-1, 1)
            y_score = clf.fit(train, y).predict_proba(train)
            y_test = preprocessing.label_binarize(y, classes=clf.classes_)
            # Compute ROC curve and ROC area for each class
            auc_lst = []
            fpr_lst = []
            tpr_lst = []
            cls_lst = []
            for i, cls in enumerate(clf.classes_):
                fpr, tpr, _ = roc_curve(y_test[:, i], y_score[:, i], drop_intermediate=False)
                roc_auc = auc(fpr, tpr)
                auc_lst.append(roc_auc)
                fpr_lst.append(fpr)
                tpr_lst.append(tpr)
                cls_lst.append(cls)
            mean_auc = sum(auc_lst) / len(auc_lst)
            report.append([col, mean_auc])
            if mean_auc >= min_auc:
                plot_data.append([col, cls_lst, fpr_lst, tpr_lst, auc_lst, mean_auc])

        plot_data = sorted(plot_data, key=lambda x:x[4], reverse=True)
        colors = get_color_pool(len(set(y)))
        for col, cls_lst, fpr_lst, tpr_lst, auc_lst, mean_auc in plot_data:
            s = figure(**plot_options)
            title = Title()
            title.text = f'{col}  Mean_AUC={mean_auc:.2f}'
            s.title = title
            url = link.format(col)
            taptool = s.select(type=TapTool)
            taptool.callback = OpenURL(url=url)

            for i, (cls, fpr, tpr, roc_auc) in enumerate(zip(cls_lst, fpr_lst, tpr_lst, auc_lst)):
                source = ColumnDataSource(data=dict(fpr=fpr, tpr=tpr))
                s.line('fpr', 'tpr', color=colors[i], legend_label=f'{cls}  AUC={roc_auc:.2f}', source=source)

            hover = HoverTool(
                tooltips=[
                    ("FPR", "@{}".format('fpr')),
                    ("TPR", "@{}".format('tpr')),
                ]
            )
            s.add_tools(hover)
            s.xaxis.axis_label = 'FPR'
            s.yaxis.axis_label = 'TPR'
            s.legend.location = 'bottom_right'
            plots.append((s, mean_auc))
    # save
    result = pd.DataFrame(report, columns=['feature', 'AUC']).sort_values(by='AUC', ascending=False)
    result.to_csv('all.feature_auc.xls', sep='\t', index=False)
    # plot all
    if plots:
        plots = [x[0] for x in sorted(plots, key=lambda x:x[1], reverse=True)]
        p = gridplot(plots, sizing_mode='stretch_{}'.format('width'), ncols=4, toolbar_location='left')
        output_file(f'top{len(plots)}.ROC.html', title="ROC")
        save(p)


def run(exp_matrix, group_info, classifier='rf',
        scale_on_row=False, scale_on_col=False,  # data preprocess
        target_rows=None, target_cols=None, pca_first=False,  # data selection
        grid_search_num:int=10, test_size=0.15, cv=10,  # grid search model
        no_feature_selection=False,  # if use boruta select feature, Need a randomforest classfier
        percent:int=90, alpha=0.05,  # for boruta feature selection
        corr_cutoff=0.75,
        slgr='yes',
        tentative=False,
        ovr=False,
        link='https://www.proteinatlas.org/search/{}'
        ):
    """
    1. 全变量模型参数优化
    网格搜索参数时有一个cv参数，在进行网格搜索前也有一次随机数据拆分，所以数据实际拆分了两次。
    也就是说，原始数据划分为3份，分别为：训练集、验证集和测试集；其中训练集用来模型训练，验证集用来调整参数，而测试集用来衡量模型表现好坏。
    本次程序设计进行N(=grid_search_num)次网格搜索，并根据模型score和测试集准确性accuracy的乘积进行排序，对乘积最高模型用ROC分析。
    2. 基于上述优化后的模型，使用borutapy进行feature筛选
    3. 对筛选的变量进行相关性分析，并对相关性较强的变量进行分组，每组派出一个代表作为最后模型使用的feature
    4. 基于最后selected features的模型参数优化
    :param exp_matrix: 表达矩阵文件，每一行表示一个基因在n个样本中的表达丰度
    :param group_info: 分组信息，第一列为样本名
    :param classifier: one of ['rf', 'svm', 'lgr']
    :param scale_on_row:
    :param scale_on_col:
    :param target_rows: 文件的第一列作为要提取的目标行，无header
    :param target_cols: 文件的第一列作为要提取的目标列，无header
    :param pca_first: 是否先使用pca进行降维，降维后的数据用于后续分析
    :param grid_search_num: 指定进行多少次网格搜索. 因为网格搜索前需对数据进行随机拆分，不同拆分可能导致搜索到的最优模型不一样。
    :param test_size: 网格调参前的数据拆分比例参数，百分比，表示使用多少比例用于参数优化后的模型
    :param cv: 网格搜素最优参数时的cv参数，
    :param no_feature_selection: 是否使用boruta进行feature筛选
    :param percent: boruta的参数，0-100，100时最严格，默认90
    :param alpha: boruta的参数，即FDR阈值，默认0.05
    :param corr_cutoff: 相关系数阈值，用于寻找相关性较强的变量
    :param slgr: 对每个变量进行一次逻辑回归分析并可视化
    :param tentative:使用borutapy筛选feature时, 如果设置该参数,则把tentative变量也包括进来
    :param ovr: 多分类时，是否采用OneVsRestClassifier策略. 目前该参数仅针对随机森林有效
        有文献报道使用随机森林时，虽然随机森林本身可以处理多分类，但该策略更有效.
    :param link: 链接feature的网址
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

    if slgr == 'yes' and classifier != 'rf':
        # 使用单变量逻辑回归，对每一个筛选出的变量进行分析并进行ROC可视化
        min_auc = 0.5 if X.shape[1] > 100 else 0
        per_logit_reg(X.T, group_info, min_auc=min_auc, link=link)

    # 如果不是随机森林模型，则先进行feature聚类并派选代表
    if classifier != 'rf':
        represents = group_collinear_vars(
            X.T, corr_cutoff=corr_cutoff, method='spearman'
        )
        X = X[[x[0] for x in represents]]

    # step2: optimize model parameter
    print('网格搜索最优模型参数')
    param_optimized_model = multi_grid_search_model(
        X, y, classifier=classifier, ovr=ovr,
        max_iter=grid_search_num, test_size=test_size, cv=cv
    )
    # step3: feature selection
    if not no_feature_selection and classifier == 'rf':
        # 目前borutapy可能仅支持rf
        print('使用BorutaPy筛选出所有重要的特征')
        if not ovr or len(set(y))==2:
            feat_selector = BorutaPy(
                param_optimized_model,
                n_estimators='auto', verbose=0, random_state=1,
                max_iter=200, perc=percent, alpha=alpha
            )
            # find all relevant features
            feat_selector.fit(X.values, y)
            # check ranking of features
            rank = pd.DataFrame({'feature':X.columns, 'rank': feat_selector.ranking_})
            rank.sort_values(by='rank').to_csv('feature_ranking.xls', sep='\t', index=False)
            # select feature
            selected = rank[rank['rank'] == 1]['feature']
            if tentative:
                selected = rank[(rank['rank'] == 2) | (rank['rank'] == 1)]['feature']
        else:
            selected = []
            ranks = []
            for i, label in enumerate(param_optimized_model.classes_):
                feat_selector = BorutaPy(
                    param_optimized_model.estimators_[i],
                    n_estimators='auto', verbose=1, random_state=1,
                    max_iter=100, perc=percent, alpha=alpha
                )
                # find all relevant features
                feat_selector.fit(X.values, y)
                # check ranking of features
                col_name = f'{label}.rank'
                rank = pd.DataFrame({'feature': X.columns, col_name: feat_selector.ranking_})
                ranks.append(rank)
                # select feature
                this_selected = rank[rank[col_name] == 1]['feature']
                if tentative:
                    this_selected = rank[(rank[col_name] == 2) | (rank[col_name] == 1)]['feature']
                selected += [x for x in this_selected if x not in selected]
            pd.concat(ranks, axis=1).to_csv('feature_ranking.xls', sep='\t', index=False)
        X = X.loc[:, selected]
        print('final selected features:', X.columns)

        # 使用单变量逻辑回归，对每一个筛选出的变量进行分析并进行ROC可视化
        if slgr == 'yes':
            per_logit_reg(X.T, group_info, min_auc=0, link=link)

        # step4: 找共线性并举出代表
        represents = group_collinear_vars(
            X.T, corr_cutoff=corr_cutoff, method='spearman'
        )

        # step5: 使用最终筛选出来的代表性feature进行最终的模型训练和评估
        X = X[[x[0] for x in represents]]
        print('基于筛选出来的feature数据，重新建模并网格搜索最优模型参数')
        best_model = multi_grid_search_model(
            X, y, classifier=classifier, ovr=ovr, prefix='final.',
            max_iter=grid_search_num, test_size=test_size, cv=cv
        )

        print(f'The final best model(with {len(represents)} features) is',
              best_model.__str__().replace("\n", '').replace(' ', ''))
        # save final data used for building model
        data = X.T
        rep_dict = dict((x[0], x[1]) for x in represents)
        data.insert(0, 'members', [rep_dict[x] for x in data.index])
        data.to_csv('final_model_based_data.xls', sep='\t')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['run', 'group_collinear_vars', 'per_logit_reg'])
