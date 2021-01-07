import matplotlib
matplotlib.use('agg')
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster import hierarchy
import fastcluster as hclust

""""
共线性问题会导致逻辑回归建模偏差甚至错误，本质原因是因为违背了变量间相互独立的假设。
虽然随机森林模型不惧怕共线性，从逻辑上不会影响最终的模型预测效果，但共线性问题也会使得随机森林模型对feature重要程度的打分不准。
所以，当我们研究哪个feature更重要时，或者要进行feature解释时，共线性问题必须要合理解决，毕竟我们不应该随机的进行feature选择。

0.假设输入的是表达矩阵，行为基因，列为样本，基因成为feature或特征
1. 使用Boruta 筛选出所有重要的特征feature，共线性问题对于该方法无影响，该方法可以筛选出所有同等重要的feature
2. 进行VIF计算，VIF衡量某个变量和剩余变量之间的共线性，该值越大，共线性越大，一般大于5时认为存在潜在共线性，大于10时，存在严重共线性。
   如果两个变量高度相关，那么他们的VIF应该都很高。
3. 进行基因之间相关系数计算，以spearman相关性为主，然后进行层级聚类，聚类到一起的一般共线性嫌疑较大，取其中一个作为代表往往不影响聚类结果。
   需要注意一点：相关性不具备完全传递性。这是一个比较棘手的问题。 例如，当我们知道A和B相关性很高，B和C的相关性很高，但是A和C之间可能几乎没有相关性。
   如果我们试图用A同时代表B和C是存在问题，用B代表A和C可能稍微好一些，但是要保留更多信息的办法可能是保留A和C
   
4. 去共线性策略思考：
    要去除共线性，不应该粗暴简单的把VIF较大的变量删除，虽然这样做也可以去除共线性，但留下来的变量未必更好。
    应该先想办法找出那些基因之间存在共线性，再考虑如何取舍变量。
    a. 先算出每一个基因的feature重要程度排名，同时也算出每一个基因的VIF，同时算出每一个基因和其他高度相关的基因
    b. 观察a的结果，思考如何找到共线性团体
5. 由于不同批次实验的rnaseq数据不具备可比较性，蛋白组数据同样面临相同的情况。
   基于这种实验得到的表达矩阵需要做什么处理才能构建一个可用于实际预测的模型？
   现在假设我们得到的表达矩阵在横向和纵向方面均具备可比较性（deseq2的标准化矩阵在纵向也即不同基因间可比较性较差，但是TPM具备纵向可比较性）
   基于当前分析经验，我们知道pca分析或聚类分析时，通常需要先按行进行scale（scale是一种伸缩变换，不是常规的标准化）才能得到的较好区分效果。
   这种scale消除了高表达基因和低表达基因之间的量纲差异，即抛弃绝对的高表达和低表达的信息，仅保留相对高低的信息。
   我们知道，我们会说rnaseq定量是相对定量，这里的相对是指不同基因之间的相对性，通常针对TPM这种指标而言，按行scale就消除了这种相对性。
   由于不同实验数据不具备可比较性，即使基于具备纵向或横向可比较性的表达矩阵，构建模型时仍需要考虑是否进行scale？
   a. 如果进行横向scale，这似乎会导致你的模型无法对新的实验数据无法进行预测，毕竟你没有办法对一个新的样本按行scale。
   可挽救的办法，考虑把这个样本混进原来的训练集后进行标准化，但这时你要考虑如对于rnaseq数据：
   从count出发，如果你加入一个新样本到训练集，标准化后，整个表的矩阵发生或多或少的改变，那你要考虑是否需要重新训练模型吗？
   b. 如果进行纵向scale，此时即每个样本各自进行scale标准化，训练好模型后，对新的数据直接纵向scale后即可预测。
   但要考虑的问题是：纵向scale后的数据适合用来进行模型训练吗？对于随机森林模型，
   
"""


def get_vif(corr_df):
    vif = np.linalg.inv(corr_df.values).diagonal()
    vif_df = pd.DataFrame(vif, index=corr_df.index, columns=['VIF'])
    return vif_df.sort_values(by='VIF', ascending=False)


def corr_cluster(corr, distance_cutoff=0, link_method='ward'):
    # 把相关系数取绝对值后进行计算，如此正相关和负相关都会聚类到一起
    condensed_distance = squareform(1 - corr.abs())
    corr_linkage = hclust.linkage(condensed_distance, method=link_method)
    # plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
    dendro = hierarchy.dendrogram(corr_linkage, ax=ax1, leaf_rotation=90)
    ax1.set_xticklabels(corr.iloc[dendro['leaves']].index)
    ax2.imshow(corr.iloc[dendro['leaves'], :].iloc[:, dendro['leaves']])
    fig.tight_layout()
    plt.savefig('spearman_corr.pdf')

    # 把距离小于distance—cutoff的叶子归为一组
    distance_desc = pd.DataFrame(corr_linkage)[2].describe()
    if distance_cutoff <= 0:
        lq = distance_desc['25%']
        print('25% distance cutoff:', lq)
    else:
        lq = distance_cutoff
    cluster_ids = hierarchy.fcluster(corr_linkage, lq, criterion='distance')
    distance_group = dict()
    ids = corr.columns
    for ind, clst_id in enumerate(cluster_ids):
        distance_group.setdefault(clst_id, list())
        distance_group[clst_id].append(ids[ind])

    for k, v in distance_group.items():
        if len(v) > 1:
            inter_corr = corr.loc[v, v]
            vif_df = get_vif(inter_corr)
            mean_corr = (inter_corr.values.flatten().sum() - len(v)) / (len(v) ** 2 - len(v))
            merged = vif_df.join(inter_corr).round(4)
            merged.index.name = f'mean_corr:{mean_corr:.4f}'
            print(merged)


def corr_pair(corr, cutoff=0.7):
    # 计算高相关性的基因
    high_corr_group = set()
    for ind, row in corr.abs().iterrows():
        group = set(row[row >= cutoff].index)
        if len(group) > 1:
            high_corr_group.add(tuple(sorted(group)))
    # merge group
    merged = set()
    for each in high_corr_group:
        tmp = set(each)
        _ = [tmp.update(set(x)) for x in high_corr_group if set(x) & tmp]
        merged.add(tuple(sorted(tmp)))
    # save
    with open('HighCorrGroup.list', 'w') as f:
        for each in merged:
            f.write(';'.join(each) + '\n')
    return merged


def view_relation(data, target_rows=None, target_cols=None, corr_method='spearman'):
    df = pd.read_csv(data, header=0, index_col=0, sep=None, engine='python')
    if target_rows is not None:
        rows = [x.strip().split()[0] for x in open(target_rows)]
        rows = [x for x in rows if x in df.index]
        df = df.loc[rows]
    if target_cols is not None:
        cols = [x.strip().split()[0] for x in open(target_cols)]
        cols = [x for x in cols if x in df.columns]
        df = df[cols]
    corr = df.T.corr(method=corr_method)
    corr_cluster(corr)
    print(get_vif(corr))
    corr_pair(corr)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
