import os
import plotly.graph_objs as go
from plotly.offline import plot as plt
import pandas as pd
import numpy as np
import scipy as scp
from scipy.cluster import hierarchy as sch
import fastcluster as hclust
import colorlover


class ClusterHeatMap():
    def __init__(self, data_file, method='average', metric="correlation",
                 out_name='clusterHeatMap.html'):
        self.data_file = data_file
        self.method = method
        self.metric = metric
        self.out_name = out_name
        outdir = os.path.dirname(out_name)
        self.outdir = outdir if outdir else os.getcwd()
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        self.data = self.process_data()
        self.ordered_genes = None
        self.ordered_samples = None
        self.layout = self.all_layout()


    def process_data(self):
        from sklearn import decomposition, preprocessing
        exp_pd = pd.read_table(self.data_file, header=0, index_col=0)
        exp_pd = exp_pd[exp_pd.sum(axis=1) > 0]
        if exp_pd.shape[0] <= 1 or exp_pd.shape[1] <=1:
            raise Exception("Data is not enough for analysis !")
        exp_pd = exp_pd[exp_pd.std(axis=1)/exp_pd.mean(axis=1)>0.5]
        exp_pd = np.log(exp_pd+1)
        exp_pd = exp_pd.apply(preprocessing.scale, axis=0)
        exp_pd = exp_pd.iloc[:300, :]
        return exp_pd

    def heatmap_xaxis(self):
        return {
            'domain': [.151, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'ticks':"",
            'anchor': 'y',
            'autorange': True,
            'constrain': 'domain'

        }

    def heatmap_yaxis(self):
        return {
            'domain': [0, 0.85],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': "",
            'anchor': 'x',
            'autorange': True,
            'constrain': 'domain',
            'scaleanchor': "y2",
        }

    def left_dendrogam_xaxis(self):
        return {
            'domain': [0, 0.15],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': "",
            'anchor': 'y2',
        }

    def left_dendrogam_yaxis(self):
        return {
            'domain': [0, 0.85],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': "",
            'scaleanchor': "y",
            'anchor': 'x2',
            'constrain': 'domain',
            'range': (-self.data.shape[0]*10, 1)
        }

    def top_dendrogram_xaxis(self):
        return {
            'domain': [0.15+0.85/2/self.data.shape[1], 1-0.85/2/self.data.shape[1]],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': "",
            'anchor': 'y3',
            'scaleanchor': 'x'
        }

    def top_dendrogram_yaxis(self):
        return {
            'domain': [.85, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': "",
            'anchor': 'x3'
        }

    def all_layout(self):
        return go.Layout(
            width=800,
            height=800,
            showlegend=False,
            hovermode='closest',
            xaxis=self.heatmap_xaxis(),
            yaxis=self.heatmap_yaxis(),
            xaxis2=self.left_dendrogam_xaxis(),
            yaxis2=self.left_dendrogam_yaxis(),
            xaxis3=self.top_dendrogram_xaxis(),
            yaxis3=self.top_dendrogram_yaxis(),
        )

    def heatmap_trace(self):
        if not self.ordered_samples:
            self.ordered_samples = range(self.data.shape[1])
        if not self.ordered_genes:
            self.ordered_genes = range(self.data.shape[0])
        heat_data = self.data.iloc[self.ordered_genes, self.ordered_samples]
        heat_map = go.Heatmap(
            x=heat_data.columns,
            y=heat_data.index,
            z=heat_data.values.tolist(),
            colorscale='YlGnBu',
            showlegend=False,
            xaxis='x',
            yaxis='y',
            name=''
        )
        return [heat_map]

    def top_dendrogram_traces(self):
        exp_pd = self.data.transpose()
        z, subcluster = self.hcluster(
            exp_pd,
            method=self.method,
            metric=self.metric,
            transpose=False,
            n_clusters=2,
            output=self.outdir,
            prefix='sample.'
        )
        self.set_link_color_palette()
        results = sch.dendrogram(
            z,
            orientation="top",
            no_plot=True,
            distance_sort=True,
            above_threshold_color='black'
        )
        self.ordered_samples = list(map(int, results['ivl']))
        icoord = scp.array(results['icoord'])
        dcoord = scp.array(results['dcoord'])
        color_list = scp.array(results['color_list'])
        # print(results['color_list'])
        trace_list = []
        for i in range(len(icoord)):
            # x and y are arrays of 4 points that make up the '∩' shapes of the dendrogram tree
            hovertext_label = None
            trace = go.Scatter(
                x=icoord[i],
                y=dcoord[i],
                mode='lines',
                marker=dict(color=color_list[i]),
                text=hovertext_label,
                hoverinfo='text',
                xaxis="x3",
                yaxis="y3"
            )
            trace_list.append(trace)
        return trace_list

    def left_dendrogram_traces(self):
        exp_pd = self.data
        z, subcluster = self.hcluster(
            exp_pd,
            method=self.method,
            metric=self.metric,
            transpose=False,
            n_clusters=2,
            output=self.outdir,
            prefix='gene.',
        )
        self.set_link_color_palette()
        results = sch.dendrogram(
            z,
            orientation="left",
            no_plot=True,
            distance_sort=True,
            above_threshold_color='black'
        )
        self.ordered_genes = list(map(int, results['ivl']))
        icoord = scp.array(results['dcoord'])*(-1)
        dcoord = scp.array(results['icoord'])*(-1)
        color_list = scp.array(results['color_list'])
        trace_list = []
        # print(min(dcoord.flatten()))
        for i in range(len(icoord)):
            # x and y are arrays of 4 points that make up the '∩' shapes of the dendrogram tree
            hovertext_label = None
            trace = go.Scatter(
                x=icoord[i],
                y=dcoord[i],
                mode='lines',
                marker=dict(color=color_list[i]),
                text=hovertext_label,
                hoverinfo='text',
                xaxis="x2",
                yaxis="y2"
            )
            trace_list.append(trace)
        return trace_list

    def hcluster(self, exp_pd, transpose=False, n_clusters=10, method='average',
                 metric='correlation', output=None, prefix=''):
        """
        'fastcluster' was used, http://www.danifold.net/fastcluster.html?section=3.
        scipy.cluster.hierarchy.linkage could also do the same thing but slower.
        However, the documentation of scipy.cluster.hierarchy.linkage is pretty good.
        >> https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html

        :param exp_pd: pandas DataFrame, expression matrix
        :param transpose: if to transpose the expression matrix
        :param n_clusters: int, optional, default: 8. The number of clusters to generate.
        :param method: methods for calculating the distance between the newly formed clusters.
            Choices: ['single', 'average', 'weighted', 'centroid', 'complete', 'median', 'ward']
        :param metric: The distance metric to use.
            Choices: ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation',
            'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski',
            'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao',
            'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule', ]
        :param output: output directory of subcluster information
        :param prefix: outfile prefix letters
        :return: tree -> tuple(tree_str, tree_list),
                 subcluster -> {0:[s1,s2], 1:[s4,s5], ...}
                 z -> cluster result of hclust.linkage
        """
        if transpose:
            exp_pd = exp_pd.transpose()
        if n_clusters > exp_pd.shape[0]:
            print("n_clusters is bigger than sample number!!")
            n_clusters = exp_pd.shape[0]
        try:
            z = hclust.linkage(exp_pd, method=method, metric=metric)
        except FloatingPointError as e:
            print("fastcluster failed as : {}".format(e))
            print('it seems that (at least) one of the vectors you want to cluster is all zeros, '
                  'so when it tries to compute the cosine distances to it there is a division by zero,'
                  ' hence nan is stored in your distance array, and that leads to your error.')
            print('Anyway, we will remove the special rows for you now.')
            check = exp_pd[exp_pd.sum(axis=1) == 0]
            if check.shape[0] >= 1:
                print('Actually, we detected that some genes have zero expression across all samples.')
                print('such as {} has zero expression across all sample'.format(check.index[0]))
                exp_pd = exp_pd[exp_pd.sum(axis=1) > 0]
            try:
                z = hclust.linkage(exp_pd, method=method, metric=metric)
            except:
                print("enhen? fastcluster failed again, we will try scipy.cluster.hierarchy")
                z = sch.linkage(exp_pd, method=method, metric=metric)
        except:
            print("fastcluster failed, we will try scipy.cluster.hierarchy")
            z = sch.linkage(exp_pd, method=method, metric=metric)
        labels = exp_pd.index
        subcluster = self.get_subcluster(z, labels, num=n_clusters)
        # write out subcluster
        if output is not None:
            if not os.path.exists(output):
                os.mkdir(output)
        else:
            output = os.getcwd()
        for k, v in subcluster.items():
            out_dir = os.path.join(output, prefix + 'subcluster_{}_{}.xls'.format(k, len(v)))
            sub = exp_pd.loc[v, :]
            if transpose:
                sub = sub.transpose()
            sub.to_csv(out_dir, sep='\t', header=True, index=True)
        # write out cluster result z
        out_dir = os.path.join(output, prefix + "linkage_result")
        pd.DataFrame(z).to_csv(out_dir, sep='\t')
        return z, subcluster

    @staticmethod
    def get_subcluster(z, labels, num=2):
        """
        get leave ids for each sub cluster
        :param zclust hclust.linkage result
        :param num: the number of sub-clusters specified.
        :param labels: leaf label from DataFrame.columns
        :return: dict with list of samples as element.
        """
        cluster = sch.cut_tree(z, num)
        tmp_pd = pd.DataFrame(cluster)
        tmp_pd['label'] = labels
        result = tmp_pd.groupby(0).groups
        subcluster = dict()
        for k in result:
            subcluster[k] = list(labels[result[k]])
        return subcluster

    @staticmethod
    def set_link_color_palette():
        colors = colorlover.scales['12']['qual']['Paired']
        sch.set_link_color_palette(colors)

    def draw(self):
        traces = list()
        traces += self.top_dendrogram_traces()
        traces += self.left_dendrogram_traces()
        traces += self.heatmap_trace()
        fig = go.Figure(data=traces, layout=self.layout)
        plt(fig, filename=self.out_name, auto_open=False)


if __name__ == '__main__':
    import sys
    p = ClusterHeatMap(sys.argv[1])
    p.draw()






