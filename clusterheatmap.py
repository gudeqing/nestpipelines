import os
import plotly.graph_objs as go
from plotly.offline import plot as plt
import pandas as pd
import numpy as np
import scipy as scp
from scipy.spatial.distance import squareform
from scipy.cluster import hierarchy as sch
import fastcluster as hclust
import colorlover


class ClusterHeatMap():
    def __init__(self, data_file=None, out_name='clusterHeatMap.html',
                 sample_cluster_method='single', sample_distance_metric="correlation",
                 gene_cluster_method='average', gene_distance_metric="euclidean",
                 cluster_gene=True, cluster_sample=True, label_gene=False,
                 only_sample_dendrogram=False,
                 only_gene_dendrogram=False,
                 do_correlation_cluster=False, corr_method='pearson',
                 sample_cluster_num=2, gene_cluster_num=10,
                 sample_group=None,
                 width=800, height=800, gene_label_size=6,
                 color_scale='YlGnBu'):

        self.scm = sample_cluster_method
        self.sdm = sample_distance_metric
        self.gcm = gene_cluster_method
        self.gdm = gene_distance_metric
        self.scn = sample_cluster_num
        self.gcn = gene_cluster_num
        self.group_dict = sample_group
        self.gene_label_size = gene_label_size
        self.do_correlation_cluster = do_correlation_cluster
        self.ordered_genes = None
        self.ordered_samples = None
        self.cluster_gene = cluster_gene
        self.cluster_sample = cluster_sample
        self.only_sample_dendrogram = only_sample_dendrogram
        self.only_gene_dendrogram = only_gene_dendrogram
        self.label_gene = label_gene
        self.height=height
        self.width=width
        self.colorscale = color_scale

        self.out_name = out_name
        outdir = os.path.dirname(out_name)
        self.outdir = outdir if outdir else os.getcwd()
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        if data_file:
            self.data_file = data_file
            self.data = self.process_data()
        else:
            # 不断测试发现，index如果为纯数字，当且仅当有基因聚类的时候将不能正常显示热图，
            # 应该是plotly的bug，推测热图自动调整绘图的过程中，会用到数字索引，奇怪的很！
            self.data = pd.DataFrame(np.random.randint(0, 20, (100, 6)),
                                     columns=list('abcdef'),
                                     index=['x'+str(x) for x in range(100)])
            self.data.to_csv('tmp.xls', header=True, index=True, sep='\t')
            self.data_file = 'tmp.xls'
            self.data = self.process_data()

        if do_correlation_cluster:
            self.data = self.data.corr(method=corr_method)

        if cluster_gene:
            self.left_dendrogram_width = 0.15
        else:
            self.left_dendrogram_width = 0
        if cluster_sample:
            self.top_dendrogram_height = 0.15
        else:
            self.top_dendrogram_height = 0

        if self.only_sample_dendrogram:
            self.left_dendrogram_width = 0
            self.top_dendrogram_height = 1
        if self.only_gene_dendrogram:
            self.left_dendrogram_width = 1
            self.top_dendrogram_height = 0

        if sample_group:
            self.group_bar_height = 0.025
            if self.top_dendrogram_height > 0:
                self.top_dendrogram_height = self.top_dendrogram_height - self.group_bar_height
        else:
            self.group_bar_height = 0

        self.layout = self.all_layout()

    def process_data(self):
        from sklearn import preprocessing
        exp_pd = pd.read_table(self.data_file, header=0, index_col=0)
        exp_pd = exp_pd[exp_pd.sum(axis=1) > 0]
        if exp_pd.shape[0] <= 1 or exp_pd.shape[1] <= 1:
            raise Exception("Data is not enough for analysis !")
        exp_pd = exp_pd[exp_pd.std(axis=1)/exp_pd.mean(axis=1) > 0.1]
        exp_pd = np.log(exp_pd+1)
        # exp_pd = exp_pd.apply(preprocessing.scale, axis=0)
        exp_pd = exp_pd.iloc[:100, :]
        return exp_pd

    def heatmap_xaxis(self):
        return {
            'domain': [self.left_dendrogram_width, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'ticks': "",
            'anchor': 'y',
            'autorange': True,
            'scaleanchor': "x3",
        }

    def heatmap_yaxis(self):
        return {
            'domain': [0, 1 - self.top_dendrogram_height - self.group_bar_height],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': self.label_gene,
            'side': 'right',
            'tickfont': dict(size=self.gene_label_size),
            'dtick': 1,
            'ticks': "",
            'anchor': 'x',
            'scaleanchor': "y2",
            'layer': 'above traces'
        }

    def left_dendrogam_xaxis(self):
        return {
            'domain': [0, self.left_dendrogram_width],
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
            'domain': [0, 1 - self.top_dendrogram_height - self.group_bar_height],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': "",
            'tickfont': dict(size=self.gene_label_size),
            'scaleanchor': "y",
            'anchor': 'x2',
            'range': (-self.data.shape[0]*10, 1)
        }

    def top_dendrogram_xaxis(self):
        return {
            'domain': [self.left_dendrogram_width, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': "",
            'anchor': 'y3',
            'scaleanchor': 'x',
            'range': (0, self.data.shape[1]*10)
        }

    def top_dendrogram_yaxis(self):
        return {
            'domain': [1 - self.top_dendrogram_height, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': "",
            'anchor': 'x3'
        }

    def group_bar_xaxis(self):
        return {
            'domain': [self.left_dendrogram_width, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': "",
            'anchor': 'y4',
            'scaleanchor': 'x',
        }

    def group_bar_yaxis(self):
        return {
            'domain': [1 - self.top_dendrogram_height - self.group_bar_height,
                       1 - self.top_dendrogram_height],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': "",
            'anchor': 'x4'
        }

    def all_layout(self):
        if self.only_sample_dendrogram:
            return go.Layout(
                width=self.width,
                height=self.height,
                autosize=True,
                showlegend=False,
                hovermode='closest',
                xaxis3=self.top_dendrogram_xaxis(),
                yaxis3=self.top_dendrogram_yaxis(),
                xaxis4=self.group_bar_xaxis(),
                yaxis4=self.group_bar_yaxis()
            )

        if self.only_gene_dendrogram:
            return go.Layout(
                width=self.width,
                height=self.height,
                showlegend=False,
                hovermode='closest',
                xaxis2=self.left_dendrogam_xaxis(),
                yaxis2=self.left_dendrogam_yaxis(),
            )

        return go.Layout(
            width=self.width,
            height=self.height,
            showlegend=False,
            hovermode='closest',
            xaxis=self.heatmap_xaxis(),
            yaxis=self.heatmap_yaxis(),
            xaxis2=self.left_dendrogam_xaxis(),
            yaxis2=self.left_dendrogam_yaxis(),
            xaxis3=self.top_dendrogram_xaxis(),
            yaxis3=self.top_dendrogram_yaxis(),
            xaxis4=self.group_bar_xaxis(),
            yaxis4=self.group_bar_yaxis()
        )

    def heatmap_trace(self):
        if not self.ordered_samples:
            self.ordered_samples = range(self.data.shape[1])
        if not self.ordered_genes:
            self.ordered_genes = range(self.data.shape[0])
        heat_data = self.data.iloc[self.ordered_genes, self.ordered_samples]
        heat_map = go.Heatmap(
            x=list(heat_data.columns),
            y=list(heat_data.index),
            z=heat_data.values,
            colorscale=self.colorscale,
            showlegend=False,
            xaxis='x',
            yaxis='y',
            name='',
            showscale=True,
            colorbar=dict(
                x=1+self.gene_label_size*(max(len(x) for x in heat_data.index))/self.width if self.label_gene else 1,
                xanchor='left',
                y=1,
                yanchor='top',
                len=0.5
            )
        )
        return [heat_map]

    def group_bar_traces(self):
        if not self.ordered_samples:
            self.ordered_samples = range(self.data.shape[1])
        ordered_samples = self.data.columns[self.ordered_samples]
        for sample in ordered_samples:
            if sample not in self.group_dict:
                self.group_dict[sample] = sample
        groups = list(set(self.group_dict.values()))
        colors = self.get_color_pool(len(groups))
        group_colors = dict(zip(groups, colors))

        sample_colors = dict()
        for sample in ordered_samples:
            group = self.group_dict[sample]
            sample_colors[sample] = group_colors[group]

        traces = list()
        ticks = range(5, self.data.shape[1]*10, 10)
        for tick, each in zip(ticks, ordered_samples):
            bar = go.Bar(
                name=each,
                x=[tick],
                y=[2],
                showlegend=False,
                xaxis='x4',
                yaxis='y4',
                marker=dict(color=sample_colors[each]),
                hoverinfo="name",
            )
            traces.append(bar)
        return traces

    def top_dendrogram_traces(self):
        exp_pd = self.data.transpose()
        z, subcluster = self.hcluster(
            exp_pd,
            method=self.scm,
            metric=self.sdm,
            transpose=False,
            n_clusters=self.scn,
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
        #
        if self.only_sample_dendrogram:
            self.layout['xaxis3']['showticklabels'] = True
            tick_values = list(range(5, exp_pd.shape[0]*10, 10))
            self.layout['xaxis3']['tickvals'] = tick_values
            self.layout['xaxis3']['ticktext'] = exp_pd.iloc[self.ordered_samples].index

        return trace_list

    def left_dendrogram_traces(self):
        exp_pd = self.data
        z, subcluster = self.hcluster(
            exp_pd,
            method=self.gcm,
            metric=self.gdm,
            transpose=False,
            n_clusters=self.gcn,
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

        if self.only_gene_dendrogram:
            self.layout['yaxis2']['showticklabels'] = True
            tick_values = list(range(5, (exp_pd.shape[0]+1)*10, 10))
            self.layout['yaxis2']['tickvals'] = [-1*x for x in tick_values]
            self.layout['yaxis2']['side'] = 'right'
            self.layout['yaxis2']['ticktext'] = exp_pd.iloc[self.ordered_genes].index

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
        if self.do_correlation_cluster:
            condensed_distance = squareform(1 - exp_pd)
            z = hclust.linkage(condensed_distance)
        else:
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

    @staticmethod
    def get_color_pool(n):
        import colorlover
        if n <= 12:
            return colorlover.scales['12']['qual']['Paired']

        from colorsys import hls_to_rgb
        color_pool = []
        for i in np.arange(60., 360., 360. / n):
            hue = i / 300.
            rand_num = np.random.random_sample()
            lightness = (50 + rand_num * 10) / 100.
            saturation = (90 + rand_num * 10) / 100.
            rgb = hls_to_rgb(hue, lightness, saturation)
            color_pool.append(tuple([int(x * 255) for x in rgb]))
        return colorlover.to_rgb(color_pool)

    def draw(self):
        traces = list()
        if self.only_sample_dendrogram:
            traces += self.top_dendrogram_traces()
            if self.group_dict:
                traces += self.group_bar_traces()
            fig = go.Figure(data=traces, layout=self.layout)
            plt(fig, filename=self.out_name, auto_open=False)
            return

        if self.only_gene_dendrogram:
            traces += self.left_dendrogram_traces()
            fig = go.Figure(data=traces, layout=self.layout)
            plt(fig, filename=self.out_name, auto_open=False)
            return

        if self.cluster_sample:
            traces += self.top_dendrogram_traces()

        if self.cluster_gene:
            traces += self.left_dendrogram_traces()

        if self.group_dict:
            traces += self.group_bar_traces()

        traces += self.heatmap_trace()

        fig = go.Figure(data=traces, layout=self.layout)
        plt(fig, filename=self.out_name, auto_open=False)


if __name__ == '__main__':
    import sys
    data_file = None
    if len(sys.argv) >= 2:
        data_file = sys.argv[1]
    p = ClusterHeatMap(
        data_file=data_file,
        cluster_sample=True,
        cluster_gene=True,
        gene_distance_metric="correlation",
        only_gene_dendrogram=False,
        do_correlation_cluster=False,
        label_gene=True,
        sample_group={'C180026R2L2': 'C', 'C180027R1L2': 'C', 'C180028R1L2': 'C'}
    )
    p.draw()

