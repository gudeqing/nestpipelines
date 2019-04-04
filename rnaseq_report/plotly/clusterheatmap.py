#! /data/users/dqgu/anaconda3/bin/python
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
__author__ = 'gudeqing'


class ClusterHeatMap(object):
    def __init__(self, data_file=None, out_name='clusterHeatMap.html',
                 sample_cluster_method='complete', sample_distance_metric="correlation",
                 gene_cluster_method='average', gene_distance_metric="euclidean",
                 cluster_gene=False, cluster_sample=False,
                 show_gene_label=False, hide_sample_label=False, show_legend=True,
                 only_sample_dendrogram=False,
                 only_gene_dendrogram=False,
                 do_correlation_cluster=False, corr_method='pearson',
                 sample_cluster_num=1, gene_cluster_num=1,
                 sample_group=None, log_base=2, log_additive=1.0, zscore_before_cluster=False,
                 lower_exp_cutoff=0.5, pass_lower_exp_num=None,
                 row_sum_cutoff=1, cv_cutoff=0., target_cols=None, target_rows=None, gene_annot=None,
                 width=1000, height=800, group_color=None, sort_cluster_by='distance',
                 gene_label_size=6, sample_label_size=10, sample_label_angle=45, outlier_k=3.0,
                 color_scale='YlGnBu', preprocess_data_func=None, transpose_data=False,
                 left_dendrogram_width=0.15, top_dendrogram_height=0.15):
        """
        cluster / correlation cluster for gene expression
        For cluster method and metric option, please refer scipy.cluster.hierarchy.linkage
        :param data_file: data file path
        :param out_name: figure file name, path info can be included
        :param sample_cluster_method: default "single"
        :param sample_distance_metric: default "correlation", correlation actually refer to pearson corr
        :param gene_cluster_method: default "average"
        :param gene_distance_metric: default "euclidean"
        :param cluster_gene: bool value indicates if to cluster gene
        :param cluster_sample: bool value indicates if to cluster sample
        :param label_gene: bool value indicates if to display gene name
        :param only_sample_dendrogram: bool value indicates if to only draw sample cluster dendrogram
        :param only_gene_dendrogram: bool value indicates if to only draw gene cluster dendrogram
        :param do_correlation_cluster: bool value indicates if to cluster sample using "corr_method" and
            display correlation heat map
        :param corr_method: correlation method, could be {'pearson', 'kendall', 'spearman'}, they are from pandas.corr
        :param sample_cluster_num: number of sample cluster to output
        :param gene_cluster_num: number of gene cluster to output
        :param sample_group: sample group dict, {'sample_name': 'group_name', ...},
            or a file with two columns [sample_name, group_name]
        :param log_base: transform data using log, value could be one of {2, 10, 1, None}, 1 means no log transformation
        :param log_additive: a small value added before doing log transformation for data
        :param zscore_before_cluster: bool indicates if to do zscore normalization, default: False.
            No effect if "do_correlation_cluster" is True
        :param lower_exp_cutoff: gene expression lower cutoff, combined with pass_lower_exp_num
        :param pass_lower_exp_num: gene with expression N times smaller than "lower_exp_cutoff" will be filtered
        :param row_sum_cutoff: gene with sum of expression lower than this cutoff will be filtered, default 1
        :param cv_cutoff: genes with cv (=mean/std) higher than will be retained
        :param target_cols: target columns to extract from data file for analysis
        :param target_rows: target rows to extract from data file for analysis
        :param gene_annot: gene annotation file, two columns, gene_id \t gene symbol
        :param width: figure width
        :param height: figure height
        :param group_color: group color dict, {'group': 'color', ...},
            or a file with two columns [group, color]
        :param sort_cluster_by: sort cluster by distance or count
        :param gene_label_size: int, gen label size, default 6
        :param outlier_k: k value for determine outlier, max color value = q3+(q3-q1)*k
        :param color_scale: pallete for heat map, refer to plotly,
            ['Blackbody', 'Bluered', 'Blues', 'Earth', 'Electric',
            'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland',
            'Rainbow', 'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd']
        :param preprocess_data_func: function provided for data filtering and transformation. Default None
        :param transpose_data: transpose raw data before future analysis
        :param left_dendrogram_width: left/sample dendrogram width, default 0.15, range(0, 1)
        :param top_dendrogram_height: top/gene dendrogram height, default 0.15, range(0, 1)
        """

        self.scm = sample_cluster_method
        self.sdm = sample_distance_metric
        self.gcm = gene_cluster_method
        self.gdm = gene_distance_metric
        self.scn = sample_cluster_num
        self.gcn = gene_cluster_num
        self.group_dict = sample_group
        self.group_color = group_color
        self.transpose_data = transpose_data
        self.sort_cluster_by = sort_cluster_by
        self.outlier_k = outlier_k
        self.gene_annot = dict(x.strip().split('\t')[:2] for x in open(gene_annot)) if gene_annot else gene_annot
        self.target_cols = [x.strip().split()[0] for x in open(target_cols)] if target_cols else None
        self.target_rows = [x.strip().split()[0] for x in open(target_rows)] if target_rows else None
        self.show_legend = show_legend
        if isinstance(sample_group, str):
            if not os.path.exists(sample_group):
                raise Exception('sample group file is not existed')
            with open(sample_group) as f:
                self.group_dict = dict(line.strip().split('\t')[:2] for line in f)
        if isinstance(group_color, str):
            if not os.path.exists(group_color):
                raise Exception('sample color file is not existed')
            with open(group_color) as f:
                self.group_color = dict(line.strip().split('\t')[:2] for line in f)
        self.gene_label_size = gene_label_size
        self.sample_label_size = sample_label_size
        self.sample_label_angle = sample_label_angle
        self.do_correlation_cluster = do_correlation_cluster
        self.ordered_genes = None
        self.ordered_samples = None
        self.cluster_gene = cluster_gene
        self.cluster_sample = cluster_sample
        self.only_sample_dendrogram = only_sample_dendrogram
        self.only_gene_dendrogram = only_gene_dendrogram
        self.label_gene = show_gene_label
        self.label_sample = not hide_sample_label
        self.height = height
        self.width = width
        self.colorscale = color_scale
        self.logbase = log_base
        self.log_additive = log_additive
        self.lower_exp_cutoff = lower_exp_cutoff
        self.pass_lower_exp_num = int(pass_lower_exp_num) if pass_lower_exp_num is not None else None
        self.row_sum_cutoff = row_sum_cutoff
        self.cv_cutoff = cv_cutoff
        self.zscore_before_cluster = zscore_before_cluster

        self.out_name = out_name
        outdir = os.path.dirname(out_name)
        self.outdir = outdir if outdir else os.getcwd()
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        if data_file:
            self.data_file = data_file
            if callable(preprocess_data_func):
                self.data = preprocess_data_func(self.data_file)
            else:
                self.data = self.process_data()
        else:
            # 不断测试发现, index如果为纯数字, 当且仅当有基因聚类的时候将不能正常显示热图, 
            # 应该是plotly的bug, 推测热图自动调整绘图的过程中, 会用到数字索引, 奇怪的很！
            print('Using random data to do test !')
            self.data = pd.DataFrame(np.random.randint(0, 20, (100, 6)),
                                     columns=list('abcdef'),
                                     index=['x'+str(x) for x in range(100)])

        if self.do_correlation_cluster:
            self.data = self.data.corr(method=corr_method)
            # we choose not use the following codes to enable more freedom of users
            # self.cluster_gene = True
            # self.cluster_sample = True
            # self.label_gene = True

        if self.cluster_gene:
            self.left_dendrogram_width = left_dendrogram_width
        else:
            self.left_dendrogram_width = 0
        if self.cluster_sample:
            self.top_dendrogram_height = top_dendrogram_height
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
        self.draw()
        if self.do_correlation_cluster:
            out_corr_file = os.path.join(outdir, 'corr.matrix.txt')
            self.data.to_csv(out_corr_file, header=True, index=True, sep='\t')

    def process_data(self):
        exp_pd = pd.read_table(self.data_file, header=0, index_col=0)
        if self.target_rows:
            exp_pd = exp_pd.loc[[x for x in self.target_rows if x in exp_pd.index], :]
        if self.target_cols:
            exp_pd = exp_pd.loc[:, [x for x in self.target_cols if x in exp_pd.columns]]
        if self.transpose_data:
            exp_pd = exp_pd.transpose()
        # exp_pd = exp_pd.applymap(lambda x: x if x <=8 else 8)
        if exp_pd.shape[0] <= 1 or exp_pd.shape[1] <= 1:
            raise Exception("Data is not enough for analysis !")
        exp_pd = exp_pd[exp_pd.sum(axis=1) > self.row_sum_cutoff]
        exp_pd = exp_pd[exp_pd.std(axis=1)/exp_pd.mean(axis=1) > self.cv_cutoff]
        pass_num_cutoff = int(exp_pd.shape[1] / 3) if self.pass_lower_exp_num is None else self.pass_lower_exp_num
        if self.lower_exp_cutoff > 0:
            pass_state = exp_pd.apply(lambda x: sum(x > self.lower_exp_cutoff), axis=1)
            exp_pd = exp_pd[pass_state >= pass_num_cutoff]
        if self.logbase == 2:
            exp_pd = np.log2(exp_pd+1)
        elif self.logbase == 10:
            exp_pd = np.log10(exp_pd+self.log_additive)
        elif not self.logbase or self.logbase == 1:
            pass
        else:
            raise Exception('log base must be one of [2, 10, 1] ')
        out_name = os.path.join(self.outdir, 'cluster.log{}.cv{}.{}outof{}over{}.preprocessed.data'.format(
            self.logbase, self.cv_cutoff, pass_num_cutoff, exp_pd.shape[1], self.lower_exp_cutoff
        ))
        exp_pd.to_csv(out_name, header=True, index=True, sep='\t')
        return exp_pd

    def heatmap_xaxis(self):
        return {
            'domain': [self.left_dendrogram_width, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'ticks': "",
            'dtick': 1,
            'anchor': 'y',
            'autorange': True,
            'scaleanchor': "x3",
            'tickangle': self.sample_label_angle,
            'tickfont': {'size': self.sample_label_size},
            'showticklabels': self.label_sample
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
            showlegend=True,
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
        out_name = os.path.join(self.outdir, 'cluster.heatmap.data')
        # plotly plot data from bottom to top, thus we have to use [::-1], or it will not match dendrogram
        heat_data = self.data.iloc[self.ordered_genes[::-1], self.ordered_samples]
        # trans gene id to gene name
        if self.gene_annot:
            heat_data.index = [self.gene_annot[x] if x in self.gene_annot else x for x in heat_data.index]
        # output data
        heat_data.to_csv(out_name, header=True, index=True, sep='\t')
        # process heat data to make color be more even
        describe = pd.Series(heat_data.values.flatten()).describe()
        upper_limit = describe["75%"] + (describe["75%"] - describe["25%"])*self.outlier_k
        upper_limit = upper_limit if upper_limit < describe['max'] else describe['max']
        lower_limit = describe["25%"] - (describe["75%"] - describe["25%"])*self.outlier_k
        lower_limit = lower_limit if lower_limit > describe['min'] else describe['min']
        heat_map = go.Heatmap(
            x=list(heat_data.columns),
            y=list(heat_data.index),
            z=heat_data.values,
            colorscale=self.colorscale,
            showlegend=False,
            xaxis='x',
            yaxis='y',
            zmin=lower_limit,
            zmax=upper_limit,
            name='',
            showscale=True,
            colorbar=dict(
                x=1+self.gene_label_size*(max(len(x) for x in heat_data.index))/self.width if self.label_gene else 1,
                xanchor='left',
                y=0,
                yanchor='bottom',
                len=0.5,
                title="log{}(X)".format(self.logbase) if self.logbase and self.logbase != 1 else ''
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
        groups = list()
        _ = [groups.append(x) for x in self.group_dict.values() if x not in groups]
        colors = self.get_color_pool(len(groups))
        group_colors = dict(zip(groups, colors))
        if self.group_color:
            for k, v in self.group_color.items():
                group_colors[k] = v

        sample_colors = dict()
        for sample in ordered_samples:
            group = self.group_dict[sample]
            sample_colors[sample] = group_colors[group]

        traces = list()
        ticks = range(5, self.data.shape[1]*10, 10)
        existed_legend = set()
        for tick, each in zip(ticks, ordered_samples):
            if each == self.group_dict[each]:
                trace_name = each
            else:
                trace_name = '{sample}({group})'.format(sample=each, group=self.group_dict[each])
            bar = go.Bar(
                name=self.group_dict[each],
                text=trace_name,
                x=[tick],
                y=[2],
                showlegend=True if sample_colors[each] not in existed_legend else False,
                xaxis='x4',
                yaxis='y4',
                marker=dict(color=sample_colors[each]),
                hoverinfo="text",
                hoverlabel=dict(namelength=-1)
            )
            traces.append(bar)
            existed_legend.add(sample_colors[each])
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
            count_sort=True if self.sort_cluster_by == 'count' else False,
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
                yaxis="y3",
                showlegend=False,
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
            count_sort=True if self.sort_cluster_by == 'count' else False,
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
                yaxis="y2",
                showlegend=False
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
        :param n_clusters: int, optional, default: 10. The number of clusters to generate.
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
            print("n_clusters is bigger than sample number",
                  "it will be forced to be sample number !")
            n_clusters = exp_pd.shape[0]

        if self.do_correlation_cluster:
            self.logbase = None
            condensed_distance = squareform(1 - exp_pd)
            z = hclust.linkage(condensed_distance)
        else:
            if self.zscore_before_cluster:
                from sklearn import preprocessing
                exp_pd = exp_pd.apply(preprocessing.scale, axis=0)
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
        if n_clusters <= 1:
            return z, None
        # write out subcluster
        labels = exp_pd.index
        subcluster = self.get_subcluster(z, labels, num=n_clusters)
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
        :param z hclust.linkage result
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
        # https://plot.ly/ipython-notebooks/color-scales/
        import colorlover
        if n <= 12:
            return colorlover.scales['12']['qual']['Paired']

        from colorsys import hls_to_rgb
        color_pool = []
        for i in np.arange(0., 360., 360. / n):
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
        if self.group_dict:
            self.layout['legend'] = dict(x=traces[-1]['colorbar']['x'])
        if not self.show_legend:
            self.layout['showlegend'] = False
        fig = go.Figure(data=traces, layout=self.layout)
        plt(fig, filename=self.out_name, auto_open=False)


if __name__ == '__main__':
    def introduce_command(func):
        import argparse
        import inspect
        import json
        import time
        if isinstance(func, type):
            description = func.__init__.__doc__
        else:
            description = func.__doc__
        parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
        func_args = inspect.getfullargspec(func)
        arg_names = func_args.args
        arg_defaults = func_args.defaults
        arg_defaults = ['None']*(len(arg_names) - len(arg_defaults)) + list(arg_defaults)
        for arg, value in zip(arg_names, arg_defaults):
            if arg == 'self':
                continue
            if value == 'None':
                parser.add_argument('-'+arg, required=True, metavar=arg)
            elif value is True:
                parser.add_argument('--'+arg, action="store_false", help='bool, default: True')
            elif value is False:
                parser.add_argument('--'+arg, action="store_true", help='bool, default: False')
            elif value is None:
                parser.add_argument('-' + arg, default=value, metavar='Default:' + str(value), )
            else:
                parser.add_argument('-' + arg, default=value, type=type(value), metavar='Default:' + str(value), )
        if func_args.varargs is not None:
            print("warning: *varargs is not supported, and will be neglected! ")
        if func_args.varkw is not None:
            print("warning: **keywords args is not supported, and will be neglected! ")
        args = parser.parse_args().__dict__
        with open("Argument_detail.json", 'w') as f:
            json.dump(args, f, indent=2, sort_keys=True)
        start = time.time()
        func(**args)
        print("total time: {}s".format(time.time() - start))

    introduce_command(ClusterHeatMap)
