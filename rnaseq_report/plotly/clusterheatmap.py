#! /data/users/dqgu/anaconda3/bin/python
import os
import plotly.graph_objs as go
from plotly.offline import plot as plt
import pandas as pd
import numpy as np
import scipy as scp
from scipy.stats import zscore
from scipy.spatial.distance import squareform
from scipy.cluster import hierarchy as sch
import fastcluster as hclust
import colorlover
__author__ = 'gudeqing'
__version_ = '3.6.7'


class ClusterHeatMap(object):
    def __init__(self, data_file=None, out_prefix='clusterHeatMap',
                 sample_cluster_method='average', sample_distance_metric="euclidean",
                 gene_cluster_method='average', gene_distance_metric="euclidean",
                 cluster_gene=False, cluster_sample=False,
                 show_gene_label=False, hide_sample_label=False, hide_legend=False,
                 only_sample_dendrogram=False,
                 only_gene_dendrogram=False,
                 sample_corr_as_heatmap=False, gene_corr_as_heatmap=False, corr_method='spearman',
                 sample_cluster_num=1, gene_cluster_num=1,
                 sample_group=None, sample_group_color=None, sample_group_is_comm=True,
                 log_base=0, log_additive=1.0, zscore='none',
                 gene_group=None, gene_group_color=None, gene_group_is_comm=False,
                 no_gene_link:bool=False, link_source="https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                 lower_exp_cutoff:float=None, pass_lower_exp_num=None,
                 row_sum_cutoff:float=None, cv_cutoff=0.,
                 target_cols=None, target_rows=None,
                 gene_names=None, sample_names=None,
                 xgap=0.05, ygap=0.05,
                 width:int=None, height:int=None, paper_bgcolor=None, plot_bgcolor=None, sort_cluster_by='distance',
                 gene_label_size:int=None, sample_label_size:int=None, sample_label_angle=45, k_outlier=3.0,
                 color_scale='RdYlBu', reverse_scale=False, preprocess_data_func=None, transpose_data=False,
                 left_dendrogram_width=0.13, top_dendrogram_height=0.13, group_bar_thickness=0.02,
                 colorbar_x:float=None, legend_x:float=None, tmp_keep=False):
        """
        A cluster/heatmap tool designed for gene expression matrix
        * note: gene name should not be pure integer

        For RnaSeq TPM data to cluster samples, best practices:
        refer  http://dx.doi.org/10.1016/j.ymeth.2017.07.023
        (1) average + kendall/spearman (not for deseq2 rlog data as data of gene is not comparable within sample.)
        (2) complete + pearson (not for deseq2 rlog data as data of gene is not comparable within sample.)
        (3) For deseq2 rlog data: I recommend: average + euclidean

        Cluster algorithm choices:
        (1) linkage mathod: ['average', 'weighted', 'centroid', 'complete', 'median', 'ward', 'single']
        (2) distance metric choices:
            ['braycurtis', 'canberra', 'chebyshev', 'cityblock',
            'correlation', 'pearson', 'spearman', 'kendall', (this are correlation metrics)
            'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski',
            'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao',
            'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule', ]
        For more cluster method and metric option detail, please refer scipy.cluster.hierarchy.linkage
        :param data_file: data file path
        :param out_prefix: figure file name prefix, path info can be included
        :param sample_cluster_method: default "average"
        :param sample_distance_metric: default "euclidean", correlation actually refer to pearson corr
        :param gene_cluster_method: default "average"
        :param gene_distance_metric: default "euclidean"
        :param cluster_gene: bool value indicates if to cluster gene
        :param cluster_sample: bool value indicates if to cluster sample
        :param show_gene_label: bool value indicates if to display gene name
        :param hide_sample_label: bool value indicate if to display sample name
        :param hide_legend: bool value indicate if to display group legends of gene and sample
        :param only_sample_dendrogram: bool value indicates if to only draw sample cluster dendrogram
        :param only_gene_dendrogram: bool value indicates if to only draw gene cluster dendrogram
        :param sample_corr_as_heatmap: bool value indicates if to cluster sample using "corr_method" and
            display correlation heat map
        :param gene_corr_as_heatmap: bool value indicates if to cluster gene using "corr_method" and
            display correlation heat map
        :param corr_method: correlation method, could be {'pearson', 'kendall'(default), 'spearman'},
            they are from pandas.corr
        :param sample_cluster_num: number of sample cluster to output
        :param gene_cluster_num: number of gene cluster to output
        :param sample_group: a file with at least two columns, first column is sample name,
            fist row is group scheme name, data among this matrix are group names. Or a pandas data frame object.
        :param sample_group_color: a file with two columns [group, color], or group color dict, {'group': 'color', ...}
        :param sample_group_is_comm: if the group file will be used in other heatmap, this value should be True,
            and this will ensure different heatmap's legend color be common. Default: True
        :param log_base: transform data using log, value could be one of {2, 10, 1, 0,None},
            1 or 0 means do no log transformation
        :param log_additive: a small value added before doing log transformation for data
        :param zscore: str indicates if to do zscore normalization by row or col or none, default: none.
            No effect if "sample/gene_corr_as_heatmap" is True
        :param gene_group: a file with at least two column, first column is gene name and fist row is group scheme,
            data in matrix are group names. or a pandas data frame object.
        :param gene_group_color: a file with two columns [group, color], or group color dict, {'group': 'color', ...}
        :param gene_group_is_comm: if the group file will be used in other heatmap, this value should be True,
            and this will ensure different heatmap's legend color be common. Default: False
        :param no_gene_link: bool to indicate if to cancel make gene link
        :param link_source: default to link to "https://www.genecards.org/cgi-bin/carddisp.pl?gene="
        :param lower_exp_cutoff: gene expression lower cutoff, combined with pass_lower_exp_num
        :param pass_lower_exp_num: gene with expression N times smaller than "lower_exp_cutoff" will be filtered，
            default: int(1/3*sample—number)
        :param row_sum_cutoff: gene with sum of expression lower than this cutoff will be filtered, default None
        :param cv_cutoff: genes with cv (=std/mean) higher than will be retained
        :param target_cols: target columns to extract from data file for analysis
        :param target_rows: target rows to extract from data file for analysis
        :param gene_names: gene annotation file, two columns, gene_id \t gene_symbol
        :param sample_names: sample annotation file, two columns, sample_id \t sample_name
        :param xgap: heatmap xgap
        :param ygap: heatmap ygap
        :param width: figure width, default to auto
        :param height: figure height, default to auto
        :param paper_bgcolor: figure background color, such as 'black'
        :param sort_cluster_by: sort cluster by distance or count
        :param gene_label_size: int, gen label size, default 7 or 10(for sample correlation)
        :param outlier_k: k value for determine outlier, max color value = q3+(q3-q1)*k
        :param color_scale: pallete for heat map, refer to plotly,
            ['Blackbody', 'Bluered', 'Blues', 'Earth', 'Electric',
            'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland',
            'Rainbow', 'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd']
            or one of colorlover scales ['RdYlBu', 'Spectral', 'RdYlGn', 'PiYG', 'PuOr', 'PRGn', 'RdGy']
        :param preprocess_data_func: function provided for data filtering and transformation. Default None
        :param transpose_data: transpose raw data before everything begins
        :param left_dendrogram_width: left/sample dendrogram width, default 0.15, range(0, 1)
        :param top_dendrogram_height: top/gene dendrogram height, default 0.15, range(0, 1)
        :param colorbar_x: colorbar x coordinate, default to self-determined
        :param legend_x: bar legend x coordinate, default to self-determined
        """

        self.scm = sample_cluster_method
        self.sdm = sample_distance_metric
        self.gcm = gene_cluster_method
        self.gdm = gene_distance_metric
        self.scn = sample_cluster_num
        self.gcn = gene_cluster_num
        self.corr_method = corr_method
        self.total_group_num = 0
        self.group_sample = sample_group
        self.sample_group_color = sample_group_color
        self.sample_group_is_comm = sample_group_is_comm
        self.group_gene = gene_group
        self.gene_group_color = gene_group_color
        self.gene_group_is_comm = gene_group_is_comm
        self.transpose_data = transpose_data
        self.sort_cluster_by = sort_cluster_by
        self.outlier_k = k_outlier
        self.link_gene = not no_gene_link
        self.link_source = link_source
        self.gene_names = dict(x.strip().split('\t')[:2] for x in open(gene_names)) if gene_names else None
        self.sample_names = dict(x.strip().split('\t')[:2] for x in open(sample_names)) if sample_names else None
        self.target_cols = [x.strip().split()[0] for x in open(target_cols)] if target_cols else None
        self.target_rows = [x.strip().split()[0] for x in open(target_rows)] if target_rows else None
        self.show_legend = not hide_legend
        self.colorbar_x = colorbar_x
        self.legend_x = legend_x
        self.keep_tmp = tmp_keep
        self.xgap = xgap
        self.ygap = ygap
        if isinstance(sample_group, str):
            if not os.path.exists(sample_group):
                raise Exception('sample group file is not existed')
            self.group_sample = pd.read_csv(sample_group, sep=None, header=0,
                                            index_col=0, engine='python', dtype=str).fillna('Unknown')
            if self.sample_names:
                new_names = [self.sample_names[x] for x in self.group_sample.index if x in self.sample_names]
                if len(new_names) == self.group_sample.shape[0]:
                    self.group_sample.index = new_names
        if isinstance(gene_group, str):
            if not os.path.exists(gene_group):
                raise Exception('sample group file is not existed')
            self.group_gene = pd.read_csv(gene_group, sep=None, header=0, index_col=0, engine='python', dtype=str).fillna('Unknown')
        if isinstance(sample_group_color, str):
            if not os.path.exists(sample_group_color):
                raise Exception('sample color file is not existed')
            with open(sample_group_color) as f:
                self.sample_group_color = dict(line.strip().split('\t')[:2] for line in f)
        if isinstance(gene_group_color, str):
            if not os.path.exists(gene_group_color):
                raise Exception('sample color file is not existed')
            with open(gene_group_color) as f:
                self.gene_group_color = dict(line.strip().split('\t')[:2] for line in f)
        self.gene_label_size = 7 if gene_label_size is None else gene_label_size
        self.sample_label_size = 9 if sample_label_size is None else sample_label_size
        self.sample_label_angle = sample_label_angle
        if sample_corr_as_heatmap and gene_corr_as_heatmap:
            raise Exception('sample_corr_as_heatmap and gene_corr_as_heatmap cannot be both True')
        self.sample_corr_as_heatmap = sample_corr_as_heatmap
        self.gene_corr_as_heatmap = gene_corr_as_heatmap
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
        self.paper_bgcolor = paper_bgcolor
        self.plot_bgcolor = plot_bgcolor
        cs_pool = colorlover.scales['11']['div']
        if color_scale in cs_pool:
            ratio = [x / 100 for x in range(1, 100, int(100 / 9))]
            ratio[0] = 0
            ratio[-1] = 1
            colors = cs_pool[color_scale][1:-1]
            self.colorscale = list(zip(ratio, colors))
        else:
            self.colorscale = color_scale
        self.reverse_scale = reverse_scale
        self.logbase = log_base
        self.log_additive = log_additive
        self.lower_exp_cutoff = lower_exp_cutoff
        self.pass_lower_exp_num = int(pass_lower_exp_num) if pass_lower_exp_num is not None else None
        self.row_sum_cutoff = row_sum_cutoff
        self.cv_cutoff = cv_cutoff
        self.zscore = zscore

        self.out_prefix = out_prefix
        outdir = os.path.dirname(self.out_prefix)
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
            gene_num = 80
            gene_names = ['xyzbacef' + str(x) for x in range(gene_num)]
            self.group_sample = pd.DataFrame(dict(
                status=['gg1', 'gg1', 'gg1', 'gg2', 'gg2', 'gg2'],
                gender=['gg1', 'gg3', 'gg1', 'gg3', 'gg2', 'gg2'],
                age=['gg4', 'gg4', 'gg1', 'gg1', 'gg1', 'gg1'],
            ),
                index=list('abcdef')
            )
            self.group_gene = pd.DataFrame(dict(
                kegg=['gene_group1' if x < 50 else 'gene_group2' for x in range(gene_num)],
                GO=['gene_group3' if x < 70 else 'gene_group1' for x in range(gene_num)],
            ),
                index=gene_names
            )
            self.cluster_gene = True
            self.cluster_sample = True
            self.scd = 'spearman'
            self.label_gene = True
            self.data = pd.DataFrame(np.random.randint(0, 20, (gene_num, 6)),
                                     columns=list('abcdef'),
                                     index=gene_names)

        if self.sample_corr_as_heatmap:
            print('calculate sample correlation')
            self.data = self.data.corr(method=corr_method)
            self.link_gene = not no_gene_link
            self.gene_label_size = 9 if gene_label_size is None else gene_label_size
            # we choose not use the following codes to enable more freedom of users
            # self.cluster_gene = True
            # self.cluster_sample = True
            # self.label_gene = True
        if self.gene_corr_as_heatmap:
            print('calculate gene correlation')
            self.data = self.data.transpose().corr(method=corr_method)
            self.link_gene = not no_gene_link
            self.sample_label_size = 7 if sample_label_size is None else sample_label_size

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

        if self.group_sample is not None:
            self.sample_bar_height = group_bar_thickness* self.group_sample.shape[1]
            if self.only_sample_dendrogram:
                self.top_dendrogram_height = 1 - self.sample_bar_height - 0.01
        else:
            self.sample_bar_height = 0

        if self.group_gene is not None:
            self.gene_bar_width = group_bar_thickness * self.group_gene.shape[1]
            if self.only_gene_dendrogram:
                self.left_dendrogram_width = 1 - self.gene_bar_width - 0.01
        else:
            self.gene_bar_width = 0

        self.layout = self.all_layout()
        self.draw()
        if self.sample_corr_as_heatmap:
            out_corr_file = os.path.join(outdir, '{}.sample.{}.corr.txt'.format(self.out_prefix, self.corr_method))
            self.data.round(4).to_csv(out_corr_file, header=True, index=True, sep='\t')
        if self.gene_corr_as_heatmap:
            out_corr_file = os.path.join(outdir, '{}.gene.{}.corr.txt'.format(self.out_prefix, self.corr_method))
            self.data.round(4).to_csv(out_corr_file, header=True, index=True, sep='\t')

    def process_data(self):
        exp_pd = pd.read_csv(self.data_file, header=0, index_col=0, sep=None, engine='python')
        exp_pd = exp_pd.fillna(0)
        if self.target_rows:
            exp_pd = exp_pd.loc[[x for x in self.target_rows if x in exp_pd.index], :]
        if self.target_cols:
            exp_pd = exp_pd.loc[:, [x for x in self.target_cols if x in exp_pd.columns]]
        if self.transpose_data:
            exp_pd = exp_pd.transpose()
        if self.target_cols or self.target_rows:
            if self.keep_tmp:
                out_name = os.path.join(self.outdir, '{}.target_raw_data'.format(self.out_prefix))
                exp_pd.round(4).to_csv(out_name, header=True, index=True, sep='\t')
        # exp_pd = exp_pd.applymap(lambda x: x if x <=8 else 8)
        if exp_pd.shape[0] <= 1 or exp_pd.shape[1] <= 1:
            raise Exception("Data is not enough for analysis !")
        if self.row_sum_cutoff is not None:
            exp_pd = exp_pd[exp_pd.sum(axis=1) > self.row_sum_cutoff]
        if self.cv_cutoff > 0:
            exp_pd = exp_pd[(exp_pd.std(axis=1)/exp_pd.mean(axis=1)).abs() > self.cv_cutoff]
        pass_num_cutoff = int(exp_pd.shape[1] / 3) if self.pass_lower_exp_num is None else self.pass_lower_exp_num
        if self.lower_exp_cutoff is not None:
            pass_state = exp_pd.apply(lambda x: sum(x > self.lower_exp_cutoff), axis=1)
            exp_pd = exp_pd[pass_state >= pass_num_cutoff]
        if self.logbase == 2:
            exp_pd = np.log2(exp_pd+self.log_additive)
        elif self.logbase == 10:
            exp_pd = np.log10(exp_pd+self.log_additive)
        elif not self.logbase or self.logbase == 1:
            pass
        else:
            raise Exception('log base must be one of [2, 10, 1] ')
        # translate sample_id to sample_name
        if self.sample_names is not None:
            exp_pd.columns = [self.sample_names[x] if x in self.sample_names else x for x in exp_pd.columns]
        if self.keep_tmp:
            out_name = os.path.join(self.outdir, '{}.log{}.cv{}.{}outof{}over{}.data'.format(
                self.out_prefix, self.logbase, self.cv_cutoff, pass_num_cutoff, exp_pd.shape[1], self.lower_exp_cutoff
            ))
            exp_pd.round(4).to_csv(out_name, header=True, index=True, sep='\t')
        # zscore data
        if not (self.sample_corr_as_heatmap or self.gene_corr_as_heatmap):
            if self.zscore == 'row':
                exp_pd = exp_pd.transpose().apply(zscore).transpose()
            elif self.zscore == 'col':
                exp_pd = exp_pd.apply(zscore)
        return exp_pd

    def heatmap_xaxis(self):
        return {
            'domain': [self.left_dendrogram_width + self.gene_bar_width, 1],
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
            'domain': [0, 1 - self.top_dendrogram_height - self.sample_bar_height],
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
            'domain': [0, 1 - self.top_dendrogram_height - self.sample_bar_height],
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
            'domain': [self.left_dendrogram_width+self.gene_bar_width, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': "",
            'anchor': 'y3',
            'scaleanchor': 'x',
            'side': 'top' if not self.only_sample_dendrogram else 'bottom',
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
            'domain': [self.left_dendrogram_width+self.gene_bar_width, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': "",
            'anchor': 'y4',
            'scaleanchor': 'x',
            'side': 'top',
            'tickangle': self.sample_label_angle,
            'tickfont': {'size': self.sample_label_size},
        }

    def group_bar_yaxis(self):
        return {
            'domain': [1 - self.top_dendrogram_height - self.sample_bar_height,
                       1 - self.top_dendrogram_height],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': True,
            'ticks': "",
            'anchor': 'x4',
            'side': 'right' if not self.only_sample_dendrogram else 'left',
        }

    def gene_bar_xaxis(self):
        return {
            'domain': [self.left_dendrogram_width, self.left_dendrogram_width + self.gene_bar_width],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': True,
            'ticks': "",
            'anchor': 'y5',
            'side': 'top' if not self.only_gene_dendrogram else 'bottom',
            'tickangle': 90
        }

    def gene_bar_yaxis(self):
        return {
            'domain':  [0, 1 - self.top_dendrogram_height - self.sample_bar_height],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': "",
            'anchor': 'x5',
            'side': 'left',
            'tickfont': dict(size=self.gene_label_size),
        }

    def all_layout(self):
        if self.height is None and self.label_gene and self.only_sample_dendrogram is False:
            self.height = self.data.shape[0]*(self.gene_label_size+3.5)
            if self.height < 200:
                self.height = 450
            elif self.height < 600:
                self.height = 600
        if (self.sample_corr_as_heatmap or self.gene_corr_as_heatmap) and self.height is None:
            self.height = self.width

        layout = go.Layout(
            plot_bgcolor=self.plot_bgcolor,
            paper_bgcolor=self.paper_bgcolor,
            width=self.width,
            height=self.height,
            autosize=True,
            showlegend=self.show_legend,
            hovermode='closest',
        )
        if self.only_sample_dendrogram:
            layout.update(
                xaxis3=self.top_dendrogram_xaxis(),
                yaxis3=self.top_dendrogram_yaxis(),
                xaxis4=self.group_bar_xaxis(),
                yaxis4=self.group_bar_yaxis(),
            )
            return layout

        if self.only_gene_dendrogram:
            layout.update(
                xaxis2=self.left_dendrogam_xaxis(),
                yaxis2=self.left_dendrogam_yaxis(),
                xaxis5=self.gene_bar_xaxis(),
                yaxis5=self.gene_bar_yaxis()
            )
            return layout

        layout.update(
            xaxis=self.heatmap_xaxis(),
            yaxis=self.heatmap_yaxis(),
            xaxis2=self.left_dendrogam_xaxis(),
            yaxis2=self.left_dendrogam_yaxis(),
            xaxis3=self.top_dendrogram_xaxis(),
            yaxis3=self.top_dendrogram_yaxis(),
            xaxis4=self.group_bar_xaxis(),
            yaxis4=self.group_bar_yaxis(),
            xaxis5=self.gene_bar_xaxis(),
            yaxis5=self.gene_bar_yaxis()
        )
        return layout

    def heatmap_trace(self):
        if not self.ordered_samples:
            self.ordered_samples = range(self.data.shape[1])
        if not self.ordered_genes:
            self.ordered_genes = range(self.data.shape[0])
        # heat_data = self.data.iloc[self.ordered_genes, self.ordered_samples]
        out_name = os.path.join(self.outdir, '{}.heatmap.data'.format(self.out_prefix))
        # plotly plot data from bottom to top, thus we have to use [::-1], or it will not match dendrogram
        heat_data = self.data.iloc[self.ordered_genes[::-1], self.ordered_samples]
        # trans gene id to gene name
        if self.gene_names:
            heat_data.index = [self.gene_names[x] if x in self.gene_names else x for x in heat_data.index]
        # output data
        if self.keep_tmp:
            heat_data.round(4).to_csv(out_name, header=True, index=True, sep='\t')
        if self.link_gene:
            heat_data.index = ["""<a href="{}{}"> {}</a>""".format(self.link_source, x, x) for x in heat_data.index]

        if self.gene_corr_as_heatmap:
            heat_data.columns = heat_data.index

        self.heat_data = heat_data
        # process heat data to make color be more even
        describe = pd.Series(heat_data.values.flatten()).describe()
        upper_limit = describe["75%"] + (describe["75%"] - describe["25%"])*self.outlier_k
        upper_limit = upper_limit if upper_limit < describe['max'] else describe['max']
        lower_limit = describe["25%"] - (describe["75%"] - describe["25%"])*self.outlier_k
        lower_limit = lower_limit if lower_limit > describe['min'] else describe['min']
        max_label_len = max(len(x) for x in self.data.index) + 2
        if self.width is None:
            colorbar_x = 1 + self.gene_label_size*max_label_len/1000 if self.label_gene else 1
        else:
            colorbar_x = 1 + self.gene_label_size*max_label_len/self.width if self.label_gene else 1
        if self.colorbar_x is not None:
            colorbar_x = self.colorbar_x
        heat_map = go.Heatmap(
            x=list(heat_data.columns),
            y=list(heat_data.index),
            z=heat_data.values,
            colorscale=self.colorscale,
            reversescale=not self.reverse_scale,
            showlegend=False,
            xaxis='x',
            yaxis='y',
            xgap=self.xgap,
            ygap=self.ygap,
            zmin=lower_limit,
            zmax=upper_limit,
            name='',
            showscale=True,
            colorbar=dict(
                x=colorbar_x,
                xanchor='left',
                y=0,
                yanchor='bottom',
                len=0.4,
                thickness=25,
                title="log{}(X)".format(self.logbase) if self.logbase and self.logbase != 1 else ''
            )
        )
        return [heat_map]

    def group_bar_traces(self):
        if not self.ordered_samples:
            self.ordered_samples = range(self.data.shape[1])
        ordered_samples = self.data.columns[self.ordered_samples]
        not_grouped = [x for x in ordered_samples if x not in self.group_sample.index]
        for each in not_grouped:
            self.group_sample.loc[each, :] = 'Unknown'
        group_df = self.group_sample.loc[ordered_samples, :]
        all_group_dict = group_df.transpose().to_dict('index')

        if not self.sample_group_is_comm:
            groups = sorted(set(group_df.values.flatten()))
        else:
            groups = sorted(set(self.group_sample.values.flatten()))

        gene_group_num = 0
        if self.group_gene is not None and self.ordered_genes:
            if sum(1 for x in self.ordered_genes if x not in self.group_gene.index) > 0:
                gene_group_num += 1
            if self.gene_group_is_comm:
                gene_group_num += len(set(self.group_gene.values.flatten()))
            else:
                target_genes = set(self.group_gene.index) & set(self.ordered_genes)
                gene_group_num += len(set(self.group_gene.loc[list(target_genes), :].values.flatten()))
        self.total_group_num = len(groups) + gene_group_num
        if self.gene_corr_as_heatmap:
            group_colors = self.gene_group_color
        else:
            colors = self.get_color_pool(self.total_group_num)
            group_colors = dict(zip(groups, colors))
            # if 'Unknown' in group_colors:
            #     group_colors['Unknown'] = 'darkgrey'
            if self.sample_group_color:
                # user defined color overrides random generated one
                for k, v in self.sample_group_color.items():
                    group_colors[k] = v
            self.sample_group_color = group_colors
        # plot
        traces = list()
        existed_legend = set()
        base = -1.01
        # if self.keep_tmp:
        with open('sample.group.colors', 'w') as f:
            for k, v in group_colors.items():
                f.write('{}\t{}\n'.format(k, v))
        for category, group_dict in all_group_dict.items():
            base += 1
            sample_colors = dict()
            for sample in ordered_samples:
                sample_colors[sample] = group_colors[group_dict[sample]]
            ticks = range(5, self.data.shape[1]*10, 10)
            for tick, each in zip(ticks, ordered_samples):
                if each == group_dict[each]:
                    trace_name = each
                else:
                    trace_name = '{sample}|{group}'.format(sample=each, group=group_dict[each])
                show_legend = True if sample_colors[each] not in existed_legend else False
                if self.gene_corr_as_heatmap:
                    show_legend = False
                bar = go.Bar(
                    name=group_dict[each],
                    legendgroup=group_dict[each],
                    text=trace_name,
                    x=[tick],
                    y=[1],
                    base=base,
                    width=10,
                    showlegend=show_legend,
                    xaxis='x4',
                    yaxis='y4',
                    marker=dict(color=sample_colors[each]),
                    hoverinfo="text",
                    hoverlabel=dict(namelength=-1)
                )
                traces.append(bar)
                existed_legend.add(sample_colors[each])
            # break
        self.layout['yaxis4']['tickvals'] = [x+0.5 for x in range(len(all_group_dict))]
        self.layout['yaxis4']['ticktext'] = list(all_group_dict.keys())
        self.layout['yaxis4']['tickfont'] = dict(size=self.sample_label_size)
        if self.only_sample_dendrogram and self.group_sample is not None:
            self.layout['margin'] = dict(l=max(len(x) for x in self.group_sample.columns)*self.sample_label_size)
            self.layout['xaxis4']['showticklabels'] = True
            self.layout['xaxis4']['side'] = 'bottom'
            self.layout['xaxis4']['tickvals'] = list(range(5, self.data.shape[1]*10, 10))
            self.layout['xaxis4']['ticktext'] = ordered_samples
        return traces

    def gene_bar_traces(self):
        if not self.ordered_genes:
            self.ordered_genes = range(self.data.shape[0])
        ordered_genes = self.data.index[self.ordered_genes[::-1]]
        not_grouped = [x for x in ordered_genes if x not in self.group_gene.index]
        for each in not_grouped:
            self.group_gene.loc[each, :] = 'Unknown'
        group_df = self.group_gene.loc[ordered_genes, :]
        all_group_dict = group_df.transpose().to_dict('index')

        if not self.gene_group_is_comm:
            groups = sorted(set(group_df.values.flatten()))
        else:
            groups = sorted(set(self.group_gene.values.flatten()))

        sample_group_num = 0
        if self.group_gene is not None and self.ordered_samples:
            if sum(1 for x in self.ordered_samples if x not in self.group_sample.index) > 0:
                sample_group_num += 1
            if self.sample_group_is_comm:
                sample_group_num += len(set(self.group_sample.values.flatten()))
            else:
                target_samples = set(self.group_sample.index) & set(self.ordered_samples)
                sample_group_num += len(set(self.group_gene.loc[list(target_samples), :].values.flatten()))
        self.total_group_num = len(groups)+sample_group_num
        if self.sample_corr_as_heatmap:
            group_colors = self.sample_group_color
        else:
            colors = self.get_color_pool(self.total_group_num)
            colors = colors[-len(groups):]
            group_colors = dict(zip(groups, colors))
            # if 'Unknown' in group_colors:
            #     group_colors['Unknown'] = 'darkgrey'
            if self.gene_group_color:
                # user defined color overrides random generated one
                for k, v in self.gene_group_color.items():
                    group_colors[k] = v
            self.gene_group_color = group_colors
        # insert link
        if self.only_gene_dendrogram:
            self.layout['yaxis5']['showticklabels'] = True
            self.layout['yaxis5']['side'] = 'right'
            self.layout['yaxis5']['dtick'] = 1
            labels = [x for x in ordered_genes]
            if self.gene_names:
                new_labels = []
                for x in labels:
                    if x in self.gene_names:
                        new_labels.append(
                            """<a href="{}{}">{}</a>""".format(self.link_source, self.gene_names[x], self.gene_names[x])
                        )
                    else:
                        new_labels.append(
                            """<a href="{}{}">{}</a>""".format(self.link_source, x, x)
                        )
                labels = new_labels
            elif self.link_gene:
                labels = ["""<a href="{}{}">{}</a>""".format(self.link_source, x, x) for x in labels]
        # plot
        traces = list()
        existed_legend = set()
        base = -1.01
        with open('gene.group.colors', 'w') as f:
            for k, v in group_colors.items():
                f.write('{}\t{}\n'.format(k, v))
        for category, group_dict in all_group_dict.items():
            base += 1
            sample_colors = dict()
            for sample in ordered_genes:
                sample_colors[sample] = group_colors[group_dict[sample]]
            for ind, each in enumerate(ordered_genes):
                if each == group_dict[each]:
                    trace_name = each
                else:
                    trace_name = '{sample}|{group}'.format(sample=each, group=group_dict[each])
                show_legend = True if sample_colors[each] not in existed_legend else False
                if self.sample_corr_as_heatmap:
                    show_legend = False
                bar = go.Bar(
                    orientation='h',
                    name=group_dict[each],
                    legendgroup=group_dict[each],
                    text=trace_name,
                    x=[1],
                    y=[each if not self.only_gene_dendrogram else labels[ind]],
                    base=base,
                    width=1,
                    showlegend=show_legend,
                    xaxis='x5',
                    yaxis='y5',
                    marker=dict(color=sample_colors[each]),
                    hoverinfo="text",
                    hoverlabel=dict(namelength=-1)
                )
                traces.append(bar)
                existed_legend.add(sample_colors[each])
            # break
        self.layout['xaxis5']['tickvals'] = [x + 0.5 for x in range(len(all_group_dict))]
        self.layout['xaxis5']['ticktext'] = list(all_group_dict.keys())
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
            prefix='{}.sample.'.format(self.out_prefix)
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
        icoord = np.array(results['icoord'])
        dcoord = np.array(results['dcoord'])
        color_list = np.array(results['color_list'])
        trace_list = []
        if self.top_dendrogram_height > 0:
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
            if self.group_sample is None:
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
            prefix='{}.gene.'.format(self.out_prefix),
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
        icoord = np.array(results['dcoord'])*(-1)
        dcoord = np.array(results['icoord'])*(-1)
        color_list = np.array(results['color_list'])
        trace_list = []
        if self.left_dendrogram_width > 0:
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
            if self.group_gene is None:
                self.layout['yaxis2']['showticklabels'] = True
                tick_values = list(range(5, (exp_pd.shape[0]+1)*10, 10))
                self.layout['yaxis2']['tickvals'] = [-1*x for x in tick_values]
                self.layout['yaxis2']['side'] = 'right'
                labels = exp_pd.iloc[self.ordered_genes].index
                if self.gene_names:
                    new_labels = []
                    for x in labels:
                        if x in self.gene_names:
                            new_labels.append(
                                """<a href="{}{}"> {}</a>""".format(self.link_source, self.gene_names[x], self.gene_names[x])
                            )
                        else:
                            new_labels.append(
                                """<a href="{}{}"> {}</a>""".format(self.link_source, x, x)
                            )
                    labels = new_labels
                elif self.link_gene:
                    labels = ["""<a href="{}{}"> {}</a>""".format(self.link_source, x, x) for x in labels]
                self.layout['yaxis2']['ticktext'] = labels

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
            Choices: ['braycurtis', 'canberra', 'chebyshev', 'cityblock',
                    'correlation', 'pearson', 'spearman', 'kendall',
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

        if self.sample_corr_as_heatmap or self.gene_corr_as_heatmap:
            self.logbase = None
            condensed_distance = squareform(1 - exp_pd.abs())
            z = hclust.linkage(condensed_distance, method=self.scm)
        else:
            if metric.lower() == 'pearson':
                metric = 'correlation'
            if metric in ['spearman', 'kendall']:
                exp_pd = exp_pd.transpose().corr(method=metric)
                exp_pd = squareform(1 - exp_pd)
                metric = 'correlation'  # no effect, just for save
            try:
                z = hclust.linkage(exp_pd, method=method, metric=metric)
            except FloatingPointError as e:
                print('it seems that (at least) one of the vectors you want to cluster is all zeros, '
                      'so when it tries to compute the cosine distances to it there is a division by zero,'
                      ' hence nan is stored in your distance array, and that leads to your error.')
                raise Exception("fastcluster failed as : {}".format(e))
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
            out_dir = os.path.join(output, prefix + 'subcluster_{}_{}.xls'.format(k,len(v)))
            sub = exp_pd.loc[v, :]
            if transpose:
                sub = sub.transpose()
            sub.round(4).to_csv(out_dir, sep='\t', header=True, index=True)
        # write out cluster result z
        out_dir = os.path.join(output, prefix + "linkage_result")
        pd.DataFrame(z).round(4).to_csv(out_dir, sep='\t')
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
        colors = colorlover.scales['7']['qual']['Set1']
        sch.set_link_color_palette(colors)

    @staticmethod
    def get_color_pool(n):
        # https://plot.ly/ipython-notebooks/color-scales/
        import colorlover
        if n <= 8:
            if n <= 3:
                n = 3
            return colorlover.scales[str(n)]['qual']['Set2']
        if n <= 12:
            return colorlover.scales[str(n)]['qual']['Set3']

        import random
        random.seed(666)

        def get_random_color(pastel_factor=0.5):
            return [(x + pastel_factor) / (1.0 + pastel_factor) for x in [random.uniform(0, 1.0) for i in [1, 2, 3]]]

        def color_distance(c1, c2):
            return sum([abs(x[0] - x[1]) for x in zip(c1, c2)])

        def generate_new_color(existing_colors, pastel_factor=0.5):
            max_distance = None
            best_color = None
            for i in range(0, 100):
                color = get_random_color(pastel_factor=pastel_factor)
                # exclude some colors
                if np.absolute(np.array(color) - np.array([1, 1, 1])).sum() < 0.08:
                    continue
                if not existing_colors:
                    return color
                best_distance = min([color_distance(color, c) for c in existing_colors])
                if not max_distance or best_distance > max_distance:
                    max_distance = best_distance
                    best_color = color
            return best_color

        color_pool = []
        for i in range(0, n):
            color_pool.append(generate_new_color(color_pool, pastel_factor=0.9))
        color_pool = [(int(x*255), int(y*255), int(z*255)) for x,y,z in color_pool]
        return colorlover.to_rgb(color_pool)

    def draw(self):
        if self.group_sample is not None or self.group_gene is not None:
            self.layout['barmode'] = 'stack'
        if self.show_legend:
            self.layout['legend'] = dict(
                tracegroupgap=2,
                traceorder='grouped',
            )
        if self.only_gene_dendrogram:
            self.layout['legend'].update(
                x=self.left_dendrogram_width,
                y=1.05,
                bgcolor='rgba(0,0,0,0)',
                orientation="h"
            )
        traces = list()
        if self.only_sample_dendrogram:
            traces += self.top_dendrogram_traces()
            if self.group_sample is not None:
                traces += self.group_bar_traces()
            fig = go.Figure(data=traces, layout=self.layout)
            plt(fig, filename=self.out_prefix+'.html', auto_open=False)
            return

        if self.only_gene_dendrogram:
            traces += self.left_dendrogram_traces()
            if self.group_gene is not None:
                traces += self.gene_bar_traces()
            fig = go.Figure(data=traces, layout=self.layout)
            plt(fig, filename=self.out_prefix+'.html', auto_open=False)
            return

        if self.cluster_sample:
            traces += self.top_dendrogram_traces()

        if self.cluster_gene:
            traces += self.left_dendrogram_traces()

        if self.group_sample is not None:
            traces += self.group_bar_traces()

        if self.group_gene is not None:
            traces += self.gene_bar_traces()

        traces += self.heatmap_trace()
        if self.show_legend:
            height = self.height or 600
            if self.total_group_num > 0.6*self.heat_data.shape[0] or \
                    (height < 601 and self.total_group_num > 9):
                legend_x = traces[-1]['colorbar'] ['x'] + 0.07
            else:
                legend_x = traces[-1]['colorbar']['x']
            self.layout['legend'].update(
                x=legend_x if self.legend_x is None else self.legend_x,
            )

        fig = go.Figure(data=traces, layout=self.layout)
        plt(fig, filename=self.out_prefix+'.html', auto_open=False)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['ClusterHeatMap'])

