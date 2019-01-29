import plotly.plotly as py
import plotly.graph_objs as go
import pandas as pd
import scipy as scp
import numpy as np
from scipy.cluster import hierarchy as sch


class ClusterHeatMap():
    def __init__(self, data_file):
        self.data_file = data_file
        self.Z = self.linkage()
        self.layout = self.layout()

    def heatmap_xaxis(self):
        return {
            'domain': [.15, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'ticks':""
        }

    def heatmap_yaxis(self):
        return {
            'domain': [0, .85],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        }

    def right_dendrogam_xaxis(self):
        return {
            'domain': [0, .15],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks':""
        }

    def right_dendrogam_yaxis(self):
        return {
            'domain': [0, .85],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        }

    def top_dendrogram_xaxis(self):
        return {
            'domain': [0.15, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        }

    def top_dendrogram_yaxis(self):
        return {
            'domain': [.85, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        }

    def layout(self):
        return go.Layout(
            width=800,
            height=800,
            showlegend=False,
            hovermode='closest',
            xaxis=self.heatmap_xaxis(),
            yaxis=self.heatmap_yaxis(),
            xaxis2=self.right_dendrogam_xaxis(),
            yaxis2=self.right_dendrogam_yaxis(),
            xaxis3=self.top_dendrogram_xaxis(),
            yaxis3=self.top_dendrogram_yaxis(),
        )

    def heatmap_traces(self):
        pass

    def top_dendrogram_traces(self):
        P = sch.dendrogram(
            self.Z,
            orientation=self.orientation,
            labels=self.labels,
            no_plot=True,
            color_threshold=color_threshold
        )
        icoord = scp.array(P['icoord'])
        dcoord = scp.array(P['dcoord'])
        ordered_labels = scp.array(P['ivl'])
        color_list = scp.array(P['color_list'])
        colors = self.get_color_dict(colorscale)
        trace_list = []
        for i in range(len(icoord)):
            # xs and ys are arrays of 4 points that make up the '∩' shapes
            # of the dendrogram tree
            if self.orientation in ['top', 'bottom']:
                xs = icoord[i]
            else:
                xs = dcoord[i]
            if self.orientation in ['top', 'bottom']:
                ys = dcoord[i]
            else:
                ys = icoord[i]
            hovertext_label = None
            trace = go.Scatter(
                type='scatter',
                x=np.multiply(self.sign[self.xaxis], xs),
                y=np.multiply(self.sign[self.yaxis], ys),
                mode='lines',
                marker=dict(color=color_list[i]),
                text=hovertext_label,
                hoverinfo='text',
                xaxis="x3",
                yaxis="y3"
            )
            trace_list.append(trace)

    def right_dendrogram_traces(self):
        pass

    def linkage(self):
        from sklearn.metrics.pairwise import pairwise_distances
        from scipy.cluster.hierarchy import linkage, dendrogram

        X = pairwise_distances(df.values, metric='euclidean')
        Z = linkage(X, 'ward')
        results = dendrogram(Z, no_plot=True)
        icoord, dcoord = results['icoord'], results['dcoord']
        labels = list(map(int, results['ivl']))
        df = df.iloc[labels]
        df_scaled = df_scaled.iloc[labels]
        pass

    def data(self):
        pass




