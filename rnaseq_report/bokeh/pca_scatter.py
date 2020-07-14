import os
import pandas as pd
import numpy as np
from sklearn import decomposition, preprocessing
from bokeh.plotting import figure, save, output_file
from bokeh.models import ColumnDataSource, CDSView, GroupFilter, HoverTool, LabelSet, Legend
from bokeh.layouts import gridplot


def pca(table, row_sum_cutoff=1, exp_cutoff=0.5, cv_cutoff=0.01, pass_exp_cutoff_num=None,
        explained_ratio=0.9, prefix='pca', no_log_transform=False, log_additive=0, no_pre_scale=False,
        x=1, y=2, group_file=None, annotate=False, text_size='6pt', marker_size=15, stretch='both'):
    # data = pd.read_csv(table, header=0, index_col=0, sep=None, engine='python')
    data = pd.read_csv(table, header=0, index_col=0, sep=None, engine='python').fillna(0)
    data = data[data.sum(axis=1) >= row_sum_cutoff]
    pass_state = data.apply(lambda x: sum(x > exp_cutoff), axis=1)
    if pass_exp_cutoff_num is None:
        data = data[pass_state >= int(data.shape[1]) / 3]
    else:
        data = data[pass_state >= pass_exp_cutoff_num]
    data = data[data.std(axis=1) / data.mean(axis=1) > cv_cutoff]
    data.round(4).to_csv('{}.filtered.data.txt'.format(prefix), header=True, index=True, sep='\t')
    if not no_log_transform:
        data = np.log(data + log_additive)
    if not no_pre_scale:
        data = data.apply(preprocessing.scale, axis=0)
    data = data.transpose()
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
    result.index.name = 'sample'
    result.columns = ['PC' + str(n) + '(' + '{:.2%}'.format(r) + ')' for n, r in _ratio]
    # result.columns = ['PC' + str(n) for n in range(1, result.shape[1] + 1)]
    out_file = '{}.xls'.format(prefix)
    result.round(4).to_csv(out_file, sep='\t', header=True, index=True)
    pc_ratio = {'PC' + str(n): r for n, r in _ratio}
    print(pc_ratio)
    # plot
    plt_data = result.iloc[:, [x-1, y-1]]
    if group_file is not None:
        group_info = pd.read_csv(group_file, sep=None, index_col=0, header=0, engine='python')
        plt_data = plt_data.join(group_info)
    pca_scatter(plt_data, out_file='{}.html'.format(prefix), annotate=annotate,
                text_size=text_size, marker_size=marker_size, stretch=stretch)
    return result


def pca_scatter(df=None, out_file='pca.html', annotate=False, text_size='6pt', marker_size=15, stretch='both'):
    # sample
    if df is None:
        df = pd.DataFrame(dict(
                    x=[1, 2, 3, 4, 5],
                    y=[2, 5, 8, 2, 7],
                    group=['g1', 'g1', 'g1', 'g2', 'g2'],
                    group2=['g3', 'g3', 'g4', 'g4', 'g4']
                ))
        df.index = ['A', 'b', 'C', 'd', 'E']
        df.index.name = 'name'

    # get group scheme
    if df.shape[1] < 3:
        df.loc[:, 'group'] = df.index

    group_scheme_list = list(df.columns[2:])

    # all markers
    marker_pool = ["asterisk", "circle", "cross", "square", "diamond", "triangle",
                   "inverted_triangle", "x", "diamond_cross", "circle_x",
                   "square_cross", "square_x", "dash", "circle_cross"]

    # match markers
    for group_scheme in group_scheme_list:
        groups = set(df[group_scheme])
        color_pool = get_color_pool(len(groups))
        if len(groups) <= len(marker_pool):
            marker_dict = dict(zip(groups, marker_pool))
            color_dict = dict(zip(groups, color_pool))
            df.loc[:, group_scheme + '_marker'] = [marker_dict[x] for x in df[group_scheme]]
            df.loc[:, group_scheme + '_color'] = [color_dict[x] for x in df[group_scheme]]
        else:
            color_dict = dict(zip(groups, color_pool))
            df.loc[:, group_scheme + '_marker'] = "circle"
            df.loc[:, group_scheme + '_color'] = [color_dict[x] for x in df[group_scheme]]

    #
    source = ColumnDataSource(df)
    plot_options = dict(
        # width=250,
        # plot_height=500,
        tools='pan,wheel_zoom,box_select, reset,save',
    )

    # plot
    plots = list()
    for ind, group_scheme in enumerate(group_scheme_list):
        if ind > 0:
            plot_options['x_range'] = plots[0].x_range
            plot_options['y_range'] = plots[0].y_range
        s = figure(**plot_options)
        hover = HoverTool(
            tooltips=[
                ("group", "@{}".format(df.index.name)),
            ]
        )
        s.add_tools(hover)
        groups = set(df[group_scheme])
        legend_items = list()
        for group in groups:
            tmp = s.scatter(
                x=df.columns[0],
                y=df.columns[1],
                marker=group_scheme+'_marker',
                source=source,
                size=marker_size,
                color=group_scheme+'_color',
                alpha=0.6,
                # legend=group,
                view=CDSView(source=source, filters=[GroupFilter(column_name=group_scheme, group=group)])
            )
            legend_items.append((group, [tmp]))
        # 如此可以保证legend在图形外面
        legend = Legend(items=legend_items,location="center")
        s.add_layout(legend, 'right')

        if annotate:
            labels = LabelSet(x=df.columns[0],
                              y=df.columns[1],
                              text=df.index.name,
                              level='glyph',
                              x_offset=5,
                              y_offset=2,
                              source=source,
                              text_font_size=text_size,
                              render_mode='canvas')
            s.add_layout(labels)

        s.legend.location = 'top_left'
        s.legend.click_policy = "hide"
        s.xaxis.axis_label = df.columns[0]
        s.yaxis.axis_label = df.columns[1]
        plots.append(s)

    p = gridplot(plots, sizing_mode='stretch_{}'.format(stretch), ncols=2)
    output_file(out_file, title="PCA Scatter")
    save(p)


def get_color_pool(n):
    # https://plot.ly/ipython-notebooks/color-scales/
    import colorlover
    if n <= 8:
        if n <= 3:
            n = 3
        return colorlover.scales[str(n)]['qual']['Set1']
    if n <= 12:
        return colorlover.scales[str(n)]['qual']['Paired']

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
            if np.absolute(np.array(color) - np.array([1, 1, 1])).sum() < 0.1:
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
    color_pool = [(int(x * 255), int(y * 255), int(z * 255)) for x, y, z in color_pool]
    color_pool = sorted(color_pool, key=lambda x: (x[0], x[1], x[2]))
    return colorlover.to_rgb(color_pool)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['pca'])
