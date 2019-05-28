import pandas as pd
pd.set_option('mode.chained_assignment','raise')
import numpy as np
import math
from scipy import stats
from bokeh.io import output_file
from bokeh.layouts import layout, column, row
from bokeh.plotting import figure, save
from bokeh.palettes import RdYlGn
from bokeh.transform import linear_cmap
from bokeh.models import (
    ColumnDataSource, CustomJS,
    CDSView, GroupFilter, HoverTool,
    LabelSet, Legend
)


def volcano_source(table, gene_symbol_ind=2, log2fc_ind=3, pvalue_ind=5,
                   fc_cutoff=2.0, pval_cutoff=0.05, limit=6, exp_ind=('6-',)):
    df_raw = pd.read_csv(table, index_col=0, header=0, sep=None, engine='python')
    df_raw.index.name = 'index'
    if gene_symbol_ind >= 2:
        df_raw.index = df_raw.index + '|' + df_raw.iloc[:, gene_symbol_ind-2]

    # prepare scatter source
    df = df_raw.iloc[:, [log2fc_ind-2, pvalue_ind-2]].copy()
    df.columns = ['log2fc', 'pvalue']
    drop_before = df.shape[0]
    df.loc[:, ['log2fc', 'pvalue']] = df[['log2fc', 'pvalue']].apply(pd.to_numeric, errors='coerce')
    df = df.dropna(axis=0)
    df = df.loc[(df['log2fc'].abs() > 0)]
    drop_after = df.shape[0]
    if drop_after-drop_before:
        print('drop', drop_after-drop_before, 'lines')
    df['regulate'] = 'NotSig'
    up = (df['log2fc'] >= math.log2(fc_cutoff)) & (df['pvalue'] <= pval_cutoff)
    down = (df['log2fc'] <= -math.log2(fc_cutoff)) & (df['pvalue'] <= pval_cutoff)
    df.loc[up, 'regulate'] = 'Up'
    df.loc[down, 'regulate'] = 'Down'
    fc_desc = df['log2fc'].describe()
    fc_upper_limit = fc_desc['75%'] + (fc_desc['75%'] - fc_desc['25%']) * limit
    fc_lower_limit = fc_desc['25%'] - (fc_desc['75%'] - fc_desc['25%']) * limit
    df['log2fc_shrink'] = df.loc[:, 'log2fc']
    df.loc[(df['log2fc_shrink'] >= fc_upper_limit), 'log2fc_shrink'] = fc_upper_limit
    df.loc[(df['log2fc_shrink'] <= fc_lower_limit), 'log2fc_shrink'] = fc_lower_limit
    df['pvalue_shrink'] = df.loc[:, 'pvalue']
    df.loc[df['pvalue'] == 0, 'pvalue_shrink'] = df['pvalue'][df['pvalue'] > 0].min()
    df['pvalue_shrink'] = -np.log10(df['pvalue_shrink'])
    p_desc = df['pvalue_shrink'].describe()
    p_limit = p_desc['75%'] + (p_desc['75%'] - p_desc['25%']) * limit
    df.loc[(df['pvalue_shrink'] >= p_limit), 'pvalue_shrink'] = p_limit
    df['color'] = 'darkgrey'
    df.loc[df['regulate']=='Up', 'color'] = 'tomato'
    df.loc[df['regulate']=='Down', 'color'] = 'mediumseagreen'

    # prepare exp bar source
    if len(exp_ind) == 1 and type(exp_ind[0]) is str:
        ind_str= exp_ind[0].split('-')
        if len(ind_str) == 1:
            start = int(ind_str[0])
            end = df_raw.shape[1]
        else:
            start = int(ind_str[0])
            end = int(ind_str[1])
        exp_bar = df_raw.iloc[:, start-2:end-1]
    else:
        exp_ind = [int(x)-2 for x in exp_ind]
        exp_bar = df_raw.iloc[:, exp_ind]

    return df, exp_bar


def plot(volcano_df, bar_df, out_file='volcano_expression.html', corr_method='pearson', corr_test=False,
         top=10, corr_pval_cutoff=0.01, label_corr_line=True, only_link_diff_gene=True):
    output_file(out_file)
    first_gene = volcano_df.index[volcano_df['regulate'] != 'NotSig' ][0]
    # filter out diff gene for analysis
    if only_link_diff_gene:
        diff_index = volcano_df.loc[volcano_df['regulate'] != 'NotSig'].index
        bar_df = bar_df.loc[diff_index]
        print('bar_df shape', bar_df.shape)
    else:
        bar_df = bar_df.loc[volcano_df.index]
    samples = list(bar_df.columns)
    point_groups = set(volcano_df['regulate'])
    point_stat = {x: x+': '+str(list(volcano_df['regulate']).count(x)) for x in point_groups}
    volcano_df['regulate'] = [point_stat[x] for x in volcano_df['regulate']]
    volcano_source = ColumnDataSource(volcano_df.round(2))
    volcano_gene_list = list(volcano_df.index)
    bar_dict_source = bar_df.round(2).transpose().to_dict('list')

    # calculate correlation between genes
    corr, corr_pval = df_corr(bar_df.transpose(), method=corr_method, test=corr_test)
    corr_dict_source = dict()
    for gene in corr.index:
        if corr_test:
            pval_pass = corr_pval.loc[gene] < corr_pval_cutoff
            tmp_corr = corr.loc[gene][pval_pass].abs().sort_values(ascending=False)
        else:
            tmp_corr = corr.loc[gene].abs().sort_values(ascending=False)
        if top+1 >= tmp_corr.shape[0]:
            ind = tmp_corr.index
        else:
            ind = tmp_corr.index[range(top+1)]
        ind = [x for x in ind if x != gene]
        if not ind:
            ind = [gene]
        top_corr = corr.loc[gene].loc[ind].round(3)
        if corr_test:
            top_corr_pval = corr_pval.loc[gene].loc[ind]
            corr_dict_source[gene] = [list(ind), list(top_corr), list(top_corr_pval)]
        else:
            corr_dict_source[gene] = [list(ind), list(top_corr)]

    # circle plot
    plot_options = dict(
        # width=250,
        # plot_height=500,
        tools='pan,wheel_zoom,box_select, reset,save',
        toolbar_location="above"
    )
    volcano = figure(**plot_options)
    volcano.circle(
        x='log2fc_shrink',
        y='pvalue_shrink',
        source=volcano_source,
        color='color',
        legend='regulate'
    )
    volcano.legend.location = 'bottom_center'
    volcano.legend.orientation = 'horizontal'
    volcano.legend.background_fill_alpha = 0.6
    volcano.xaxis.axis_label = 'log2(FoldChange)'
    volcano.yaxis.axis_label = '-log10(Pvalue)'

    # corr bar
    corr_bar = figure(**plot_options, x_range=corr_dict_source[first_gene][0])
    if corr_test:
        dynamic2 = ColumnDataSource({
            'x': corr_dict_source[first_gene][0],
            'y': corr_dict_source[first_gene][1],
            'pval': corr_dict_source[first_gene][2]
        })
    else:
        dynamic2 = ColumnDataSource({
            'x': corr_dict_source[first_gene][0],
            'y': corr_dict_source[first_gene][1],
        })
    mapper = linear_cmap(palette=RdYlGn[11], field_name='y', low=-1, high=1)
    corr_bar_sr = corr_bar.vbar(x='x', top='y', width=0.5, source=dynamic2, color=mapper)
    corr_bar.xaxis.major_label_orientation = math.pi / 4
    corr_bar.yaxis.axis_label = 'Correlation'
    corr_bar.xgrid.grid_line_color = None
    corr_bar.title.text = first_gene

    # corr line
    corr_line = figure(**plot_options, x_axis_location='above', y_axis_location='right')
    dynamic3 = ColumnDataSource({
        'x': [0]*len(samples),
        'y': [0]*len(samples),
        'samples': samples
    })
    # corr_line.line(x='x', y='y', source=dynamic3, color='tomato')
    corr_line.circle(x='x', y='y', size=12, source=dynamic3, alpha=0.5)
    from bokeh.models import LinearAxis
    # 为了能再callback时更新label信息
    corr_line.xaxis.visible = None
    corr_line.yaxis.visible = None
    corr_line_xaxis = LinearAxis(axis_label="gene1")
    corr_line_yaxis = LinearAxis(axis_label="gene2")
    corr_line.add_layout(corr_line_xaxis, 'below')
    corr_line.add_layout(corr_line_yaxis, 'left')
    corr_line.title.text = 'gene1 vs gene2'
    corr_line.add_tools(HoverTool(
        # tooltips=[('sample', '@samples'),('gene1_expr', '@x'), ('gene2_expr', '@y')])
        tooltips=[('sample', '@samples')])
    )
    if label_corr_line:
        labels = LabelSet(x='x',
                          y='y',
                          text='samples',
                          level='glyph',
                          x_offset=0,
                          y_offset=0,
                          source=dynamic3,
                          text_font_size='6pt',
                          render_mode='canvas')
        corr_line.add_layout(labels)

    # bar plot
    exp_bar = figure(**plot_options, x_range=samples)
    dynamic = ColumnDataSource({
        'x': samples,
        'y': bar_df.loc[first_gene]
    })
    describe = pd.Series(bar_df.values.flatten()).describe()
    upper_limit = describe["75%"] + (describe["75%"] - describe["25%"]) * 3
    upper_limit = upper_limit if upper_limit < describe['max'] else describe['max']
    lower_limit = describe["25%"] - (describe["75%"] - describe["25%"]) * 3
    lower_limit = lower_limit if lower_limit > describe['min'] else describe['min']
    print('Exp Color range: ', lower_limit, upper_limit)
    mapper = linear_cmap(palette=RdYlGn[11], field_name='y', low=lower_limit, high=upper_limit)
    exp_bar.vbar(x='x', top='y', width=0.5, source=dynamic, color=mapper)
    exp_bar.xaxis.major_label_orientation = math.pi / 4
    exp_bar.yaxis.axis_label = 'Expression'
    exp_bar.xgrid.grid_line_color = None
    exp_bar.title.text = first_gene

    # interaction
    js_args = {
        'exp': bar_dict_source,
        'dynamic': dynamic,
        'exp_bar_title': exp_bar.title,
        'genes': volcano_gene_list,
        'samples': samples,
        'exp_bar_x': exp_bar.x_range,
        'corr': corr_dict_source,
        'dynamic2': dynamic2,
        'corr_bar_x': corr_bar.x_range,
        'title2': corr_bar.title
    }
    if corr_test:
        js_code = """
                    var ind = cb_data.index['1d'].indices[0];
                    var gene = genes[ind];
                    
                    var d = dynamic2.data;
                    d['x'] = corr[gene][0];
                    d['y'] = corr[gene][1];
                    d['pval'] = corr[gene][2];
                    corr_bar_x.factors = corr[gene][0];
                    title2.text = gene;
                    dynamic2.change.emit();
                    
                    var d3 = dynamic.data;
                    d3['x'] = samples;
                    d3['y'] = exp[gene];
                    exp_bar_x.factors = samples;
                    exp_bar_title.text = gene;
                    dynamic.change.emit();
                """
    else:
        js_code = """
                    var ind = cb_data.index['1d'].indices[0];
                    var gene = genes[ind];

                    var d = dynamic2.data;
                    d['x'] = corr[gene][0];
                    d['y'] = corr[gene][1];
                    corr_bar_x.factors = corr[gene][0];
                    title2.text = gene;
                    dynamic2.change.emit();

                    var d3 = dynamic.data;
                    d3['x'] = samples;
                    d3['y'] = exp[gene];
                    exp_bar_x.factors = samples;
                    exp_bar_title.text = gene;
                    dynamic.change.emit();
                """
    on_hover = CustomJS(args=js_args, code=js_code)

    # define tools
    hover = HoverTool(
        tooltips=[
            ('gene', '@index'),
            ('log2fc', '@log2fc'),
            ('pvalue', '@pvalue')
        ],
        callback=on_hover,
        point_policy="snap_to_data",
    )
    volcano.add_tools(hover)

    # exp_bar.yaxis.axis_label = first_gene
    bar_hover = HoverTool(
        tooltips=[
            ('sample', '@x'),
            ('expr', '@y'),
        ],
    )
    exp_bar.add_tools(bar_hover)

    # corr hover
    # corr_hover = HoverTool(
    #     tooltips=[
    #         ('gene', '@x'),
    #         ('corr', '@y'),
    #     ],
    # )
    # corr_bar.add_tools(corr_hover)

    # corr line between two genes
    js_args = {
        'exp': bar_dict_source,
        'dynamic': dynamic3,
        'corr_line_title': corr_line.title,
        'corr_bar': corr_bar_sr.data_source,
        'corr_bar_title': corr_bar.title,
        'corr_line_xrange': corr_line.x_range,
        'corr_line_xaxis': corr_line_xaxis,
        'corr_line_yaxis': corr_line_yaxis,
    }
    js_code = """
        var ind = cb_data.index['1d'].indices[0];
        var data = corr_bar.data;
        var gene = data.x[ind];
        var correlation = data.y[ind];
        var gene2 = corr_bar_title.text;

        var d3 = dynamic.data;
        d3['x'] = exp[gene];
        d3['y'] = exp[gene2];
        corr_line_xaxis.axis_label = gene;
        corr_line_yaxis.axis_label = gene2;
        corr_line_title.text = gene + ' vs ' + gene2 + ' = ' + correlation;
        dynamic.change.emit();
    """
    on_hover = CustomJS(args=js_args, code=js_code)

    # define tools
    if corr_test:
        hover = HoverTool(
            tooltips=[
                ('gene', '@x'),
                ('corr', '@y'),
                ('pvalue', '@pval'),
            ],
            callback=on_hover,
        )
    else:
        hover = HoverTool(
            tooltips=[
                ('gene', '@x'),
                ('corr', '@y'),
            ],
            callback=on_hover,
        )
    corr_bar.add_tools(hover)

    # layout
    lout = layout([exp_bar, volcano], [corr_bar, corr_line], sizing_mode='stretch_width')
    save(lout)


def volcano(table, gene_symbol_ind=1, log2fc_ind=2, pvalue_ind=4, exp_ind:list=('6-',),
            fc_cutoff=2.0, pval_cutoff=0.05, limit=5, corr_method='pearson', top=20,
            corr_test=True, corr_pval_cutoff=0.001, label_corr_line=True, link_all_gene=False,
            out_file='volcano_expression.html'):
    volcano_data, bar_data = volcano_source(
        table, gene_symbol_ind=gene_symbol_ind,
        log2fc_ind=log2fc_ind,
        pvalue_ind=pvalue_ind,
        fc_cutoff=fc_cutoff,
        pval_cutoff=pval_cutoff,
        limit=limit,
        exp_ind=exp_ind,
    )
    plot(volcano_data, bar_data,
         out_file=out_file,
         corr_method=corr_method, top=top,
         corr_pval_cutoff=corr_pval_cutoff,
         corr_test=corr_test,
         label_corr_line=label_corr_line,
         only_link_diff_gene=not link_all_gene
         )


def correlation_function(func_name):
    func_name = func_name.lower()
    if func_name == "pearson":
        return stats.pearsonr
    elif func_name == "spearman":
        return stats.spearmanr
    elif func_name == "biserial":
        return stats.pointbiserialr
    elif func_name == 'kendall':
        return stats.kendalltau
    else:
        raise Exception('{} is not in [pearson, spearman, biserial, kendall]')


def df_corr(df, method="pearson", test=False):
    if not test:
        return df.corr(method), None

    corr_func = correlation_function(method)
    coef_matrix = np.zeros((df.shape[1], df.shape[1]))
    pval_matrix = np.zeros((df.shape[1], df.shape[1]))

    for i in range(df.shape[1]):
        v = df[df.columns[i]]
        for j in range(i, df.shape[1]):
            corr, pval = corr_func(v, df[df.columns[j]])
            coef_matrix[i, j] = corr
            coef_matrix[j, i] = corr
            pval_matrix[j, i] = pval
            pval_matrix[i, j] = pval
    df_coef = pd.DataFrame(coef_matrix, columns=df.columns, index=df.columns)
    df_pval = pd.DataFrame(pval_matrix, columns=df.columns, index=df.columns)
    return df_coef, df_pval


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['volcano'])
