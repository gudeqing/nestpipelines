import pandas as pd
import numpy as np
import math
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
                   fc_cutoff=2.0, pval_cutoff=0.05, limit=3):
    df = pd.read_csv(table, index_col=0, header=0, sep=None, engine='python')
    if gene_symbol_ind >= 2:
        df = df.iloc[:, [gene_symbol_ind-1, log2fc_ind-1, pvalue_ind-1]]
        df.columns = ['gene_symbol', 'log2fc', 'pvalue']
    else:
        df = df.iloc[:, [log2fc_ind-1, pvalue_ind-1]]
        df.columns = ['log2fc', 'pvalue']
        df['gene_symbol'] = df.index
    drop_before = df.shape[0]
    df.loc[:, ['log2fc', 'pvalue']] = df[['log2fc', 'pvalue']].apply(pd.to_numeric, errors='coerce')
    df = df.dropna(axis=0)
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
    return df


def bar_source(table, exp_ind=None):
    df = pd.read_csv(table, index_col=0, header=0, sep=None, engine='python')
    df = df.round(2)
    if exp_ind is not None:
        df = df.iloc[:, exp_ind]
    else:
        df = df.iloc[:, 6:]
    return df


def plot(volcano_df, bar_df, out_file='volcano_expression.html', corr_method='pearson', top=10):
    output_file(out_file)
    if volcano_df.index[0] != volcano_df['gene_symbol'][0]:
        first_gene = volcano_df.index[0] + ' | ' + volcano_df['gene_symbol'][0]
    else:
        first_gene = volcano_df.index[0]
    bar_df = bar_df.loc[volcano_df.index]
    samples = list(bar_df.columns)
    point_groups = set(volcano_df['regulate'])
    point_stat = {x: x+': '+str(list(volcano_df['regulate']).count(x)) for x in point_groups}
    volcano_df['regulate'] = [point_stat[x] for x in volcano_df['regulate']]
    volcano_source = ColumnDataSource(volcano_df)
    gene_list = list(bar_df.index)
    bar_dict_source = bar_df.transpose().to_dict('list')
    # calculate correlation between genes
    corr = bar_df.transpose().corr(method=corr_method)
    corr_dict_source = dict()
    for gene in corr.index:
        tmp_corr = corr.loc[gene].abs().sort_values(ascending=False)
        ind = tmp_corr.index[range(top+1)]
        ind = [x for x in ind if x != gene]
        top_corr = corr.loc[gene].loc[ind]
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
    volcano.legend.background_fill_alpha = 0.5
    volcano.xaxis.axis_label = 'log2(FoldChange)'
    volcano.yaxis.axis_label = '-log10(Pvalue)'

    # corr bar
    corr_bar = figure(**plot_options, x_range=corr_dict_source[first_gene][0])
    dynamic2 = ColumnDataSource({
        'x': corr_dict_source[first_gene][0],
        'y': corr_dict_source[first_gene][1],
    })
    mapper = linear_cmap(palette=RdYlGn[11], field_name='y', low=-1, high=1)
    corr_bar.vbar(x='x', top='y', width=0.5, source=dynamic2, color=mapper)
    corr_bar.xaxis.major_label_orientation = math.pi / 4
    corr_bar.yaxis.axis_label = 'Correlation'
    corr_bar.xgrid.grid_line_color = None
    corr_bar.title.text = first_gene

    # bar plot
    exp_bar = figure(**plot_options, x_range=samples)
    dynamic = ColumnDataSource({
        'x': samples,
        'y': bar_df.loc[first_gene]
    })
    mapper = linear_cmap(palette=RdYlGn[11], field_name='y', low=1, high=1000)
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
        'genes': gene_list,
        'samples': samples,
        'exp_bar_x': exp_bar.x_range,
        'corr': corr_dict_source,
        'dynamic2': dynamic2,
        'corr_bar_x': corr_bar.x_range,
        'title2': corr_bar.title
    }
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
            ('gene', '@gene_symbol'),
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
    corr_hover = HoverTool(
        tooltips=[
            ('gene', '@x'),
            ('corr', '@y'),
        ],
    )
    corr_bar.add_tools(corr_hover)

    # layout
    lout = layout([exp_bar, volcano], corr_bar, sizing_mode='stretch_width')
    save(lout)


def volcano(table, gene_symbol_ind=1, log2fc_ind=2, pvalue_ind=4, exp_ind:list=None,
            fc_cutoff=2.0, pval_cutoff=0.05, limit=3, corr_method='pearson', top=10,
            out_file='volcano_expression.html'):
    volcano_data = volcano_source(table,
                                  gene_symbol_ind=gene_symbol_ind,
                                  log2fc_ind=log2fc_ind,
                                  pvalue_ind=pvalue_ind,
                                  fc_cutoff=fc_cutoff,
                                  pval_cutoff=pval_cutoff,
                                  limit=limit
                                  )
    bar_data = bar_source(table, exp_ind=exp_ind)
    plot(volcano_data, bar_data, out_file=out_file, corr_method=corr_method, top=top)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['volcano'])
