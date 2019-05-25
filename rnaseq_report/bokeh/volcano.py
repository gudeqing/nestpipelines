import pandas as pd
import numpy as np
import math
from bokeh.io import output_file
from bokeh.layouts import row
from bokeh.plotting import figure, save
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
    df['color'] = 'grey'
    df.loc[df['regulate']=='Up', 'color'] = 'red'
    df.loc[df['regulate']=='Down', 'color'] = 'green'
    return df


def bar_source(table, exp_ind=None):
    df = pd.read_csv(table, index_col=0, header=0, sep=None, engine='python')
    if exp_ind is not None:
        df = df.iloc[:, exp_ind]
    else:
        df = df.iloc[:, 6:]
    return df


def plot(volcano_source, bar_source, out_file='volcano_expression.html'):
    index_name = bar_source.index.name
    output_file(out_file)
    if volcano_source.index[0] != volcano_source['gene_symbol'][0]:
        first_gene = volcano_source.index[0] + ' | ' + volcano_source['gene_symbol'][0]
    else:
        first_gene = volcano_source.index[0]
    bar_source = bar_source.loc[volcano_source.index]
    samples = list(bar_source.columns)
    point_groups = set(volcano_source['regulate'])
    point_stat = {x: x+': '+str(list(volcano_source['regulate']).count(x)) for x in point_groups}
    volcano_source['regulate'] = [point_stat[x] for x in volcano_source['regulate']]
    volcano_source = ColumnDataSource(volcano_source)
    dynamic = ColumnDataSource({'x': samples, 'y': bar_source.iloc[0]})
    gene_list = list(bar_source.index)
    bar_source = bar_source.transpose().to_dict('list')

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

    # bar plot
    exp_bar = figure(**plot_options, x_range=samples)
    exp_bar.vbar(x='x', top='y', width=0.5, source=dynamic, color='green')
    exp_bar.xaxis.major_label_orientation = math.pi / 4
    exp_bar.yaxis.axis_label = 'Expression'
    # exp_bar.title.text = first_gene

    # interaction
    js_args = {
        'bar': bar_source,
        'dynamic': dynamic,
        'title': exp_bar.title,
        'genes': gene_list
    }
    js_code = """
        var ind = cb_data.index['1d'].indices[0];
        var d3 = dynamic.data;
        console.log(ind);
        var gene = genes[ind];
        title.text = gene;
        d3['y'] = bar[gene];
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
        callback=on_hover
    )
    volcano.add_tools(hover)

    # exp_bar.yaxis.axis_label = first_gene
    bar_hover = HoverTool(
        tooltips=[
            ('sample', '@x'),
        ],
    )
    exp_bar.add_tools(bar_hover)

    # layout
    lout = row(exp_bar, volcano, sizing_mode='stretch_both')
    save(lout)


def volcano(table, gene_symbol_ind=1, log2fc_ind=2, pvalue_ind=4, exp_ind:list=None,
            fc_cutoff=2.0, pval_cutoff=0.05, limit=3, out_file='volcano_expression.html'):
    volcano_data = volcano_source(table,
                                  gene_symbol_ind=gene_symbol_ind,
                                  log2fc_ind=log2fc_ind,
                                  pvalue_ind=pvalue_ind,
                                  fc_cutoff=fc_cutoff,
                                  pval_cutoff=pval_cutoff,
                                  limit=limit
                                  )
    bar_data = bar_source(table, exp_ind=exp_ind)
    plot(volcano_data, bar_data, out_file=out_file)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['volcano'])
