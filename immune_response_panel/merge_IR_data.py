import pandas as pd
import numpy as np


def merge_data(files:list, index_cols:list=None, new_col_name=None, out='merged.table.xls'):
    if index_cols is None:
        index_cols = ['Gene', 'Target', 'ENTREZ_GENE_ID', 'NCBI_NAME', 'NCBI_ACCESSION', 'GENE_FUNCTION']
    table = pd.read_csv(files[0], index_col=None, header=0, sep=None, engine='python')
    table.columns = [x.strip().strip('"') for x in table.columns]
    table.set_index(index_cols, inplace=True)
    for each in files[1:]:
        each_table = pd.read_csv(each, index_col=None, header=0, sep=None, engine='python')
        each_table.columns = [x.strip().strip('"') for x in each_table.columns]
        each_table.set_index(index_cols, inplace=True)
        table = table.join(each_table, how='outer')
    table.columns = [x.strip().strip('"').split('R', 1)[0] for x in table.columns]
    new_name_df = pd.read_csv(new_col_name, index_col=None, header=None, sep=None, engine='python')
    new_name_dict = dict(zip(new_name_df[0], new_name_df[1]))
    new_names = []
    for col in table.columns:
        if col in new_name_dict:
            new_names.append(new_name_dict[col])
        else:
            new_names.append(col)
            print(f'{col} is not found in new name dict')
    table.columns = new_names
    table = table[[x for x in new_name_dict.values() if x in table.columns]]
    print(f'merged table size: {table.shape}')
    table.to_csv(out, sep='\t')


def housekeeping_normalize(count_matrix, housekeeping=None, out='hk.log2.norm.count.xls'):
    if housekeeping is None:
        housekeeping = [
            "GUSB","HMBS","G6PD","POLR2A","ABCF1","TFRC","LRP1","TBP","SDHA","LMNA","TUBB"
        ]
    else:
        housekeeping = [x.strip() for x in open(housekeeping)]
    data = pd.read_csv(count_matrix, index_col=0, header=0, sep=None, engine='python')
    housekeeping_data = data.loc[housekeeping]
    logsum = np.log2(housekeeping_data+1).sum(axis=0)
    factor = logsum/len(housekeeping)
    print('norm_factor:\n', factor)
    normalized = np.log2(data+1) - factor + np.log2(10**6)
    normalized.to_csv(out, sep='\t')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['merge_data', 'housekeeping_normalize'])


