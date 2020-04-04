import pandas as pd


def de50(files:list):
    for each in files:
        a = pd.read_csv(each, header=0)
        cols = ['Total Counts', 'Frequency']
        a = a[cols]
        sf = 0
        st = 0
        total = a.shape[0]
        for i in range(a.shape[0]):
            sf += a.iloc[i, 1]
            st += 1
            if sf >= 0.5:
                #print(each)
                #print(st, total)
                print(each, st/total)
                break


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['de50'])

