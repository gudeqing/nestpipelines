import matplotlib
matplotlib.use('agg')
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def bar3dplot(data, out='bar3d.png', width=12, height=9):
    data = pd.read_csv(data, header=0, index_col=0, sep=None, engine='python')
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_subplot(111, projection='3d')
    x = list()
    y = list()
    dz = list()
    colors = list()
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            x.append(j)
            y.append(i)
            height = data.iloc[i, j]
            dz.append(height)
            if height < 0.005:
                colors.append('grey')
            elif height < 0.01:
                colors.append('lightgreen')
            elif height < 0.02:
                colors.append('yellow')
            elif height < 0.03:
                colors.append('pink')
            else:
                colors.append('red')
    z = [0]*len(dz)
    dx = [0.5]*len(x)
    dy = [0.5]*len(y)
    ax.bar3d(x, y, z, dx, dy, dz, alpha=0.8, zsort='average', shade=True, color=colors)
    # grid
    # ax.xaxis._axinfo["grid"].update({"linewidth":1, "color" : 'k'})
    ax.set_xticks([x-0.5 for x in range(data.shape[1])])
    ax.set_yticks([x+0.6 for x in range(data.shape[0])])
    ax.set_xticklabels(data.columns, size=6, rotation=90)
    ax.set_yticklabels(data.index, size=7, rotation=90)
    ax.set_zlabel('frequency')
    # plot legend
    lg1 = plt.Rectangle((-1, 0), 1, 1, fc="grey")
    lg2 = plt.Rectangle((-1, 0), 1, 1, fc="lightgreen")
    lg3 = plt.Rectangle((-1, 0), 1, 1, fc="yellow")
    lg4 = plt.Rectangle((-1, 0), 1, 1, fc="pink")
    lg5 = plt.Rectangle((-1, 0), 1, 1, fc="red")
    ax.legend([lg1, lg2, lg3, lg4, lg5], ['<0.005', '<0.01', '<0.02', '<0.03', '>=0.03'])
    # adjust tick label to axis distance
    ax.tick_params(axis='both', which='major', pad=0.1)

    plt.autoscale(enable=True, axis='both', tight=True)
    plt.savefig(out, dpi=300)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['bar3dplot'], log=False)
