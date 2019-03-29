import pygraphviz as pgv

h1 = '启动程序，收集输入参数'
h2 = '有无参数输入'
h31 = '模拟数据并生成表达矩阵热图\n程序终止'
h2_h31 = '否'
h2_h32 = '是'
h32 = '根据参数内容进行分析决策'
h41 = "是否对基因聚类"
h42 = "是否对样本聚类"
h43 = "是否仅做样本相关性聚类\n并展示相关性热图"
h44 = "是否仅展示样本聚类树\n不画热图"
h45 = "是否仅展示基因聚类树\n不画热图"
h5 = "收集决策"
h61 = "数据预处理"
h62 = "设定画图布局"
h7 = "分析和可视化"
h8 = "输出结果，程序终止"

G = pgv.AGraph(strict=True, directed=True, rankdir='TB')
G.add_node(h1, shape='octagon')
G.add_node(h2, shape='diamond')
G.add_edge(h1, h2)
G.add_node(h31, shape='rectangle')
G.add_edge(h2, h31, label=h2_h31)
G.add_node(h32, shape='rectangle')
G.add_edge(h2, h32, label=h2_h32)
G.add_node(h41, shape='rectangle')
G.add_node(h42, shape='rectangle')
G.add_node(h43, shape='rectangle')
G.add_node(h44, shape='rectangle')
G.add_node(h45, shape='rectangle')
G.add_edge(h32, h41)
G.add_edge(h32, h42)
G.add_edge(h32, h43)
G.add_edge(h32, h44)
G.add_edge(h32, h45)
G.add_edge(h41, h5)
G.add_edge(h42, h5)
G.add_edge(h43, h5)
G.add_edge(h44, h5)
G.add_edge(h45, h5)
G.add_node(h5, shape='rectangle')
G.add_node(h61, shape='rectangle')
G.add_node(h62, shape='rectangle')
G.add_node(h7, shape='rectangle')
G.add_node(h8, shape='oval')
G.add_edge(h5, h61)
G.add_edge(h5, h62)
G.add_edge(h61, h7)
G.add_edge(h62, h7)
G.add_edge(h7, h8)

img_file = 'workflow.svg'
img_fmt = 'svg'
G.draw(path=img_file, format=img_fmt, prog='dot')



