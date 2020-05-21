import pygraphviz as pgv
import os


class StateGraph(object):
    def __init__(self, state):
        self.state = state
        self.graph = pgv.AGraph(directed=True, rankdir='TD')

    def _add_nodes(self):
        for node, cmd_info in self.state.items():
            status = cmd_info['state']
            node_detail = node.split('_', 1)
            if status == 'success':
                color = '#7FFF00'
            elif status == 'failed':
                color = '#FFD700'
            elif status == 'running':
                color = '#9F79EE'
            elif status == 'queueing':
                color = '#87CEFF'
            elif status == 'killed':
                color = 'red'
            else:
                color = '#A8A8A8'
            used_time = cmd_info['used_time']
            if isinstance(used_time, str):
                if used_time == 'unknown':
                    pass
                else:
                    try:
                        float(used_time)
                        node_detail.append(used_time+'s')
                    except ValueError:
                        node_detail.append(used_time)
            elif float(used_time) <= 0:
                pass
            else:
                node_detail.append(str(used_time) + 's')
            self.graph.add_node(
                node,
                # 谷歌浏览器可以正常显示tooltip
                tooltip=cmd_info['cmd'].replace(' ', '\n'),
                shape="box",
                style="rounded, filled",
                fillcolor=color,
                color="mediumseagreen",
                label='\n'.join(node_detail)
            )

    def _add_edges(self):
        for target in self.state:
            sources = self.state[target]['depend'].strip()
            if sources:
                sources = sources.split(',')
                edges = zip(sources, [target]*len(sources))
                if self.state[target]['state'] == 'success':
                    color = 'green'
                elif self.state[target]['state'] == 'running':
                    color = '#836FFF'
                else:
                    color = '#4D4D4D'
                self.graph.add_edges_from(edges, color=color)
            else:
                self.graph.add_edge('Input', target, color='green')

    def draw(self, img_file='state.svg'):
        self._add_nodes()
        self._add_edges()
        img_fmt = os.path.splitext(img_file)[1][1:]
        self.graph.draw(path=img_file, format=img_fmt, prog='dot')


state = dict()
fields = ['name', 'state', 'used_time', 'mem', 'cpu', 'pid', 'depend', 'cmd']
value1 = ['去除低质量序列', 'running', '', '', '', '', '', '', '']
value2 = ['对原始数据进行质量统计', 'running', '', '', '', '', '', '', '']
value3 = ['比对到基因组', 'running', '', '', '', '', '去除低质量序列', '']
value4 = ['基于比对结果的质量统计', 'running', '', '', '', '', '比对到基因组', '']
value5 = ['表达定量', 'running', '', '', '', '', '去除低质量序列', '']
value6 = ['差异基因分析', 'running', '', '', '', '', '表达定量', '']
value7 = ['GO富集分析', 'running', '', '', '', '', '差异基因分析', '']
value8 = ['KEGG富集分析', 'running', '', '', '', '', '差异基因分析', '']
for each in [value1, value2, value3, value4, value5, value6, value7, value8]:
    state[each[0]] = dict(zip(fields[1:], each[1:]))
# print(state)
StateGraph(state).draw(img_file='workflow.svg')

