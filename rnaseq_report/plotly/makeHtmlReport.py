import os
import os.path as path
from os.path import join
import re
import yaml
import shutil
from jinja2 import Template
from collections import OrderedDict as dict
import pandas as pd


def make_report_cfg(result_dir, exclude_dirs: list=None, image_formats=('html', 'png', 'xls', 'svg'),
                    out='report.yml'):
    exclude_dirs = exclude_dirs if exclude_dirs else list()
    modules = os.listdir(result_dir)
    modules = [x for x in modules if x not in exclude_dirs and path.isdir(x) and x != 'html.utils']
    if not modules:
        print('No directory found!')
    try:
        modules_order = [(x, int(x.split('.', 1)[0])) for x in modules]
        modules_order.sort(key=lambda x: x[1])
        modules = [x[0] for x in modules_order]
    except:
        print('warn: cannot determine report module order; If order needed, please add "order." before directory name')
    cfg_dict = dict()
    for module in modules:
        slid_dir = join(result_dir, module)
        slides = [x for x in os.listdir(slid_dir) if path.isdir(join(slid_dir, x))]
        try:
            slide_order = [(x, int(x.split('.', 1)[0])) for x in slides]
            slide_order.sort(key=lambda x: x[1])
            slides = [x[0] for x in slide_order]
        except:
            pass
        if not slides:
            print('{} has no sub-diretory and make no slider for this module'.format(module))
            continue
        else:
            cfg_dict[module] = dict()
        for ind, slide in enumerate(slides):
            images = [x for x in os.listdir(join(slid_dir, slide)) if x.endswith(image_formats)]
            if not images:
                print('No image found to make a slide for {}'.format(slide))
                continue
            else:
                cfg_dict[module]['slider_' + str(ind)] = {'name': slide}
            images = [x for x in images if not x.endswith(('xls'))] + [x for x in images if x.endswith(('xls'))]

            # delete existed html for *xls
            possible_html = [x[:-3]+'html' for x in images if x.endswith(('xls'))]
            images = [x for x in images if x not in possible_html]

            cfg_dict[module]['slider_' + str(ind)]['content'] = dict()
            content = cfg_dict[module]['slider_' + str(ind)]['content']
            for ind, img in enumerate(images):
                content['image_' + str(ind+1)] = dict()
                img_info = content['image_' + str(ind+1)]
                desc_file = join(slid_dir, slide, img + '.describe.txt')
                desc = open(desc_file).read().strip() if path.exists(desc_file) else ''
                img_info['name'] = str(ind+1)+ ': ' + str(img.split('.')[0])
                img_info['desc'] = desc.replace('\n', '<br>')
                img_info['path'] = join(slid_dir, slide, img)
                img_info['frmt'] = img.rsplit('.', 1)[1]
                if img_info['frmt'] in ['xls']:
                    # img_info['name'] += "<object>(<a href='{}'>Download Source Table</a>)</object>".format(
                    #     path.relpath(path.abspath(img_info['path']), start=path.abspath(slid_dir))
                    # )
                    use_cols = None
                    if path.exists(img_info['path']+'.describe.txt'):
                        with open(img_info['path']+'.describe.txt') as f:
                            desc = f.readlines()
                        if desc[0].startswith('display_columns'):
                            use_cols = [x.strip() for x in desc[0].split('display_columns:', 1)[1].strip().split(';')]
                            desc = desc[1:]
                        img_info['desc'] = '<br>'.join(desc)
                    table_img = table2html(
                        [img_info['path']],
                        use_cols=use_cols,
                        title=path.basename(img_info['path'])
                    )
                    img_info['path'] = table_img[0]
                    img_info['frmt'] = 'html'

    with open(out, 'w') as f:
        yaml.dump(cfg_dict, f)
    return out


def parse_report_cfg(cfg_file):
    with open(cfg_file, encoding='utf-8') as f:
        cfg_dict = yaml.load(f)
    return cfg_dict


def make_slider(images:list, image_ids:list=None, image_desc:list=None, template="templates/slide.jinja2",
                out='slider.html', link_images=False):
    out_dir = os.path.dirname(out)
    if not out_dir:
        out_dir = os.getcwd()
    else:
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
    if not images:
        return

    if link_images:
        link_dir = os.path.join(out_dir, os.path.basename(out)[:-5]+'.images')
        if not os.path.exists(link_dir):
            os.mkdir(link_dir)
        tmp_list = list()
        for each in images:
            each = os.path.abspath(each)
            target = os.path.join(link_dir, os.path.basename(each))
            if os.path.exists(target):
                continue
            os.symlink(each, target)
            tmp_list.append(target)
        images = tmp_list

    if not image_ids:
        image_ids = [os.path.basename(x).split('.')[0] for x in images]
    if len(image_ids) < len(images):
        raise Exception('number of image ids is less than the number of images')

    if not image_desc:
        image_desc = ['']*len(images)
    if len(image_desc) == 1:
        image_desc = image_desc*len(images)
    else:
        if len(image_desc) < len(images):
            raise Exception('number of image desc is less than the number of images')

    template = Template(open(template).read())
    img_info_dict = dict()
    img_cls = ["slide"]*len(images)
    img_cls[0] = "slide slide--current"
    for img, img_id, desc, tag in zip(images, image_ids, image_desc, img_cls):
        img_info_dict[img_id] = dict()
        img_info_dict[img_id]['path'] = os.path.relpath(os.path.abspath(img), start=os.path.abspath(out_dir))
        img_info_dict[img_id]['cls'] = tag
        img_info_dict[img_id]['description'] = desc
    # print(img_info_dict)
    content = template.render(img_info_dict=img_info_dict)
    with open(out, 'w', encoding='utf-8') as f:
        f.write(content)
    return out


def table_head_div():
    return """
<!DOCTYPE html>
<html lang="en" class="no-js">
<head>
    <title>table</title>
    <meta charset="UTF-8" />
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style type="text/css">
        body{background:white}
        a{text-decoration:none;}
        #myInput {
        background-image: url('https://static.runoob.com/images/mix/searchicon.png');
        background-position: 10px 12px; /* 定位搜索按钮 */
        background-repeat: no-repeat; /* 不重复图片 */
        padding: 12px 12px 12px 40px;
        border: 1px solid #ddd;
        margin-bottom: 1px;
        }

        #myTable {
            border-collapse: collapse;
            width: 100%;
            border: 1px solid #ddd;
            font-size: 10px;
        }

        #myTable th, #myTable td {
            text-align: left;
            padding: 12px;
        }

        #myTable tr {
            /* 表格添加边框 */
            border-bottom: 1px solid #ddd;
        }

        #myTable tr.header, #myTable tr:hover {
            /* 表头及鼠标移动过 tr 时添加背景 */
            background-color: #f1f1f1;
        }
    </style>
    <script type="text/javascript">
        function onSearch(obj){//js函数开始
          setTimeout(function(){//因为是即时查询，需要用setTimeout进行延迟，让值写入到input内，再读取
            var storeId = document.getElementById('myTable');//获取table的id标识
            var rowsLength = storeId.rows.length;//表格总共有多少行
            var key = obj.value;//获取输入框的值
            var searchCol = obj.name;//要搜索的哪一列
            for(var i=1;i<rowsLength;i++){//按表的行数进行循环，本例第一行是标题，所以i=1，从第二行开始筛选（从0数起）
              var searchText = storeId.rows[i].cells[searchCol].innerHTML;//取得table行，列的值
              if(searchText.match(key)){//用match函数进行筛选，如果input的值，即变量 key的值为空，返回的是ture，
                storeId.rows[i].style.display='';//显示行操作，
              }else{
                storeId.rows[i].style.display='none';//隐藏行操作
              }
            }
          },200);//200为延时时间
        }
    </script>
</head>
    """


def table2html(table_file: list, use_cols: list=None, use_rows: list=None, top=500, title='',
                           search_cols:list=None, header:list=None, index_col:list=None, transpose=False):
    pd.set_option('display.max_colwidth', 200)
    results = []
    header = 0 if header is None else [int(x) for x in header]
    index_col = 0 if index_col is None else [int(x) for x in index_col]
    for table in table_file:
        data = pd.read_csv(table, header=header, index_col=index_col, sep='\t')
        if transpose:
            data = data.transpose()
        use_cols = use_cols if use_cols is not None else data.columns
        use_rows = use_rows if use_rows is not None else data.index
        data = data.loc[use_rows, use_cols]
        top = data.shape[0] if use_rows is None else top
        df = data.iloc[:top, :]
        link_dict = dict()
        ind = 0

        def proc(x):
            nonlocal ind
            ind += 1
            if type(x) == float:
                if x < 0.001:
                    return format(x, '.2e')
                else:
                    return round(x, 3)
            else:
                if type(x) == str or type(x) == bytes:
                    if re.match(r'http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+]|[!*(),]|(?:%[0-9a-fA-F][0-9a-fA-F]))+', x):
                        link_key = 'href_link_{}'.format(ind)
                        link_dict[link_key] = x
                        return link_key
                    else:
                        return x
                else:
                    return x

        df = df.applymap(proc)
        table_div = df.to_html(buf=None, header=True, index=True, na_rep='NaN',
                               formatters=None, float_format=None, sparsify=None, index_names=False, justify=None,
                               bold_rows=True, classes=None, escape=True, max_rows=None, max_cols=None,
                               show_dimensions=True, notebook=False, decimal='.', border=None)

        table_div = table_div.replace('table', 'table id="myTable"', 1)
        for k, v in link_dict.items():
            table_div = table_div.replace(k, '<a href="{}" target=_blank>link</a>'.format(v), 1)
        final_html = table_head_div() + '\n<body>'
        search_cols = [0] if search_cols is None else search_cols
        for each in search_cols:
            final_html += '\n<input name="{}" type="text" id="myInput" onkeydown="onSearch(this)" ' \
                          'value="" placeholder="搜素第{}列"/>'.format(each, int(each)+1)
        if top < data.shape[0]:
            final_html += '&nbsp;Only show the first {} rows'.format(top)
        if title:
            final_html += '&nbsp;<a href="{}" target=_blank>{}</a>'.format(title, path.basename(title))
        final_html += '\n'+ table_div + '\n</body>'
        out_name = table.rsplit('.', 1)[0] + '.html'
        with open(out_name, 'w', encoding='utf-8') as f:
            f.write(final_html)
        results.append(out_name)
    return results


def make_report(cfg_from, report_dir=None, link_images=False, exclude_dirs:list=None,
                image_formats:tuple=('html', 'png', 'xls', 'svg')):
    """
    make html slider
    :param cfg_from: 结果目录或生成报告的配置文件, 如果提供结果目录, 则自动从结果目录生成Report configuration文件
    :param report_dir: 报告输出目录，默认为当前目录
    :param link_images: 若设置，则把图片软链接到报告目录，默认不做软链接
    :param exclude_dirs: 排除指定目录，使其不参与报告制作，当cfg_from为结果目录时，该参数有效
    :param image_formats: 指定可以用于制作报告的图片格式，后续将根据文件后缀是否与之匹配决定提取的图片，当cfg_from为结果目录时，该参数有效
    :return:
    """
    report_dir = os.getcwd() if report_dir is None else report_dir
    if path.isdir(cfg_from):
        cfg_file = make_report_cfg(cfg_from, image_formats=image_formats, exclude_dirs=exclude_dirs,
                                   out=join(report_dir, 'report.yml'))
    else:
        cfg_file = cfg_from
    cfg_dict = parse_report_cfg(cfg_file)
    index_html_path = join(report_dir, 'index.html')
    html_utils_dir = os.path.dirname(os.path.dirname(__file__))
    index_template = join(html_utils_dir, 'templates', 'index.jinja2')
    slide_template = join(html_utils_dir, 'templates', 'slide.jinja2')
    report_dict = dict()
    for module in cfg_dict:
        report_dict[module] = dict()
        for slider in cfg_dict[module]:
            slider = cfg_dict[module][slider]
            img_name = []
            img_path = []
            img_desc = []
            img_frmt = []
            out_name = join(report_dir, module, slider['name']+'.html')
            for image in slider['content']:
                image = slider['content'][image]
                img_name.append(image['name'])
                img_path.append(image['path'])
                img_desc.append(image['desc'])
                img_frmt.append(image['frmt'])
            make_slider(images=img_path, image_ids=img_name, image_desc=img_desc,
                        template=slide_template, link_images=link_images, out=out_name)
            report_dict[module][slider['name']] = out_name
    # format report use relative path, this is important
    for section in report_dict:
        for each in report_dict[section]:
            # print(section)
            tmp_path = path.abspath(report_dict[section][each])
            report_dict[section][each] = path.relpath(tmp_path, start=path.abspath(report_dir))
    # format html
    if not path.exists(join(report_dir, 'html.utils')):
        shutil.copytree(join(html_utils_dir, 'html.utils'), join(report_dir, 'html.utils'))
    template = Template(open(index_template).read())
    content = template.render(result=report_dict)
    with open(index_html_path, 'w', encoding='utf-8') as f:
        f.write(content)


if __name__ == '__main__':
    class Func2Command(object):
        def __init__(self, callable_dict):
            self.call(callable_dict)

        @staticmethod
        def introduce_command(func):
            import argparse
            import inspect
            import json
            import time
            import sys
            if isinstance(func, type):
                description = func.__init__.__doc__
            else:
                description = func.__doc__
            if '-h' not in sys.argv or '--help' in sys.argv or '-help' in sys.argv:
                description = None
            if description:
                _ = [print(x.strip()) for x in description.split('\n') if x.strip()]
                parser = argparse.ArgumentParser(add_help=False)
            else:
                parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=description)
            func_args = inspect.getfullargspec(func)
            arg_names = func_args.args
            if not arg_names:
                func()
                return
            arg_defaults = func_args.defaults
            if not arg_defaults:
                arg_defaults = []
            arg_defaults = ['None']*(len(arg_names) - len(arg_defaults)) + list(arg_defaults)
            sig = inspect.signature(func)
            for arg, value in zip(arg_names, arg_defaults):
                if arg == 'self':
                    continue
                arg_type = sig.parameters[arg].annotation
                if value == 'None':
                    if arg_type in [list, tuple, set]:
                        parser.add_argument('-' + arg, nargs='+', required=True, metavar=arg)
                    elif arg_type in [int, str, float]:
                        parser.add_argument('-' + arg, type=arg_type, required=True, metavar=arg)
                    else:
                        parser.add_argument('-'+arg, required=True, metavar=arg)
                elif type(value) == bool:
                    if value:
                        parser.add_argument('--'+arg, action="store_false", help='default: True')
                    else:
                        parser.add_argument('--'+arg, action="store_true", help='default: False')
                elif value is None:
                    if arg_type in [list, tuple, set]:
                        parser.add_argument('-' + arg, nargs='+', required=False, metavar=arg)
                    elif arg_type in [int, str, float]:
                        parser.add_argument('-' + arg, type=arg_type, required=False, metavar=arg)
                    else:
                        parser.add_argument('-'+arg, required=False, metavar='Default:' + str(value))
                else:
                    if arg_type in [list, tuple, set] or (type(value) in [list, tuple, set]):
                        default_value = ' '.join(str(x) for x in value)
                        if type(value) in [list, tuple]:
                            one_value = value[0]
                        else:
                            one_value = value.pop()
                        parser.add_argument('-' + arg, default=value, nargs='+', type=type(one_value),
                                            metavar='Default:'+default_value, )
                    else:
                        parser.add_argument('-' + arg, default=value, type=type(value),
                                            metavar='Default:' + str(value), )
            if func_args.varargs is not None:
                print("warning: *varargs is not supported, and will be neglected! ")
            if func_args.varkw is not None:
                print("warning: **keywords args is not supported, and will be neglected! ")
            args = parser.parse_args().__dict__
            # try:
            #     with open("Argument_detail.json", 'w') as f:
            #         json.dump(args, f, indent=2, sort_keys=True)
            # except IOError:
            #     print('Current Directory is not writable, thus argument log is not written !')
            start = time.time()
            func(**args)
            print("total time: {}s".format(time.time() - start))

        def call(self, callable_dict):
            import sys
            excludes = ['introduce_command', 'Func2Command']
            _ = [callable_dict.pop(x) for x in excludes if x in callable_dict]
            if len(callable_dict) >= 2:
                if len(sys.argv) <= 1:
                    print("The tool has the following sub-commands: ")
                    _ = [print(x) for x in callable_dict]
                    return
                sub_cmd = sys.argv[1]
                sys.argv.remove(sub_cmd)

                if sub_cmd in callable_dict:
                    self.introduce_command(callable_dict[sub_cmd])
                else:
                    print('sub-command: {} is not defined'.format(sub_cmd))
            elif len(callable_dict) == 1:
                self.introduce_command(callable_dict.pop(list(callable_dict.keys())[0]))

    callable_dict = {x: y for x, y in locals().items() if callable(y) and len(x)>8}
    _ = [callable_dict.pop(x) for x in {'Func2Command'} if x in callable_dict]
    Func2Command(callable_dict)
